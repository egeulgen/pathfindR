# =============================================================================
# Greedy active subnetwork search.
#
# Design notes on mutability
# --------------------------
# R passes all vectors (logical, integer, numeric) by value. Any function that
# needs to mutate shared state across recursive calls must store that state in
# an *environment*, which is the only reference-semantic object in base R.
#
# The expand / removal routines therefore receive a single `state` environment
# that holds every mutable field:
#   state$size              – current component size (integer)
#   state$zsum              – current z-score sum (double)
#   state$score             – current calibrated score (double)
#   state$best_score        – best score seen so far (double)
#   state$comp_members      – logical[n_nodes], TRUE = in component
#   state$removable_vec     – logical[n_nodes], TRUE = leaf candidate
#   state$node2predecessor  – integer[n_nodes], predecessor on the DFS path
#   state$node2dep_count    – integer[n_nodes], number of dependent children
# =============================================================================

#' Build the within-depth reachability vector via iterative BFS
#'
#' Returns a logical vector (not modifies in place) so that the result is
#' usable from any calling context without relying on \code{<<-}, which
#' resolves to *lexical* parent environments rather than the calling frame.
#'
#' @param nbr_idx  List of integer neighbour-index vectors.
#' @param n_nodes  Total node count (length of returned vector).
#' @param start    Integer node ID to expand from.
#' @param depth    Maximum hop distance.
#'
#' @return Logical vector of length \code{n_nodes}.
.greedy_init_max_depth_idx <- function(nbr_idx, n_nodes, start, depth) {
  within_vec <- logical(n_nodes)
  within_vec[start] <- TRUE
  if (depth == 0L) {
    return(within_vec)
  }

  frontier <- start
  dist <- integer(n_nodes) # 0 everywhere; start treated as distance 0

  while (length(frontier) > 0L) {
    next_frontier <- integer(0L)
    for (cur in frontier) {
      d_next <- dist[cur] + 1L
      if (d_next <= depth) {
        for (nb in nbr_idx[[cur]]) {
          if (!within_vec[nb]) {
            within_vec[nb] <- TRUE
            dist[nb] <- d_next
            next_frontier <- c(next_frontier, nb)
          }
        }
      }
    }
    frontier <- next_frontier
  }
  within_vec
}

#' Recursive greedy expansion
#'
#' All mutable per-seed state lives in the \code{state} environment so that
#' changes made in recursive calls are immediately visible to the parent call.
#' \code{.score_subnetwork} is inlined here (~28x faster than a function call
#' in a tight loop, as measured by micro-benchmark).
#'
#' @param nbr_idx      List of integer neighbour-index vectors.
#' @param z_vec        Numeric z-score vector aligned to node indices.
#' @param sc_means     Numeric vector of MC means (indexed by subnetwork size).
#' @param sc_stds      Numeric vector of MC stds  (indexed by subnetwork size).
#' @param within_vec   Logical vector of allowed nodes, or NULL for no limit.
#' @param search_depth Depth budget restored on every improvement.
#' @param depth        Remaining depth budget for this call.
#' @param state        Environment holding all mutable component state.
#' @param last_added   Integer node ID most recently added.
#'
#' @return Logical; TRUE if the best score improved anywhere in this subtree.
.greedy_expand_idx <- function(nbr_idx, z_vec, sc_means, sc_stds,
                               within_vec, search_depth, depth,
                               state, last_added) {
  improved <- FALSE

  if (state$score > state$best_score) {
    depth <- search_depth
    improved <- TRUE
    state$best_score <- state$score
  }

  if (depth > 0L) {
    any_improved <- FALSE
    state$removable_vec[last_added] <- FALSE
    dependent_count <- 0L

    for (nb in nbr_idx[[last_added]]) {
      within_ok <- is.null(within_vec) || within_vec[nb]
      if (within_ok && !state$comp_members[nb]) {
        # --- add nb ---
        state$size <- state$size + 1L
        state$comp_members[nb] <- TRUE
        state$zsum <- state$zsum + z_vec[nb]
        n <- state$size
        # inlined .score_subnetwork (single-node → 0, else calibrated)
        state$score <- if (n == 1L) {
          0
        } else {
          (state$zsum / sqrt(n) - sc_means[n]) / sc_stds[n]
        }
        state$removable_vec[nb] <- TRUE

        this_improved <- .greedy_expand_idx(
          nbr_idx, z_vec, sc_means, sc_stds,
          within_vec, search_depth, depth - 1L,
          state, nb
        )

        if (!this_improved) {
          # --- remove nb (backtrack) ---
          state$size <- state$size - 1L
          state$comp_members[nb] <- FALSE
          state$zsum <- state$zsum - z_vec[nb]
          n <- state$size
          state$score <- if (n == 1L) {
            0
          } else {
            (state$zsum / sqrt(n) - sc_means[n]) / sc_stds[n]
          }
          state$removable_vec[nb] <- FALSE
        } else {
          dependent_count <- dependent_count + 1L
          any_improved <- TRUE
          state$node2predecessor[nb] <- last_added
        }
      }
    }

    improved <- improved || any_improved
    if (dependent_count > 0L) {
      state$removable_vec[last_added] <- FALSE
      state$node2dep_count[last_added] <- dependent_count
    }
  }

  improved
}

#' Greedy removal pass
#'
#' Tries removing each leaf node; keeps the removal if the score improves.
#' All state is mutated in place through the \code{state} environment.
#'
#' @param state    Environment holding all mutable component state.
#' @param z_vec    Numeric z-score vector aligned to node indices.
#' @param sc_means Numeric vector of MC means.
#' @param sc_stds  Numeric vector of MC stds.
#'
#' @return Invisibly NULL.
.greedy_removal_idx <- function(state, z_vec, sc_means, sc_stds) {
  snapshot <- which(state$removable_vec)
  for (cur in snapshot) {
    old_size <- state$size
    old_zsum <- state$zsum
    old_score <- state$score

    state$size <- old_size - 1L
    state$comp_members[cur] <- FALSE
    state$zsum <- old_zsum - z_vec[cur]
    n <- state$size
    state$score <- if (n == 1L) {
      0
    } else {
      (state$zsum / sqrt(n) - sc_means[n]) / sc_stds[n]
    }

    if (state$score > state$best_score) {
      state$best_score <- state$score
      pred <- state$node2predecessor[cur]
      if (pred > 0L) {
        dep <- state$node2dep_count[pred] - 1L
        if (dep == 0L) {
          state$removable_vec[pred] <- TRUE
        } else {
          state$node2dep_count[pred] <- dep
        }
      }
    } else {
      # put back
      state$size <- old_size
      state$comp_members[cur] <- TRUE
      state$zsum <- old_zsum
      state$score <- old_score
    }
  }
  invisible(NULL)
}

#' Filter subnetworks by mutual overlap
#'
#' Walks the score-sorted list, keeping a subnetwork and discarding any later
#' one whose Jaccard overlap with a kept one exceeds the threshold.
#'
#' @param comps           List of subnetwork objects (score-sorted, descending).
#' @param overlap_threshold Numeric threshold in (0, 1).
#' @param subnetwork_num  Maximum number of subnetworks to return.
#'
#' @return Filtered list of subnetwork objects.
.greedy_filter <- function(comps, overlap_threshold, subnetwork_num) {
  n <- length(comps)
  if (n == 0L) {
    return(list())
  }

  node_lists <- lapply(comps, function(s) s$nodes)
  all_nodes <- unique(unlist(node_lists, use.names = FALSE))
  id_map <- stats::setNames(seq_along(all_nodes), all_nodes)
  ids <- lapply(node_lists, function(ns) unname(id_map[ns]))
  sizes <- lengths(ids)

  mark <- logical(length(all_nodes))
  to_delete <- logical(n)
  filtered <- vector("list", 0L)

  i <- 1L
  while (i <= n && length(filtered) < subnetwork_num) {
    if (!to_delete[i]) {
      s1_ids <- ids[[i]]
      s1_size <- sizes[i]
      filtered[[length(filtered) + 1L]] <- comps[[i]]

      mark[s1_ids] <- TRUE
      if (i < n) {
        for (j in seq.int(i + 1L, n)) {
          if (!to_delete[j]) {
            common <- sum(mark[ids[[j]]])
            denom <- if (s1_size < sizes[j]) s1_size else sizes[j]
            if (common / denom > overlap_threshold) to_delete[j] <- TRUE
          }
        }
      }
      mark[s1_ids] <- FALSE
    }
    i <- i + 1L
  }

  filtered
}

#' Run the greedy active subnetwork search
#'
#' @param network  A network from \code{.build_network()}.
#' @param sc       A score context from \code{.build_score_context()}.
#' @param params   A list of run parameters.
#' @param verbose  Logical; emit progress messages.
#'
#' @return A list of subnetwork objects.
.greedy_search <- function(network, sc, params, verbose = FALSE) {
  nodes <- network$nodes
  n_nodes <- length(nodes)
  if (n_nodes == 0L) {
    return(list())
  }

  max_depth <- params$gr_max_depth
  search_depth <- params$gr_search_depth

  # Aligned z-score vector (integer-indexed, avoids name lookup)
  z_vec <- as.numeric(sc$z[nodes])
  sc_means <- sc$means
  sc_stds <- sc$stds

  # Build neighbour index (integer IDs, pre-filtered by z-score)
  node_map <- stats::setNames(seq_len(n_nodes), nodes)
  nbr_idx <- lapply(seq_len(n_nodes), function(i) {
    nbs <- unname(node_map[network$nbr[[nodes[i]]]])
    nbs[z_vec[nbs] > -1.0] # prune strongly negative neighbours
  })

  # Seed selection: top-5% z-score nodes (computed once)
  z_cutoff <- stats::quantile(z_vec, probs = 0.95, na.rm = TRUE)
  promising_seeds <- which(z_vec >= z_cutoff)

  if (length(promising_seeds) > n_nodes * 0.20) {
    promising_seeds <- order(z_vec, decreasing = TRUE)[
      seq_len(max(1L, as.integer(n_nodes * 0.05)))
    ]
  } else if (length(promising_seeds) < 50L && n_nodes >= 50L) {
    promising_seeds <- order(z_vec, decreasing = TRUE)[seq_len(50L)]
  }

  n_seeds <- length(promising_seeds)
  best_score_per_node <- rep(-Inf, n_nodes)
  best_comp_per_node <- vector("list", n_nodes)

  percent <- 0L
  for (i in seq_along(promising_seeds)) {
    seed_id <- promising_seeds[i]

    if (verbose) {
      new_pct <- (100L * (i - 1L)) %/% n_seeds
      if (new_pct > percent) {
        percent <- new_pct
        message(percent, "% (", i, "/", n_seeds, " seeds)")
      }
    }

    within_vec <- if (max_depth == 0L) {
      NULL
    } else {
      .greedy_init_max_depth_idx(nbr_idx, n_nodes, seed_id, max_depth)
    }

    # All mutable per-seed state lives in one environment so recursive
    # calls see each other's mutations immediately.
    state <- new.env(parent = emptyenv())
    state$size <- 1L
    state$zsum <- z_vec[seed_id]
    state$score <- 0 # single-node always scores 0
    state$best_score <- -Inf
    state$comp_members <- logical(n_nodes)
    state$removable_vec <- logical(n_nodes)
    state$node2predecessor <- integer(n_nodes)
    state$node2dep_count <- integer(n_nodes)

    state$comp_members[seed_id] <- TRUE
    state$node2dep_count[seed_id] <- 1L

    .greedy_expand_idx(
      nbr_idx, z_vec, sc_means, sc_stds,
      within_vec, search_depth, search_depth,
      state, seed_id
    )

    .greedy_removal_idx(state, z_vec, sc_means, sc_stds)

    final_idx <- which(state$comp_members)
    final_score <- state$best_score

    if (length(final_idx) >= 2L && final_score > 0) {
      comp_obj <- list(nodes = nodes[final_idx], score = final_score)
      for (ni in final_idx) {
        if (best_score_per_node[ni] < final_score) {
          best_score_per_node[ni] <- final_score
          best_comp_per_node[[ni]] <- comp_obj
        }
      }
    }
  }
  if (verbose) message("100%")

  # Collect, deduplicate, sort
  comps <- best_comp_per_node[!vapply(best_comp_per_node, is.null, logical(1))]
  if (length(comps) == 0L) {
    return(list())
  }

  sig <- vapply(
    comps,
    function(c) paste(sort.int(match(c$nodes, nodes)), collapse = " "),
    character(1)
  )
  comps <- comps[!duplicated(sig)]

  scores <- vapply(comps, function(c) c$score, numeric(1))
  comps <- comps[order(scores, decreasing = TRUE)]

  .greedy_filter(comps, params$gr_overlap_threshold, params$gr_subnetwork_num)
}
