# =============================================================================
# Greedy active subnetwork search.
#
# Design notes on mutability and performance
# ------------------------------------------
# R passes all vectors by value. To share mutable state across recursive calls
# we split state into two categories:
#
#   vstate (environment)  – the four large per-node vectors that must be visible
#                           to every recursive frame simultaneously:
#                             $comp_members       logical[N]
#                             $removable_vec      logical[N]
#                             $node2predecessor   integer[N]
#                             $node2dep_count     integer[N]
#
#   scalars (function args / return values)  – size, zsum, score, best_score.
#     These are passed down as arguments and the only one that needs to bubble
#     UP the call stack is best_score, which is returned.  This avoids `env$`
#     read/write overhead on the four hot scalars (~2x faster than storing them
#     in the environment alongside the vectors).
#
# =============================================================================

#' Build the within-depth reachability vector via iterative BFS
#'
#' Returns a logical vector rather than modifying in-place, because \code{<<-}
#' resolves through *lexical* parent environments, not the call stack, and
#' therefore does not reach a variable in the calling function's local frame.
#'
#' @param nbr_idx  List of integer neighbour-index vectors.
#' @param n_nodes  Total node count.
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
  dist <- integer(n_nodes)

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
#' The four hot scalars (size, zsum, score, best_score) are passed as function
#' arguments and returned, avoiding \code{env$} overhead on every add/remove.
#' The four large vectors live in \code{vstate} (an environment) so mutations
#' are immediately visible to all recursive frames.
#'
#' \code{.score_subnetwork} is inlined: single-node → 0, else calibrated.
#'
#' @param nbr_idx      List of integer neighbour-index vectors.
#' @param z_vec        Numeric z-score vector (integer-indexed).
#' @param sc_means     Monte-Carlo means (indexed by subnetwork size).
#' @param sc_stds      Monte-Carlo stds  (indexed by subnetwork size).
#' @param within_vec   Logical reachability mask, or NULL for no limit.
#' @param search_depth Depth budget restored on every improvement.
#' @param depth        Remaining depth budget for this call.
#' @param size         Current component size.
#' @param zsum         Current z-score sum.
#' @param score        Current calibrated score.
#' @param best_score   Best score seen anywhere in this subtree so far.
#' @param vstate       Environment holding the four mutable node vectors.
#' @param last_added   Integer node ID most recently added.
#'
#' @return A list with \code{$improved} (logical) and \code{$best_score}
#'   (the highest score seen in this subtree).
.greedy_expand_idx <- function(nbr_idx, z_vec, sc_means, sc_stds,
                               within_vec, search_depth, depth,
                               size, zsum, score, best_score,
                               vstate, last_added) {
  improved <- FALSE

  if (score > best_score) {
    depth <- search_depth
    improved <- TRUE
    best_score <- score
  }

  if (depth > 0L) {
    any_improved <- FALSE
    vstate$removable_vec[last_added] <- FALSE
    dep_count <- 0L

    for (nb in nbr_idx[[last_added]]) {
      within_ok <- is.null(within_vec) || within_vec[nb]
      if (within_ok && !vstate$comp_members[nb]) {
        # --- add nb ---
        new_size <- size + 1L
        new_zsum <- zsum + z_vec[nb]
        new_score <- (new_zsum / sqrt(new_size) - sc_means[new_size]) / sc_stds[new_size]

        vstate$comp_members[nb] <- TRUE
        vstate$removable_vec[nb] <- TRUE

        res <- .greedy_expand_idx(
          nbr_idx, z_vec, sc_means, sc_stds,
          within_vec, search_depth, depth - 1L,
          new_size, new_zsum, new_score, best_score,
          vstate, nb
        )
        best_score <- res$best_score

        if (!res$improved) {
          # backtrack
          vstate$comp_members[nb] <- FALSE
          vstate$removable_vec[nb] <- FALSE
        } else {
          dep_count <- dep_count + 1L
          any_improved <- TRUE
          vstate$node2predecessor[nb] <- last_added
        }
      }
    }

    improved <- improved || any_improved
    if (dep_count > 0L) {
      vstate$removable_vec[last_added] <- FALSE
      vstate$node2dep_count[last_added] <- dep_count
    }
  }

  list(improved = improved, best_score = best_score)
}

#' Greedy removal pass
#'
#' Tries removing each removable leaf; keeps the removal if the score
#' improves. All state lives in \code{vstate}; the four scalars are passed
#' and the updated best_score is returned.
#'
#' @param vstate   Environment holding mutable node vectors.
#' @param z_vec    Numeric z-score vector.
#' @param sc_means Monte-Carlo means.
#' @param sc_stds  Monte-Carlo stds.
#' @param size     Current component size.
#' @param zsum     Current z-score sum.
#' @param score    Current calibrated score.
#' @param best_score Best score seen so far.
#'
#' @return Updated best_score (numeric).
.greedy_removal_idx <- function(vstate, z_vec, sc_means, sc_stds,
                                size, zsum, score, best_score) {
  snapshot <- which(vstate$removable_vec)
  for (cur in snapshot) {
    new_size <- size - 1L
    new_zsum <- zsum - z_vec[cur]
    new_score <- if (new_size <= 1L) {
      0
    } else {
      (new_zsum / sqrt(new_size) - sc_means[new_size]) / sc_stds[new_size]
    }

    if (new_score > best_score) {
      best_score <- new_score
      size <- new_size
      zsum <- new_zsum
      score <- new_score
      vstate$comp_members[cur] <- FALSE
      pred <- vstate$node2predecessor[cur]
      if (pred > 0L) {
        dep <- vstate$node2dep_count[pred] - 1L
        if (dep == 0L) {
          vstate$removable_vec[pred] <- TRUE
        } else {
          vstate$node2dep_count[pred] <- dep
        }
      }
    }
    # else: leave cur in the component (no change needed)
  }
  best_score
}

#' Filter subnetworks by mutual overlap
#'
#' @param comps            Score-sorted list of subnetwork objects (desc).
#' @param overlap_threshold Numeric threshold in (0, 1).
#' @param subnetwork_num   Maximum number of subnetworks to return.
#'
#' @return Filtered list.
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
  set.seed(params$seed)
  nodes <- network$nodes
  n_nodes <- length(nodes)
  if (n_nodes == 0L) {
    return(list())
  }

  max_depth <- params$gr_max_depth
  search_depth <- params$gr_search_depth

  z_vec <- as.numeric(sc$z[nodes])
  sc_means <- sc$means
  sc_stds <- sc$stds

  # Build neighbour index from the integer edge list — avoids repeated
  # string-keyed hash lookups
  el <- igraph::as_edgelist(network$g, names = FALSE)
  both <- rbind(el, el[, 2:1, drop = FALSE])
  both <- both[order(both[, 1L]), , drop = FALSE]
  row_starts <- c(1L, which(diff(both[, 1L]) > 0L) + 1L, nrow(both) + 1L)

  nbr_idx <- vector("list", n_nodes)
  for (i in seq_len(n_nodes)) {
    s <- row_starts[i]
    e <- row_starts[i + 1L] - 1L
    nb <- if (s <= e) both[s:e, 2L] else integer(0L)
    nbr_idx[[i]] <- nb[z_vec[nb] > -1.0] # prune strongly negative neighbours
  }

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

  # Pre-allocated zero-filled templates for per-seed vector resets
  lv_template <- logical(n_nodes)
  iv_template <- integer(n_nodes)

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

    # vstate: environment holding only the four per-node vectors.
    # Scalars (size, zsum, score, best_score) are kept as local variables.
    vstate <- new.env(parent = emptyenv())
    vstate$comp_members <- lv_template # template copy, not fresh alloc
    vstate$removable_vec <- lv_template
    vstate$node2predecessor <- iv_template
    vstate$node2dep_count <- iv_template

    vstate$comp_members[seed_id] <- TRUE
    vstate$node2dep_count[seed_id] <- 1L

    seed_zsum <- z_vec[seed_id]
    seed_score <- 0 # single-node always scores 0

    res <- .greedy_expand_idx(
      nbr_idx, z_vec, sc_means, sc_stds,
      within_vec, search_depth, search_depth,
      1L, seed_zsum, seed_score, -Inf,
      vstate, seed_id
    )

    final_best <- .greedy_removal_idx(
      vstate, z_vec, sc_means, sc_stds,
      # size/zsum/score after expand: reconstruct from comp_members
      # (we don't track them through expand's return value for removal —
      #  removal operates on whatever is in vstate$comp_members)
      sum(vstate$comp_members),
      sum(z_vec[vstate$comp_members]),
      0, # score recomputed inside removal
      res$best_score
    )

    final_idx <- which(vstate$comp_members)

    if (length(final_idx) >= 2L && final_best > 0) {
      comp_obj <- list(nodes = nodes[final_idx], score = final_best)
      for (ni in final_idx) {
        if (best_score_per_node[ni] < final_best) {
          best_score_per_node[ni] <- final_best
          best_comp_per_node[[ni]] <- comp_obj
        }
      }
    }
  }
  if (verbose) message("100%")

  # Collect, deduplicate (by sorted integer index), sort by score
  comps <- best_comp_per_node[!vapply(best_comp_per_node, is.null, logical(1))]
  if (length(comps) == 0L) {
    return(list())
  }

  sig <- vapply(
    comps,
    function(c) paste0(sort.int(match(c$nodes, nodes)), collapse = "\x01"),
    character(1)
  )
  comps <- comps[!duplicated(sig)]

  scores <- vapply(comps, function(c) c$score, numeric(1))
  comps <- comps[order(scores, decreasing = TRUE)]

  .greedy_filter(comps, params$gr_overlap_threshold, params$gr_subnetwork_num)
}
