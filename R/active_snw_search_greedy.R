# =============================================================================
# Greedy active subnetwork search. Each node is used as a seed and the current
# subnetwork is recursively expanded through neighbouring nodes whenever the
# calibrated score improves. A subsequent pruning step removes dispensable
# nodes, after which duplicate and highly overlapping subnetworks are filtered
# to yield the final ranked candidates.
# =============================================================================

#' Build the within-depth reachability vector via iterative BFS
#'
#' @param nbr_idx  List of integer neighbour-index vectors.
#' @param n_nodes  Total node count.
#' @param start    Integer seed node ID.
#' @param depth    Maximum hop distance.
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
    nf <- integer(0L)
    for (cur in frontier) {
      d_next <- dist[cur] + 1L
      if (d_next <= depth) {
        for (nb in nbr_idx[[cur]]) {
          if (!within_vec[nb]) {
            within_vec[nb] <- TRUE
            dist[nb] <- d_next
            nf <- c(nf, nb)
          }
        }
      }
    }
    frontier <- nf
  }
  within_vec
}

#' Recursive greedy expansion
#'
#' Scalars are passed as arguments and returned to avoid env$ lookup overhead.
#' Returns size and zsum alongside improved/best_score so the caller can pass
#' them directly to the removal pass (avoids two O(N) sum() scans per seed).
#'
#' @param nbr_idx      Integer neighbour-index list.
#' @param z_vec        Numeric z-score vector.
#' @param sc_means     Monte-Carlo means (by subnetwork size).
#' @param sc_stds      Monte-Carlo stds  (by subnetwork size).
#' @param within_vec   Logical reachability mask or NULL.
#' @param search_depth Depth budget restored on improvement.
#' @param depth        Remaining depth budget.
#' @param size         Current component size.
#' @param zsum         Current z-score sum.
#' @param score        Current calibrated score.
#' @param best_score   Best score seen in this subtree.
#' @param vstate       Environment with four mutable node vectors.
#' @param last_added   Integer node ID most recently added.
#' @return List: improved, best_score, size, zsum.
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
          vstate$comp_members[nb] <- FALSE
          vstate$removable_vec[nb] <- FALSE
        } else {
          # Accept: update size/zsum to reflect the child's final state
          size <- res$size
          zsum <- res$zsum
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

  list(improved = improved, best_score = best_score, size = size, zsum = zsum)
}

#' Greedy removal pass
#'
#' @param vstate     Environment with mutable node vectors.
#' @param z_vec      Numeric z-score vector.
#' @param sc_means   MC means.
#' @param sc_stds    MC stds.
#' @param size       Current component size (from expand return value).
#' @param zsum       Current z-score sum (from expand return value).
#' @param best_score Best score so far.
#' @return Updated best_score.
.greedy_removal_idx <- function(vstate, z_vec, sc_means, sc_stds,
                                size, zsum, best_score) {
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
  }
  best_score
}

#' Filter candidates by mutual overlap using an inverted node index
#'
#' Instead of an O(n²) double loop, builds a node→candidate inverted index so
#' each kept component only compares against the small set of candidates that
#' share at least one node with it.  For the tiny components produced at
#' depth=1 this reduces comparisons from O(n²) to O(n × d) where d is the
#' mean number of candidates sharing a node (~9x faster at N=12 000).
#'
#' @param idx_list  List of sorted integer index vectors, one per candidate.
#' @param scores    Numeric scores aligned to \code{idx_list} (desc order).
#' @param threshold Overlap threshold in (0, 1).
#' @param limit     Maximum candidates to keep.
#' @param n_nodes   Total node count.
#' @return Integer vector of kept positions into \code{idx_list}.
.greedy_filter_int <- function(idx_list, scores, threshold, limit, n_nodes) {
  n <- length(idx_list)
  if (n == 0L) {
    return(integer(0L))
  }

  sizes <- lengths(idx_list)

  # Inverted index: node i -> integer vector of candidate positions containing i
  inv <- vector("list", n_nodes)
  for (i in seq_len(n)) {
    for (ni in idx_list[[i]]) {
      inv[[ni]] <- c(inv[[ni]], i)
    }
  }

  to_delete <- logical(n)
  kept <- integer(0L)

  i <- 1L
  while (i <= n && length(kept) < limit) {
    if (!to_delete[i]) {
      s1 <- idx_list[[i]]
      s1_sz <- sizes[i]
      kept <- c(kept, i)

      # Candidates that share at least one node with s1
      partners <- unique(unlist(inv[s1], use.names = FALSE))
      partners <- partners[partners > i & !to_delete[partners]]

      for (j in partners) {
        common <- length(intersect(s1, idx_list[[j]]))
        denom <- if (s1_sz < sizes[j]) s1_sz else sizes[j]
        if (common / denom > threshold) to_delete[j] <- TRUE
      }
    }
    i <- i + 1L
  }

  kept
}

#' Run the greedy active subnetwork search
#'
#' @param network  A network from \code{build_network()}.
#' @param sc       A score context from \code{build_score_context()}.
#' @param params   A list of run parameters.
#' @param verbose  Logical; emit progress messages.
#' @return A list of subnetwork objects.
.greedy_search <- function(network, sc, params, verbose = FALSE) {
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

  # Build integer neighbour index from the edge list
  el <- igraph::as_edgelist(network$g, names = FALSE)
  both <- rbind(el, el[, 2:1, drop = FALSE])
  both <- both[order(both[, 1L]), , drop = FALSE]
  row_starts <- c(1L, which(diff(both[, 1L]) > 0L) + 1L, nrow(both) + 1L)

  nbr_idx <- vector("list", n_nodes)
  for (i in seq_len(n_nodes)) {
    s <- row_starts[i]
    e <- row_starts[i + 1L] - 1L
    nb <- if (s <= e) both[s:e, 2L] else integer(0L)
    nbr_idx[[i]] <- nb
  }

  # Every node is used as a seed (matches the reference algorithm). Restricting
  # seeds to high-scoring nodes is much faster but changes the result set, so it
  # is deliberately not done here, can be revisited at a later date
  promising_seeds <- seq_len(n_nodes)

  n_seeds <- length(promising_seeds)
  lv_template <- logical(n_nodes)
  iv_template <- integer(n_nodes)

  # Keep the best-scoring component for each NODE (node2best), then
  # collect the distinct components afterwards.
  node2best <- vector("list", n_nodes)

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

    vstate <- new.env(parent = emptyenv())
    vstate$comp_members <- lv_template
    vstate$removable_vec <- lv_template
    vstate$node2predecessor <- iv_template
    vstate$node2dep_count <- iv_template

    vstate$comp_members[seed_id] <- TRUE
    vstate$node2dep_count[seed_id] <- 1L

    seed_zsum <- z_vec[seed_id]

    res <- .greedy_expand_idx(
      nbr_idx, z_vec, sc_means, sc_stds,
      within_vec, search_depth, search_depth,
      1L, seed_zsum, 0, -Inf,
      vstate, seed_id
    )

    # Pass size/zsum from expand directly — no O(N) sum() scans needed
    final_best <- .greedy_removal_idx(
      vstate, z_vec, sc_means, sc_stds,
      res$size, res$zsum, res$best_score
    )

    final_idx <- which(vstate$comp_members)

    if (length(final_idx) >= 2L && final_best > 0) {
      for (nd in final_idx) {
        ex <- node2best[[nd]]
        if (is.null(ex) || final_best > ex$score) {
          node2best[[nd]] <- list(idx = final_idx, score = final_best)
        }
      }
    }
  }
  if (verbose) message("100%")

  # Collect the distinct best-per-node components
  seen <- new.env(hash = TRUE, parent = emptyenv())
  for (nd in seq_len(n_nodes)) {
    c <- node2best[[nd]]
    if (!is.null(c)) {
      key <- paste0(c$idx, collapse = "\x01")
      if (is.null(seen[[key]])) {
        seen[[key]] <- list(nodes = nodes[c$idx], idx = c$idx, score = c$score)
      }
    }
  }
  candidates <- as.list(seen)
  if (length(candidates) == 0L) {
    return(list())
  }

  scores <- vapply(candidates, function(x) x$score, numeric(1))
  ord <- order(scores, decreasing = TRUE)
  candidates <- candidates[ord]
  scores <- scores[ord]

  idx_list <- lapply(candidates, function(x) x$idx)
  kept <- .greedy_filter_int(
    idx_list, scores,
    params$gr_overlap_threshold,
    params$gr_subnetwork_num,
    n_nodes
  )

  lapply(candidates[kept], function(x) list(nodes = x$nodes, score = x$score))
}
