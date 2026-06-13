# =============================================================================
# Greedy active subnetwork search. For every node it grows a high-scoring
# component by repeatedly adding the most promising neighbours (within a depth
# limit) and then attempts to prune the component again. The best component
# seen for each node is kept, and the resulting set of components is filtered
# by score, size and mutual overlap.
# =============================================================================

#' Mark every node within \code{depth} hops of \code{current} using integer IDs
#'
#' @param nbr_idx A list of integer vectors representing node neighbors.
#' @param within_vec A logical vector acting as a set of reachable node IDs.
#' @param current The integer node ID currently being expanded.
#' @param depth Remaining depth budget.
#'
#' @return Invisibly \code{NULL}; \code{within_vec} is modified in place via assignment.
.greedy_init_max_depth_idx <- function(nbr_idx, within_vec, current, depth) {
  within_vec[current] <<- TRUE
  if (depth > 0) {
    for (nb in nbr_idx[[current]]) {
      if (!within_vec[nb]) {
        .greedy_init_max_depth_idx(nbr_idx, within_vec, nb, depth - 1L)
      }
    }
  }
  invisible(NULL)
}

#' Recursive greedy expansion of a component using flat vectors
#'
#' @param nbr_idx A list of integer vectors of neighbors.
#' @param sc A score context.
#' @param within_vec Either NULL or a logical vector of allowed nodes.
#' @param search_depth The depth budget restored on every improvement.
#' @param depth The remaining depth budget.
#' @param comp_state A list environment holding mutable running state:
#'   \code{size}, \code{zsum}, \code{score}, and \code{best_score}.
#' @param comp_members A logical vector tracking if a node ID is in the component.
#' @param last_added The integer node ID added most recently.
#' @param removable_vec A logical vector acting as a set of removable node IDs.
#' @param node2predecessor An integer vector tracking the path.
#' @param node2dependent_count An integer vector tracking dependency counts.
#'
#' @return Logical; whether the best score improved on this branch.
.greedy_expand_idx <- function(nbr_idx, sc, within_vec, search_depth, depth,
                               comp_state, comp_members, last_added, removable_vec,
                               node2predecessor, node2dependent_count) {
  improved <- FALSE

  if (comp_state$score > comp_state$best_score) {
    depth <- search_depth
    improved <- TRUE
    comp_state$best_score <- comp_state$score
  }

  if (depth > 0) {
    any_improved <- FALSE
    removable_vec[last_added] <- FALSE
    dependent_count <- 0L

    for (new_neighbor in nbr_idx[[last_added]]) {
      within_ok <- is.null(within_vec) || within_vec[new_neighbor]
      if (within_ok && !comp_members[new_neighbor]) {
        # Inlined component add step
        comp_state$size <- comp_state$size + 1L
        comp_members[new_neighbor] <- TRUE
        comp_state$zsum <- comp_state$zsum + sc$z_vec[new_neighbor]
        comp_state$score <- .score_subnetwork(sc, comp_state$size, comp_state$zsum, TRUE)
        removable_vec[new_neighbor] <- TRUE

        this_improved <- .greedy_expand_idx(
          nbr_idx, sc, within_vec, search_depth, depth - 1L,
          comp_state, comp_members, new_neighbor, removable_vec,
          node2predecessor, node2dependent_count
        )

        if (!this_improved) {
          # Inlined component remove step
          comp_state$size <- comp_state$size - 1L
          comp_members[new_neighbor] <- FALSE
          comp_state$zsum <- comp_state$zsum - sc$z_vec[new_neighbor]
          comp_state$score <- .score_subnetwork(sc, comp_state$size, comp_state$zsum, TRUE)
          removable_vec[new_neighbor] <- FALSE
        } else {
          dependent_count <- dependent_count + 1L
          any_improved <- TRUE
          node2predecessor[new_neighbor] <- last_added
        }
      }
    }

    improved <- improved || any_improved
    if (dependent_count > 0L) {
      removable_vec[last_added] <- FALSE
      node2dependent_count[last_added] <- dependent_count
    }
  }

  improved
}

#' Greedy removal pass using vector indexing
#'
#' @return Invisibly \code{NULL}.
.greedy_removal_idx <- function(comp_state, comp_members, removable_vec,
                                node2predecessor, node2dependent_count, sc) {
  snapshot <- which(removable_vec)
  for (current in snapshot) {
    # Inlined component remove step
    comp_state$size <- comp_state$size - 1L
    comp_members[current] <- FALSE
    comp_state$zsum <- comp_state$zsum - sc$z_vec[current]
    comp_state$score <- .score_subnetwork(sc, comp_state$size, comp_state$zsum, TRUE)

    if (comp_state$score > comp_state$best_score) {
      comp_state$best_score <- comp_state$score
      predecessor <- node2predecessor[current]
      if (predecessor > 0L) {
        dep_count <- node2dependent_count[predecessor] - 1L
        if (dep_count == 0L) {
          removable_vec[predecessor] <- TRUE
        } else {
          node2dependent_count[predecessor] <- dep_count
        }
      }
    } else {
      # Inlined component put back step
      comp_state$size <- comp_state$size + 1L
      comp_members[current] <- TRUE
      comp_state$zsum <- comp_state$zsum + sc$z_vec[current]
      comp_state$score <- .score_subnetwork(sc, comp_state$size, comp_state$zsum, TRUE)
    }
  }
  invisible(NULL)
}

#' Filter subnetworks by mutual overlap
#'
#' Walks the score-sorted list keeping a subnetwork and discarding any later
#' subnetwork whose overlap with a kept one exceeds the threshold, up to a
#' maximum number of results.
#'
#' @param comps components
#' @param ovelap_threshold overlap threshold
#' @param subnetwork_num number of subnetworks
#' @return filtered components
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
  while (i < n && length(filtered) < subnetwork_num) {
    if (!to_delete[i]) {
      s1_ids <- ids[[i]]
      s1_size <- sizes[i]
      filtered[[length(filtered) + 1L]] <- comps[[i]]

      mark[s1_ids] <- TRUE
      for (j in (i + 1L):n) {
        if (!to_delete[j]) {
          common <- sum(mark[ids[[j]]])
          size <- if (s1_size < sizes[j]) s1_size else sizes[j]
          if (common / size > overlap_threshold) {
            to_delete[j] <- TRUE
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
#' @param network A network from \code{.build_network()}.
#' @param sc A score context.
#' @param params A list of run parameters.
#' @param verbose Logical; emit progress messages.
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

  # --- Integer Mapping Initialization ---
  node_map <- stats::setNames(seq_len(n_nodes), nodes)
  # Only allow neighbors that aren't actively penalizing the subnetwork score
  nbr_idx <- lapply(seq_len(n_nodes), function(i) {
    all_nbrs <- unname(node_map[network$nbr[[nodes[i]]]])
    # Filter out neighbors with terrible z-scores (taken here as z < -1) to prevent searching dead space
    all_nbrs[sc$z_vec[all_nbrs] > -1.0]
  })
  sc$z_vec <- as.numeric(sc$z[nodes]) # Faster aligned vector access

  # Fast structure to store the best component score seen for each integer node ID
  best_score_per_node <- rep(-Inf, n_nodes)
  best_comp_per_node <- vector("list", n_nodes)

  percent <- 0L
  for (seed_id in seq_len(n_nodes)) {
    # --- Fast Percentile Seed Selection ---
    # Target the top 5% highest z-scores (95th percentile)
    percentile_target <- 0.95
    z_cutoff <- stats::quantile(sc$z_vec, probs = percentile_target, na.rm = TRUE)

    # Identify seeds matching the threshold
    promising_seeds <- which(sc$z_vec >= z_cutoff)

    # Safety fallbacks:
    # 1. Ensure we don't accidentally select thousands of identical zero/negative ties
    # 2. Ensure we have at least a few seeds if the network is tiny
    if (length(promising_seeds) > (n_nodes * 0.20)) {
      # If a massive tie-break group blows past 20% of the network, take a hard top slice instead
      promising_seeds <- order(sc$z_vec, decreasing = TRUE)[1:as.integer(n_nodes * 0.05)]
    } else if (length(promising_seeds) < 50L && n_nodes >= 50L) {
      promising_seeds <- order(sc$z_vec, decreasing = TRUE)[1:50L]
    }

    percent <- 0L
    n_seeds <- length(promising_seeds)
    for (i in seq_along(promising_seeds)) {
      seed_id <- promising_seeds[i]
      seed_name <- nodes[seed_id]

      if (verbose) {
        new_percent <- (100L * (seed_id - 1L)) %/% n_nodes
        if (new_percent > percent) {
          percent <- new_percent
          message(percent, "% of selected seeds checked (", i, "/", n_seeds, ")")
        }
      }

      if (max_depth == 0L) {
        within_vec <- NULL
      } else {
        within_vec <- logical(n_nodes)
        .greedy_init_max_depth_idx(nbr_idx, within_vec, seed_id, max_depth)
      }

      # Tracking states initialized as flat structures
      comp_state <- list(
        size = 1L,
        zsum = sc$z_vec[seed_id],
        score = .score_subnetwork(sc, 1L, sc$z_vec[seed_id], TRUE),
        best_score = -Inf
      )

      comp_members <- logical(n_nodes)
      comp_members[seed_id] <- TRUE
      removable_vec <- logical(n_nodes)
      node2predecessor <- integer(n_nodes)
      node2dependent_count <- integer(n_nodes)
      node2dependent_count[seed_id] <- 1L

      .greedy_expand_idx(
        nbr_idx, sc, within_vec, search_depth, search_depth,
        comp_state, comp_members, seed_id, removable_vec,
        node2predecessor, node2dependent_count
      )

      .greedy_removal_idx(
        comp_state, comp_members, removable_vec,
        node2predecessor, node2dependent_count, sc
      )

      # If the constructed component beats previous paths holding these nodes, capture it
      final_nodes_idx <- which(comp_members)
      final_score <- comp_state$best_score

      if (length(final_nodes_idx) > 0L) {
        comp_obj <- list(nodes = nodes[final_nodes_idx], score = final_score)
        for (node_idx in final_nodes_idx) {
          if (best_score_per_node[node_idx] < final_score) {
            best_score_per_node[node_idx] <- final_score
            best_comp_per_node[[node_idx]] <- comp_obj
          }
        }
      }
    }
  }
  if (verbose) message("100%")

  # Extract valid components and clean up null items
  comps <- best_comp_per_node[!vapply(best_comp_per_node, is.null, logical(1))]

  # Collapse identical components
  if (length(comps) > 0L) {
    sig <- vapply(
      comps,
      function(c) paste0(sort(c$nodes), collapse = "\u0001"),
      character(1)
    )
    comps <- comps[!duplicated(sig)]
  }

  if (length(comps) == 0L) {
    return(list())
  }

  scores <- vapply(comps, function(c) c$score, numeric(1))
  comps <- comps[order(scores, decreasing = TRUE)]

  comps <- comps[vapply(comps, function(c) c$score, numeric(1)) > 0]
  comps <- comps[vapply(comps, function(c) length(c$nodes), integer(1)) >= 2L]

  .greedy_filter(comps, params$gr_overlap_threshold, params$gr_subnetwork_num)
}
