# =============================================================================
# Greedy active subnetwork search. For every node it grows a high-scoring
# component by repeatedly adding the most promising neighbours (within a depth
# limit) and then attempts to prune the component again. The best component
# seen for each node is kept, and the resulting set of components is filtered
# by score, size and mutual overlap.
# =============================================================================

#' Mark every node within \code{depth} hops of \code{current}
#'
#' Populates the \code{within} environment (used as a hash set) with the names
#' of nodes that lie within the allowed depth of the starting node.
#'
#' @param network A network from \code{.build_network()}.
#' @param within An environment used as a set of reachable node names.
#' @param current The node currently being expanded.
#' @param depth Remaining depth budget.
#'
#' @return Invisibly \code{NULL}; \code{within} is modified in place.
.greedy_init_max_depth <- function(network, within, current, depth) {
  assign(current, TRUE, envir = within)
  if (depth > 0) {
    for (nb in network$nbr[[current]]) {
      if (!exists(nb, envir = within, inherits = FALSE)) {
        .greedy_init_max_depth(network, within, nb, depth - 1L)
      }
    }
  }
  invisible(NULL)
}

#' Recursive greedy expansion of a component
#'
#' Adds neighbouring nodes as long as doing so improves the best score seen
#' from the current starting point; the search depth is refreshed each time an
#' improvement is found.
#'
#' @param network A network from \code{.build_network()}.
#' @param gs A greedy-state environment holding \code{bestScore},
#'   \code{node2Predecessor} and \code{node2DependentCount}.
#' @param within Either \code{NULL} (no depth limit) or an environment of nodes
#'   within the allowed depth.
#' @param search_depth The depth budget restored on every improvement.
#' @param depth The remaining depth budget.
#' @param comp A mutable component environment.
#' @param last_added The node added most recently.
#' @param removable An environment used as a set of removable nodes.
#'
#' @return Logical; whether the best score improved on this branch.
.greedy_expand <- function(network, gs, within, search_depth, depth,
                           comp, last_added, removable) {
  improved <- FALSE

  if (comp$score > gs$bestScore) {
    depth <- search_depth
    improved <- TRUE
    gs$bestScore <- comp$score
  }

  if (depth > 0) {
    any_improved <- FALSE
    if (exists(last_added, envir = removable, inherits = FALSE)) {
      rm(list = last_added, envir = removable)
    }
    dependent_count <- 0L

    for (new_neighbor in network$nbr[[last_added]]) {
      within_ok <- is.null(within) ||
        exists(new_neighbor, envir = within, inherits = FALSE)
      if (within_ok && !.component_contains(comp, new_neighbor)) {
        .component_add(comp, new_neighbor)
        assign(new_neighbor, TRUE, envir = removable)

        this_improved <- .greedy_expand(
          network, gs, within, search_depth,
          depth - 1L, comp, new_neighbor,
          removable
        )
        if (!this_improved) {
          .component_remove(comp, new_neighbor)
          if (exists(new_neighbor, envir = removable, inherits = FALSE)) {
            rm(list = new_neighbor, envir = removable)
          }
        } else {
          dependent_count <- dependent_count + 1L
          any_improved <- TRUE
          assign(new_neighbor, last_added, envir = gs$node2Predecessor)
        }
      }
    }

    improved <- improved || any_improved
    if (dependent_count > 0L) {
      if (exists(last_added, envir = removable, inherits = FALSE)) {
        rm(list = last_added, envir = removable)
      }
      assign(last_added, dependent_count, envir = gs$node2DependentCount)
    }
  }

  improved
}

#' Greedy removal pass
#'
#' Tries to remove each removable node; a removal is accepted when it improves
#' the best score, otherwise the node is put back.
#'
#' @param gs A greedy-state environment.
#' @param comp A mutable component environment.
#' @param removable An environment used as a set of removable nodes.
#'
#' @return Invisibly \code{NULL}.
.greedy_removal <- function(gs, comp, removable) {
  snapshot <- ls(removable)
  for (current in snapshot) {
    .component_remove(comp, current)
    score <- comp$score
    if (score > gs$bestScore) {
      gs$bestScore <- score
      if (exists(current, envir = gs$node2Predecessor, inherits = FALSE)) {
        predecessor <- get(current, envir = gs$node2Predecessor)
        dependent_count <- get(predecessor, envir = gs$node2DependentCount) - 1L
        if (dependent_count == 0L) {
          assign(predecessor, TRUE, envir = removable)
        } else {
          assign(predecessor, dependent_count, envir = gs$node2DependentCount)
        }
      }
    } else {
      .component_add(comp, current)
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
#' @param comps A score-sorted list of subnetwork objects.
#' @param overlap_threshold Overlap above which a subnetwork is discarded.
#' @param subnetwork_num Maximum number of subnetworks to keep.
#'
#' @return A filtered list of subnetwork objects.
.greedy_filter <- function(comps, overlap_threshold, subnetwork_num) {
  n <- length(comps)
  if (n == 0L) {
    return(list())
  }

  to_delete <- logical(n)
  filtered <- list()

  i <- 1L
  while (i < n && length(filtered) < subnetwork_num) {
    if (!to_delete[i]) {
      s1 <- comps[[i]]
      s1_nodes <- s1$nodes
      s1_size <- length(s1_nodes)
      filtered[[length(filtered) + 1L]] <- s1

      for (j in (i + 1L):n) {
        if (!to_delete[j]) {
          s2 <- comps[[j]]
          common <- length(intersect(s1_nodes, s2$nodes))
          size <- min(s1_size, length(s2$nodes))
          overlap <- common / size
          if (overlap > overlap_threshold) {
            to_delete[j] <- TRUE
          }
        }
      }
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

  node2best <- new.env(parent = emptyenv())

  percent <- 0L
  for (node_no in seq_len(n_nodes)) {
    seed <- nodes[node_no]

    if (verbose) {
      new_percent <- (100L * (node_no - 1L)) %/% n_nodes
      if (new_percent > percent) {
        percent <- new_percent
        message(percent, "% of seeds checked")
      }
    }

    if (max_depth == 0L) {
      within <- NULL
    } else {
      within <- new.env(parent = emptyenv())
      .greedy_init_max_depth(network, within, seed, max_depth)
    }

    comp <- .component_new(sc, seed)
    gs <- new.env(parent = emptyenv())
    gs$bestScore <- -Inf
    gs$node2Predecessor <- new.env(parent = emptyenv())
    gs$node2DependentCount <- new.env(parent = emptyenv())
    assign(seed, 1L, envir = gs$node2DependentCount)
    removable <- new.env(parent = emptyenv())

    .greedy_expand(
      network, gs, within, search_depth, search_depth,
      comp, seed, removable
    )
    .greedy_removal(gs, comp, removable)

    comp_nodes <- comp$nodes
    comp_score <- comp$score
    for (node in comp_nodes) {
      old <- if (exists(node, envir = node2best, inherits = FALSE)) {
        get(node, envir = node2best)
      } else {
        NULL
      }
      if (is.null(old) || old$score < comp_score) {
        assign(node, list(nodes = comp_nodes, score = comp_score),
          envir = node2best
        )
      }
    }
  }
  if (verbose) message("100%")

  keys <- ls(node2best)
  comps <- lapply(keys, function(k) get(k, envir = node2best))

  # collapse identical components (the same component is referenced by each
  # of its member nodes)
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
