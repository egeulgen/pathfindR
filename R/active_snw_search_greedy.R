# =============================================================================
# Greedy active subnetwork search. Each node is used as a seed and the current
# subnetwork is recursively expanded through neighbouring nodes whenever the
# calibrated score improves. A subsequent pruning step removes dispensable
# nodes, after which duplicate and highly overlapping subnetworks are filtered
# to yield the final ranked candidates.
# =============================================================================

#' Run the greedy active subnetwork search
#'
#' @param network  A network from \code{build_network()}.
#' @param score_context A score context from \code{build_score_context()}.
#' @param params   A list of run parameters.
#' @param verbose  Logical; emit progress messages.
#' @return A list of subnetwork objects.
.greedy_search <- function(network, score_context, params, verbose = FALSE) {
  nodes <- network$nodes
  n_nodes <- length(nodes)
  if (n_nodes == 0L) {
    return(list())
  }

  max_depth <- as.integer(params$gr_max_depth)
  search_depth <- as.integer(params$gr_search_depth)
  ol_threshold <- as.numeric(params$gr_overlap_threshold)
  max_output <- as.integer(params$gr_subnetwork_num)

  z_vec <- as.numeric(score_context$z[nodes])
  sc_means <- as.numeric(score_context$means)
  sc_stds <- as.numeric(score_context$stds)

  # Build the neighbour index sorted by ascending z-score. The greedy expansion
  # tries neighbours in order and keeps the first one whose addition improves
  # the score. Placing low-z (non-significant) neighbours first means they are
  # evaluated and quickly rejected, so the expansion reaches high-z neighbours
  # while best_score is still low -- giving them the best chance to register as
  # score-improving and be included. This is order-independent of the arbitrary
  # edge-list order from igraph and materially increases overlap with the
  # original Java implementation of active subnetwork search
  # (pathfindR <= v2.X.Y).
  nbr_raw <- lapply(igraph::adjacent_vertices(network$g, igraph::V(network$g)), as.integer)
  nbr_idx <- lapply(nbr_raw, function(ids) ids[order(z_vec[ids], decreasing = FALSE)])


  if (verbose) message("Running greedy search in C++...")
  candidates <- run_greedy_search(
    nbr_idx            = nbr_idx,
    z_vec              = z_vec,
    sc_means           = sc_means,
    sc_stds            = sc_stds,
    node_names         = nodes,
    max_depth          = max_depth,
    search_depth       = search_depth,
    n_nodes            = n_nodes,
    overlap_threshold  = ol_threshold,
    subnetwork_num     = max_output
  )
  return(candidates)
}
