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

  # Neighbour iteration order is taken directly from build_network()'s
  # reconstruction of Java's per-node HashSet iteration order (network$nbr_idx,
  # 1-based ids aligned to `nodes`). This replaces the earlier heuristic that
  # sorted neighbours by ascending z-score: to match the Java reference exactly
  # the greedy expansion must visit neighbours in the same order Java does, which
  # is the HashSet order, not a score-based order. `nodes` is already in Java's
  # networkNodeList order, so seeds are processed in the same order too.
  nbr_idx <- network$nbr_idx
  if (is.null(nbr_idx)) {
    stop("network$nbr_idx is missing; rebuild the network with build_network().")
  }

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
