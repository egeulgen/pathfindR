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

  # Directly hand Rcpp the pre-built neighbor ID lists from network$nbr
  # Converting strings to integer indices based on name2id map
  nbr_idx <- lapply(nodes, function(nd) {
    as.integer(network$name2id[network$nbr[[nd]]])
  })

  z_vec <- as.numeric(score_context$z[nodes])
  sc_means <- as.numeric(score_context$means)
  sc_stds <- as.numeric(score_context$stds)

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
