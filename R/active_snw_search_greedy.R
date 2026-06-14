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

  # Build integer neighbour index directly from the igraph edge list.
  el <- igraph::as_edgelist(network$g, names = FALSE)
  both <- rbind(el, el[, 2:1, drop = FALSE])
  both <- both[order(both[, 1L]), , drop = FALSE]
  rs <- c(1L, which(diff(both[, 1L]) > 0L) + 1L, nrow(both) + 1L)
  nbr_idx <- lapply(seq_len(n_nodes), function(i) {
    s <- rs[i]
    e <- rs[i + 1L] - 1L
    if (s <= e) both[s:e, 2L] else integer(0L)
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
