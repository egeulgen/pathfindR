# =============================================================================
# Simulated annealing active subnetwork search. Nodes are switched on/off; at
# each step a random node is toggled and the resulting set of subnetworks is
# compared rank-by-rank against the current one, accepting the change with a
# temperature-dependent probability.
# =============================================================================

#' Run the simulated annealing active subnetwork search
#'
#' The candidate solution starts either with all positive-z-score nodes on
#' (when \code{params$start_with_all_positives}) or with each node switched on
#' independently with probability \code{params$gene_init_prob}. A random node is
#' toggled each step and the resulting subnetworks are compared rank-by-rank
#' against the current ones, accepting the change with a temperature-dependent
#' probability. The temperature decays geometrically from \code{sa_initial_temp}
#' to \code{sa_final_temp} over \code{sa_iterations} steps.
#'
#' The search loop runs in C++ (\code{run_simulated_annealing}) and reproduces
#' the Java reference exactly: the same \code{java.util.Random} stream (seeded
#' with \code{params$seed}), the same initial-state draw order (node order), the
#' same \code{nextInt}-based node toggling, and the same acceptance walk over the
#' score-ranked subnetworks (advance to the next pair when a worse move passes
#' its Boltzmann draw; accept on a >0.001 improvement; otherwise reject).
#'
#' @param network A network from \code{build_network()} (provides the CSR
#'   adjacency \code{csr_offsets} / \code{csr_nbrs}).
#' @param sc A score context from \code{build_score_context()}.
#' @param params A list of run parameters.
#' @param verbose Logical; emit progress messages.
#'
#' @return A list of subnetwork objects (all connected components of the final
#'   solution, sorted by score descending).
.simulated_annealing <- function(network, sc, params, verbose = FALSE) {
  nodes <- network$nodes
  N <- length(nodes)
  if (N == 0L) {
    return(list())
  }

  if (verbose) message("Running simulated annealing in C++...")

  z_vec <- as.numeric(sc$z[nodes])
  means <- as.numeric(sc$means)
  stds <- as.numeric(sc$stds)
  offsets <- network$csr_offsets
  nbrs <- network$csr_nbrs

  # Entire SA loop runs in C++ with the java.util.Random replica, returning the
  # final accepted "on" mask (logical, aligned to node order).
  on_vec <- run_simulated_annealing(
    csr_offsets               = offsets,
    csr_nbrs                  = nbrs,
    z                         = z_vec,
    means                     = means,
    stds                      = stds,
    n_nodes                   = N,
    gene_init_prob            = as.numeric(params$gene_init_prob),
    start_with_all_positives  = isTRUE(params$start_with_all_positives),
    sa_initial_temp           = as.numeric(params$sa_initial_temp),
    sa_final_temp             = as.numeric(params$sa_final_temp),
    sa_iterations             = as.integer(params$sa_iterations),
    seed                      = as.integer(params$seed)
  )

  on <- stats::setNames(as.logical(on_vec), nodes)

  if (verbose) message("100%")

  # Materialise the final accepted solution's subnetworks once.
  .sort_subnetworks_desc(.find_subnetworks(network, sc, nodes[on]))
}
