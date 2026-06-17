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
#' independently with probability \code{params$gene_init_prob}. The temperature
#' decays geometrically from \code{sa_initial_temp} to \code{sa_final_temp}
#' over \code{sa_iterations} steps.
#'
#' @param network A network from \code{build_network()}.
#' @param sc A score context from \code{build_score_context()}.
#' @param params A list of run parameters.
#' @param verbose Logical; emit progress messages.
#'
#' @return A list of subnetwork objects (all connected components of the final
#'   solution, sorted by score descending).
.simulated_annealing <- function(network, sc, params, verbose = FALSE) {
  set.seed(params$seed)
  nodes <- network$nodes
  N <- length(nodes)
  if (N == 0L) {
    return(list())
  }

  z <- sc$z[nodes]
  on <- stats::setNames(logical(N), nodes)

  if (isTRUE(params$start_with_all_positives)) {
    on[z > 0] <- TRUE
  } else {
    on[stats::runif(N) < params$gene_init_prob] <- TRUE
  }

  # Precomputed inputs for the C++ component scorer (built once).
  z_vec <- as.numeric(z)
  means <- as.numeric(sc$means)
  stds <- as.numeric(sc$stds)
  offsets <- network$csr_offsets
  nbrs <- network$csr_nbrs

  # Current solution is tracked as its sorted-descending component scores; this
  # is all the rank-by-rank acceptance test below needs. The actual node
  # membership is reconstructed only once, at the very end. component_scores_sorted()
  # reproduces .sort_subnetworks_desc(.find_subnetworks(...)) score-for-score.
  cur_scores <- component_scores_sorted(on, offsets, nbrs, z_vec, means, stds)

  T0 <- params$sa_initial_temp
  Tf <- params$sa_final_temp
  total <- params$sa_iterations
  Temp <- T0
  temp_step <- 1 - (Tf / T0)^(1 / total)

  percent <- 0L
  for (iter in seq_len(total)) {
    if (verbose) {
      new_percent <- (100L * (iter - 1L)) %/% total
      if (new_percent > percent) {
        percent <- new_percent
        message(percent, "%")
      }
    }

    idx <- sample.int(N, 1L)
    on[idx] <- !on[idx]

    new_scores <- component_scores_sorted(on, offsets, nbrs, z_vec, means, stds)

    # If the current solution is empty, always accept an improving move.
    if (length(cur_scores) == 0L) {
      if (length(new_scores) > 0L) {
        cur_scores <- new_scores
      } else {
        on[idx] <- !on[idx]   # revert: neither state has subnetworks
      }
      Temp <- Temp * (1 - temp_step)
      next
    }

    decision <- FALSE
    keep <- FALSE
    m <- min(length(cur_scores), length(new_scores))
    k <- 1L
    while (!decision && k <= m) {
      delta <- new_scores[k] - cur_scores[k]
      if (delta > 0.001) {
        keep <- TRUE
        decision <- TRUE
      } else if (stats::runif(1) < exp(delta / Temp)) {
        # Accept worse move with Boltzmann probability exp(delta/T).
        # delta <= 0 here, so this probability is in (0, 1] and falls
        # as temperature decreases — the correct SA behaviour.
        keep <- TRUE
        decision <- TRUE
      } else {
        keep <- FALSE
        decision <- TRUE
      }
      k <- k + 1L
    }

    if (keep) {
      cur_scores <- new_scores
    } else {
      on[idx] <- !on[idx] # revert the toggle
    }

    Temp <- Temp * (1 - temp_step)
  }
  if (verbose) message("100%")

  # Materialise the final accepted solution's subnetworks once.
  .sort_subnetworks_desc(.find_subnetworks(network, sc, nodes[on]))
}
