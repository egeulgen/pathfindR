#' Active subnetwork search
#'
#' Searches a molecular interaction network for connected subnetworks of genes
#' that are jointly enriched for low experimental p-values ("active modules").
#' Gene p-values are converted to z-scores, subnetwork scores are calibrated
#' against a Monte-Carlo background, and one of three search strategies is used
#' to find high-scoring connected subnetworks.
#'
#' @param network A network list as returned by \code{.build_network()}.
#' @param score_context A score context list as returned by \code{.build_score_context()}.
#' @param method Search strategy: \code{"GR"} (greedy, the default),
#'   \code{"SA"} (simulated annealing) or \code{"GA"} (genetic algorithm).
#' @param params A fully-formed params list controlling the search, as built
#'   by \code{get_active_subnetworks()}.
#' @param verbose Logical; emit progress messages.
#'
#' @return A list of subnetwork objects sorted by score descending, each with
#'   elements \code{nodes} (character vector of gene names) and \code{score}.
#'   The list is empty when no positive-scoring subnetwork is found.
#'
#' @seealso \code{\link{get_active_subnetworks}} which builds \code{network},
#'   \code{score_context} and \code{params} and calls this function.
#'
#' @export
active_subnetwork_search <- function(network, score_context, method = c("GR", "SA", "GA"),
                                     params, verbose = FALSE) {
  method <- match.arg(method)

  if (length(network$nodes) == 0L) {
    return(list())
  }

  if (verbose) {
    message("Searching subnetworks")
  }
  subnetworks <- switch(method,
    GR = .greedy_search(network, score_context, params, verbose),
    SA = .simulated_annealing(network, score_context, params, verbose),
    GA = .genetic_algorithm(network, score_context, params, verbose)
  )
  if (verbose) {
    message("Finished searching subnetworks")
  }

  subnetworks <- .sort_subnetworks_desc(subnetworks)
  if (length(subnetworks) > 0L) {
    keep <- vapply(subnetworks, function(s) s$score, numeric(1)) > 0
    subnetworks <- subnetworks[keep]
  }

  return(subnetworks)
}
