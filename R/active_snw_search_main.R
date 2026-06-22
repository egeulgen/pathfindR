#' Active subnetwork search
#'
#' Searches a molecular interaction network for connected subnetworks of genes
#' that are jointly enriched for low experimental p-values ("active modules").
#' Gene p-values are converted to z-scores, subnetwork scores are calibrated
#' against a Monte-Carlo background, and one of three search strategies is used
#' to find high-scoring connected subnetworks.
#'
#' @param network A network list as returned by \code{build_network()}.
#' @param score_context A score context list as returned by \code{build_score_context()}.
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
#' @examples
#' \dontrun{
#' # Write a minimal SIF file
#' sif <- data.frame(
#'   source = c("A", "A", "B", "C", "D"),
#'   type   = "interacts",
#'   target = c("B", "C", "C", "D", "E")
#' )
#' sif_path <- tempfile(fileext = ".sif")
#' write.table(sif, sif_path,
#'   sep = "\t", row.names = FALSE, col.names = FALSE,
#'   quote = FALSE
#' )
#'
#' experiment <- data.frame(
#'   gene   = c("A", "B", "C", "D", "E"),
#'   pvalue = c(0.001, 0.002, 0.001, 0.5, 0.6)
#' )
#'
#' params <- list(
#'   start_with_all_positives = FALSE,
#'   gene_init_prob           = 0.1,
#'   p_for_nonsignificant     = 0.5,
#'   sa_initial_temp          = 1.0,
#'   sa_final_temp            = 0.01,
#'   sa_iterations            = 10000L,
#'   ga_population_size       = 400L,
#'   ga_iterations            = 200L,
#'   ga_crossover_rate        = 1,
#'   ga_mutation_rate         = 0,
#'   gr_max_depth             = 1L,
#'   gr_search_depth          = 1L,
#'   gr_overlap_threshold     = 0.5,
#'   gr_subnetwork_num        = 1000,
#'   seed                     = 1234L
#' )
#'
#' network <- build_network(sif_path)
#' sc <- build_score_context(network, experiment, params)
#' snws <- active_subnetwork_search(network, sc, method = "GR", params = params)
#' }
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
