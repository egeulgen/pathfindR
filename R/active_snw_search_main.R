#' Active subnetwork search
#'
#' Searches a molecular interaction network for connected subnetworks of genes
#' that are jointly enriched for low experimental p-values ("active modules").
#' Gene p-values are converted to z-scores, subnetwork scores are calibrated
#' against a Monte-Carlo background, and one of three search strategies is used
#' to find high-scoring connected subnetworks.
#'
#' @param pin_path Character string specifying the path to the SIF file.
#'   The file should be tab-delimited and contain at least 3 columns.
#' @param experiment A data frame of gene / p-value pairs. Columns named
#'   \code{gene} and \code{pvalue} are used if present, otherwise the first two
#'   columns are taken as gene and p-value. Gene names are upper-cased and, if a
#'   gene appears more than once, its smallest p-value is kept.
#' @param method Search strategy: \code{"GR"} (greedy, the default),
#'   \code{"SA"} (simulated annealing) or \code{"GA"} (genetic algorithm).
#' @param start_with_all_positives Logical. For \code{SA}/\code{GA}, initialise
#'   the candidate solution with all positive-z-score nodes switched on instead
#'   of random initialisation.
#' @param gene_init_prob Probability of switching a gene on in the random
#'   initial solution for \code{SA}/\code{GA}.
#' @param sa_initial_temp,sa_final_temp,sa_iterations Simulated-annealing
#'   schedule: initial temperature, final temperature and number of iterations.
#' @param ga_population_size,ga_iterations Genetic-algorithm population size and
#'   maximum number of generations.
#' @param ga_crossover_rate,ga_mutation_rate Genetic-algorithm crossover
#'   probability and per-gene mutation rate.
#' @param gr_max_depth Greedy search: maximum distance from the seed at which
#'   nodes may be added (\code{0} = no limit).
#' @param gr_search_depth Greedy search: depth budget restored on every
#'   improvement.
#' @param gr_overlap_threshold Greedy search: subnetworks overlapping a
#'   higher-scoring one by more than this fraction are discarded.
#' @param gr_subnetwork_num Greedy search: maximum number of subnetworks to
#'   return.
#' @param seed Seed for the random number generator, used both for the
#'   Monte-Carlo score calibration and for the search, so that runs are
#'   reproducible.
#' @param verbose Logical; emit progress messages.
#'
#' @return A list of character vectors, each holding the gene names of one
#'   high-scoring subnetwork, ordered from highest to lowest score. The list is
#'   empty when no positive-scoring subnetwork is found.
#'
#' @examples
#' adjacency <- list(
#'   A = c("B", "C"),
#'   B = c("A", "C"),
#'   C = c("A", "B", "D"),
#'   D = c("C", "E"),
#'   E = "D"
#' )
#' experiment <- data.frame(
#'   gene   = c("A", "B", "C", "D", "E"),
#'   pvalue = c(0.001, 0.002, 0.001, 0.5, 0.6)
#' )
#' active_subnetwork_search(pin_path, experiment, method = "GR")
#'
#' @export
active_subnetwork_search <- function(pin_path,
                                     experiment,
                                     method = c("GR", "SA", "GA"),
                                     start_with_all_positives = FALSE,
                                     gene_init_prob = 0.1,
                                     sa_initial_temp = 1.0,
                                     sa_final_temp = 0.01,
                                     sa_iterations = 10000L,
                                     ga_population_size = 400L,
                                     ga_iterations = 200L,
                                     ga_crossover_rate = 1,
                                     ga_mutation_rate = 0,
                                     gr_max_depth = 1L,
                                     gr_search_depth = 1L,
                                     gr_overlap_threshold = 0.5,
                                     gr_subnetwork_num = 1000,
                                     seed = 1234L,
                                     verbose = FALSE) {
  method <- match.arg(method)

  params <- list(
    start_with_all_positives = isTRUE(start_with_all_positives),
    gene_init_prob           = gene_init_prob,
    p_for_nonsignificant     = 0.5,
    sa_initial_temp          = sa_initial_temp,
    sa_final_temp            = sa_final_temp,
    sa_iterations            = as.integer(sa_iterations),
    ga_population_size       = as.integer(ga_population_size),
    ga_iterations            = as.integer(ga_iterations),
    ga_crossover_rate        = ga_crossover_rate,
    ga_mutation_rate         = ga_mutation_rate,
    gr_max_depth             = as.integer(gr_max_depth),
    gr_search_depth          = as.integer(gr_search_depth),
    gr_overlap_threshold     = gr_overlap_threshold,
    gr_subnetwork_num        = gr_subnetwork_num,
    seed                     = seed
  )

  network <- .build_network(pin_path)
  if (length(network$nodes) == 0L) {
    return(list())
  }

  sc <- .build_score_context(network, experiment, params)

  subnetworks <- switch(method,
    GR = .greedy_search(network, sc, params, verbose),
    SA = .simulated_annealing(network, sc, params, verbose),
    GA = .genetic_algorithm(network, sc, params, verbose)
  )

  subnetworks <- .sort_subnetworks_desc(subnetworks)
  if (length(subnetworks) > 0L) {
    keep <- vapply(subnetworks, function(s) s$score, numeric(1)) > 0
    subnetworks <- subnetworks[keep]
  }

  return(subnetworks)
}
