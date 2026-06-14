#' Get Active Subnetworks
#'
#' Performs active subnetwork search and filters identified active subnetworks
#' before returning final subnetworks.
#'
#' @param input_for_search input the input data that active subnetwork search uses. The input
#' must be a data frame containing at least these 2 columns: \describe{
#'   \item{GENE}{Gene Symbol}
#'   \item{P_VALUE}{p value obtained through a test, e.g. differential expression/methylation}
#' }
#' @inheritParams return_pin_path
#' @inheritParams filter_active_subnetworks
#' @param search_method algorithm to use when performing active subnetwork
#'  search. Options are greedy search (GR), simulated annealing (SA) or genetic
#'  algorithm (GA) for the search (default = 'GR').
#' @param seed_for_stochastic_methods seed for reproducibility while running active subnetwork search (applies for GR and SA)
#' @param verbose boolean value indicating whether to print messages (default=FALSE)
#' @param start_with_all_positives if TRUE: in GA, adds an individual with all positive
#'  nodes. In SA, initializes candidate solution with all positive nodes. (default = FALSE)
#' @param gene_init_prob For SA and GA, probability of adding a gene in initial solution (default = 0.1)
#' @param sa_initial_temp Initial temperature for SA (default = 1.0)
#' @param sa_final_temp Final temperature for SA (default = 0.01)
#' @param sa_iterations Iteration number for SA (default = 10000)
#' @param ga_population_size Population size for GA (default = 400)
#' @param ga_iterations Iteration number for GA (default = 200)
#' @param ga_crossover_rate Applies crossover with the given probability in GA (default = 1, i.e. always perform crossover)
#' @param ga_mutation_rate For GA, applies mutation with given mutation rate (default = 0, i.e. mutation off)
#' @param gr_max_depth Sets max depth in greedy search, 0 for no limit (default = 1)
#' @param gr_search_depth Search depth in greedy search (default = 1)
#' @param gr_overlap_threshold Overlap threshold for results of greedy search (default = 0.5)
#' @param gr_subnetwork_num Number of subnetworks to be presented in the results (default = 1000)
#'
#' @return A list of genes in every identified active subnetwork that has a score greater than
#' the `score_quan_thr`th quantile and that has at least `sig_gene_thr` affected genes.
#'
#' @export
#'
#' @examples
#' \donttest{
#' processed_df <- example_pathfindR_input[1:15, -2]
#' colnames(processed_df) <- c("GENE", "P_VALUE")
#' GR_snws <- get_active_subnetworks(processed_df, pin_name_path = "KEGG")
#' }
get_active_subnetworks <- function(
  input_for_search, pin_name_path = "Biogrid",
  score_quan_thr = 0.8, sig_gene_thr = 0.02, search_method = "GR",
  seed_for_stochastic_methods = 1234, verbose = FALSE, start_with_all_positives = FALSE, gene_init_prob = 0.1,
  sa_initial_temp = 1, sa_final_temp = 0.01, sa_iterations = 10000, ga_population_size = 400, ga_iterations = 10000,
  ga_crossover_rate = 1, ga_mutation_rate = 0, gr_max_depth = 1, gr_search_depth = 1, gr_overlap_threshold = 0.5,
  gr_subnetwork_num = 1000
) {
  if (!is.data.frame(input_for_search)) {
    stop("`input_for_search` should be data frame")
  }
  cnames <- c("GENE", "P_VALUE")
  if (any(!cnames %in% colnames(input_for_search))) {
    stop("`input_for_search` should contain the columns ", paste(dQuote(cnames),
      collapse = ","
    ))
  }

  valid_mets <- c("GR", "SA", "GA")
  if (!search_method %in% valid_mets) {
    stop("`search_method` should be one of ", paste(dQuote(valid_mets), collapse = ", "))
  }

  if (!is.logical(verbose)) {
    stop("`verbose` should be either TRUE or FALSE")
  }

  if (!is.logical(start_with_all_positives)) {
    stop("`start_with_all_positives` should be either TRUE or FALSE")
  }

  pin_path <- return_pin_path(pin_name_path)

  input_for_search$GENE <- base::toupper(input_for_search$GENE)
  experiment_df <- input_for_search[, c("GENE", "P_VALUE")]

  ############ Run active Subnetwork Search
  identified_active_snws <- active_subnetwork_search(
    pin_path = pin_path,
    experiment = experiment_df,
    method = search_method,
    start_with_all_positives = start_with_all_positives,
    gene_init_prob = gene_init_prob,
    sa_initial_temp = sa_initial_temp,
    sa_final_temp = sa_final_temp,
    sa_iterations = sa_iterations,
    ga_population_size = ga_population_size,
    ga_iterations = ga_iterations,
    ga_crossover_rate = ga_crossover_rate,
    ga_mutation_rate = ga_mutation_rate,
    gr_max_depth = gr_max_depth,
    gr_search_depth = gr_search_depth,
    gr_overlap_threshold = gr_overlap_threshold,
    gr_subnetwork_num = gr_subnetwork_num,
    seed = seed_for_stochastic_methods,
    verbose = verbose
  )

  ############ Parse and filter active subnetworks
  filtered_snws <- filter_active_subnetworks(
    active_snws = identified_active_snws, sig_genes_vec = experiment_df$GENE,
    score_quan_thr = score_quan_thr, sig_gene_thr = sig_gene_thr
  )

  if (is.null(filtered_snws)) {
    snws <- list()
  } else {
    snws <- filtered_snws$subnetworks
  }
  message(paste0("Found ", length(snws), " active subnetworks\n\n"))

  return(snws)
}
