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
#' @param seedForRandom seed for reproducibility while running active subnetwork search (applies for GR and SA)
#' @param silent_option boolean value indicating whether to print the messages
#' to the console (FALSE) or not (TRUE, this will print to a temp. file) during
#' active subnetwork search (default = TRUE). This option was added because
#' during parallel runs, the console messages get disorderly printed.
#' @param use_all_positives if TRUE: in GA, adds an individual with all positive
#'  nodes. In SA, initializes candidate solution with all positive nodes. (default = FALSE)
#' @param geneInitProbs For SA and GA, probability of adding a gene in initial solution (default = 0.1)
#' @param saTemp0 Initial temperature for SA (default = 1.0)
#' @param saTemp1 Final temperature for SA (default = 0.01)
#' @param saIter Iteration number for SA (default = 10000)
#' @param gaPop Population size for GA (default = 400)
#' @param gaIter Iteration number for GA (default = 200)
#' @param gaThread Number of threads to be used in GA (default = 5)
#' @param gaCrossover Applies crossover with the given probability in GA (default = 1, i.e. always perform crossover)
#' @param gaMut For GA, applies mutation with given mutation rate (default = 0, i.e. mutation off)
#' @param grMaxDepth Sets max depth in greedy search, 0 for no limit (default = 1)
#' @param grSearchDepth Search depth in greedy search (default = 1)
#' @param grOverlap Overlap threshold for results of greedy search (default = 0.5)
#' @param grSubNum Number of subnetworks to be presented in the results (default = 1000)
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
  seedForRandom = 1234, silent_option = TRUE, use_all_positives = FALSE, geneInitProbs = 0.1,
  saTemp0 = 1, saTemp1 = 0.01, saIter = 10000, gaPop = 400, gaIter = 10000, gaThread = 5,
  gaCrossover = 1, gaMut = 0, grMaxDepth = 1, grSearchDepth = 1, grOverlap = 0.5,
  grSubNum = 1000
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

  # search_method
  valid_mets <- c("GR", "SA", "GA")
  if (!search_method %in% valid_mets) {
    stop("`search_method` should be one of ", paste(dQuote(valid_mets), collapse = ", "))
  }

  # silent_option
  if (!is.logical(silent_option)) {
    stop("`silent_option` should be either TRUE or FALSE")
  }

  # use_all_positives
  if (!is.logical(use_all_positives)) {
    stop("`use_all_positives` should be either TRUE or FALSE")
  }

  pin_path <- return_pin_path(pin_name_path)

  input_for_search$GENE <- base::toupper(input_for_search$GENE)
  experiment_df <- input_for_search[, c("GENE", "P_VALUE")]

  ############ Run active Subnetwork Search
  identified_active_snws <- active_subnetwork_search(
    pin_path = pin_path,
    experiment = experiment_df,
    method = search_method,
    start_with_all_positives = use_all_positives,
    gene_init_prob = geneInitProbs,
    sa_initial_temp = saTemp0,
    sa_final_temp = saTemp1,
    sa_iterations = saIter,
    ga_population_size = gaPop,
    ga_iterations = gaIter,
    ga_crossover_rate = gaCrossover,
    ga_mutation_rate = gaMut,
    gr_max_depth = grMaxDepth,
    gr_search_depth = grSearchDepth,
    gr_overlap_threshold = grOverlap,
    gr_subnetwork_num = grSubNum,
    seed = seedForRandom,
    verbose = !silent_option
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
