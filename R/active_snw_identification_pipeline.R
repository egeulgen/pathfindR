#' Get Active Subnetworks
#'
#' Performs active subnetwork search and filters identified active subnetworks
#' before returning final subnetworks.
#'
#' @param significant_genes vector of significant genes for the experiment
#' @inheritParams filter_active_subnetworks
#' @param search_method algorithm to use when performing active subnetwork
#'  search. Options are greedy search (GR), simulated annealing (SA) or genetic
#'  algorithm (GA) for the search (default = 'GR').
#' @param seed_for_stochastic_methods seed for reproducibility while running active subnetwork search
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
#' @param network Prebuilt network object as returned by \code{build_network()},
#'   built once by \code{active_snw_enrichment_wrapper} and passed to every
#'   iteration to avoid redundant PIN file I/O.
#' @param score_context Prebuilt score context as returned by \code{build_score_context()},
#'   built once by \code{active_snw_enrichment_wrapper} and passed to every
#'   iteration to avoid redundant Monte-Carlo calibration.
#'
#' @return A list of genes in every identified active subnetwork that has a score greater than
#' the `score_quan_thr`th quantile and that has at least `sig_gene_thr` affected genes.
#'
#' @export
#'
#' @examples
#' \donttest{
#' experiment_df <- example_pathfindR_input[1:15, c(1, 3)]
#' colnames(experiment_df) <- c("gene", "pvalue")
#' pin_path <- return_pin_path("KEGG")
#' network <- build_network(pin_path)
#' score_context <- build_score_context(
#'   network,
#'   experiment_df,
#'   list(p_for_nonsignificant = 0.5, seed = 1234L)
#' )
#' GR_snws <- get_active_subnetworks(
#'   significant_genes = experiment_df$gene,
#'   network = network,
#'   score_context = score_context
#' )
#' }
get_active_subnetworks <- function(
  significant_genes,
  network,
  score_context,
  score_quan_thr = 0.8,
  sig_gene_thr = 0.02,
  search_method = "GR",
  seed_for_stochastic_methods = 1234,
  verbose = FALSE,
  start_with_all_positives = FALSE,
  gene_init_prob = 0.1,
  sa_initial_temp = 1,
  sa_final_temp = 0.01,
  sa_iterations = 10000,
  ga_population_size = 400,
  ga_iterations = 200,
  ga_crossover_rate = 1,
  ga_mutation_rate = 0,
  gr_max_depth = 1,
  gr_search_depth = 1,
  gr_overlap_threshold = 0.5,
  gr_subnetwork_num = 1000
) {
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

  ############ Run active Subnetwork Search
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
    seed                     = seed_for_stochastic_methods
  )
  identified_active_snws <- active_subnetwork_search(
    network = network,
    score_context = score_context,
    method = search_method,
    params = params,
    verbose = verbose
  )

  ############ Parse and filter active subnetworks
  filtered_snws <- filter_active_subnetworks(
    active_snws = identified_active_snws, sig_genes_vec = significant_genes,
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
