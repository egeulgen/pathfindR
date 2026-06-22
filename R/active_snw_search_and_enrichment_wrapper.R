#' Wrapper for Active Subnetwork Search + Enrichment over Single/Multiple Iteration(s)
#'
#' @param input_processed processed input data frame
#' @param pin_path path/to/PIN/file
#' @param gset_list list for gene sets.
#' @param disable_parallel boolean to indicate whether to disable parallel runs
#'  via \code{foreach} (default = FALSE)
#' @inheritParams run_pathfindR
#' @inheritParams get_active_subnetworks
#' @inheritParams enrichment_analyses
#' @param iterations number of iterations for active subnetwork search and
#'  enrichment analyses (Default = 10)
#' @param n_processes optional argument for specifying the number of processes
#'  used by foreach. If not specified, the function determines this
#'  automatically (Default == NULL. Gets set to 1 for Genetic Algorithm)
#'
#' @return Data frame of combined pathfindR enrichment results
active_snw_enrichment_wrapper <- function(
  input_processed,
  pin_path,
  gset_list,
  enrichment_threshold,
  list_active_snw_genes,
  adj_method = "bonferroni",
  search_method = "GR",
  disable_parallel = FALSE,
  start_with_all_positives = FALSE,
  iterations = 10,
  n_processes = NULL,
  score_quan_thr = 0.8,
  sig_gene_thr = 0.02,
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
  gr_subnetwork_num = 1000,
  verbose = FALSE
) {
  message("## Performing Active Subnetwork Search and Enrichment")
  ############ Argument checks Active Subnetwork Search Method
  valid_mets <- c("GR", "SA", "GA")
  if (!search_method %in% valid_mets) {
    stop("`search_method` should be one of ", paste(dQuote(valid_mets), collapse = ", "))
  }

  ## If search_method is GA, set iterations as 1
  if (search_method == "GA") {
    warning("`iterations` is set to 1 because `search_method = \"GA\"`", call. = FALSE)
    iterations <- 1
  }

  if (!is.null(n_processes)) {
    if (!is.numeric(n_processes)) {
      stop("`n_processes` should be either NULL or a positive integer")
    }
    if (n_processes < 1) {
      stop("`n_processes` should be >= 1")
    }
  }

  # calculate the number of processes, if necessary
  if (is.null(n_processes)) {
    n_processes <- parallel::detectCores() - 1
  }

  ## If iterations < n_processes, set n_processes to iterations
  if (iterations < n_processes & iterations != 1) {
    message("`n_processes` is set to `iterations` because `iterations` < `n_processes`")
    n_processes <- iterations
  }

  if (!is.logical(start_with_all_positives)) {
    stop("`start_with_all_positives` should be either TRUE or FALSE")
  }

  if (!is.logical(verbose)) {
    stop("`verbose` should be either TRUE or FALSE")
  }

  if (!is.logical(disable_parallel)) {
    stop("`disable_parallel` should be either TRUE or FALSE")
  }

  if (!is.numeric(iterations)) {
    stop("`iterations` should be a positive integer")
  }
  if (iterations < 1) {
    stop("`iterations` should be >= 1")
  }

  gene_init_prob <- 0.1
  if (iterations > 1) {
    gene_init_prob <- seq(from = 0.01, to = 0.2, length.out = iterations)
  }

  ## Build the network once, before any iteration starts.  The network depends
  ## only on the PIN file, so it is genuinely identical across iterations and the
  ## (expensive) SIF parse + Java-order reconstruction is done a single time.
  message("## Building network (once for all iterations)")
  network <- build_network(pin_path)

  ## Seed-independent experiment data frame
  ## built once and passed to every iteration.
  experiment_df <- data.frame(
    gene = base::toupper(input_processed$GENE),
    pvalue = input_processed$P_VALUE
  )

  if (iterations == 1) {
    combined_res <- single_iter_wrapper(
      i = NULL,
      pin_path,
      network,
      experiment_df,
      score_quan_thr,
      sig_gene_thr,
      search_method,
      verbose,
      start_with_all_positives,
      gene_init_prob,
      sa_initial_temp,
      sa_final_temp,
      sa_iterations,
      ga_population_size,
      ga_iterations,
      ga_crossover_rate,
      ga_mutation_rate,
      gr_max_depth,
      gr_search_depth,
      gr_overlap_threshold,
      gr_subnetwork_num,
      gset_list,
      adj_method,
      enrichment_threshold,
      list_active_snw_genes
    )
  } else {
    if (!disable_parallel) {
      cl <- parallel::makeCluster(n_processes, setup_strategy = "sequential")
      doParallel::registerDoParallel(cl)
      `%dopar%` <- foreach::`%dopar%`
      combined_res <- foreach::foreach(
        i = 1:iterations, .combine = rbind,
        .packages = "pathfindR"
      ) %dopar% {
        single_iter_wrapper(
          i,
          pin_path,
          network,
          experiment_df,
          score_quan_thr,
          sig_gene_thr,
          search_method,
          verbose,
          start_with_all_positives,
          gene_init_prob,
          sa_initial_temp,
          sa_final_temp,
          sa_iterations,
          ga_population_size,
          ga_iterations,
          ga_crossover_rate,
          ga_mutation_rate,
          gr_max_depth,
          gr_search_depth,
          gr_overlap_threshold,
          gr_subnetwork_num,
          gset_list,
          adj_method,
          enrichment_threshold,
          list_active_snw_genes
        )
      }
      parallel::stopCluster(cl)
    } else {
      combined_res <- c()
      for (i in 1:iterations) {
        current_res <- single_iter_wrapper(
          i,
          pin_path,
          network,
          experiment_df,
          score_quan_thr,
          sig_gene_thr,
          search_method,
          verbose,
          start_with_all_positives,
          gene_init_prob,
          sa_initial_temp,
          sa_final_temp,
          sa_iterations,
          ga_population_size,
          ga_iterations,
          ga_crossover_rate,
          ga_mutation_rate,
          gr_max_depth,
          gr_search_depth,
          gr_overlap_threshold,
          gr_subnetwork_num,
          gset_list,
          adj_method,
          enrichment_threshold,
          list_active_snw_genes
        )
        combined_res <- rbind(combined_res, current_res)
      }
    }
  }
  return(combined_res)
}

#' Active Subnetwork Search + Enrichment Analysis Wrapper for a Single Iteration
#'
#' @param i current iteration index (default = \code{NULL})
#' @param experiment_df input experiment data frame
#' @inheritParams get_active_subnetworks
#' @inheritParams enrichment_analyses
#' @inheritParams active_snw_enrichment_wrapper
#'
#' @return Data frame of enrichment results using active subnetwork search results
single_iter_wrapper <- function(
  i = NULL,
  pin_path,
  network,
  experiment_df,
  score_quan_thr,
  sig_gene_thr,
  search_method,
  verbose,
  start_with_all_positives,
  gene_init_prob,
  sa_initial_temp,
  sa_final_temp,
  sa_iterations,
  ga_population_size,
  ga_iterations,
  ga_crossover_rate,
  ga_mutation_rate,
  gr_max_depth,
  gr_search_depth,
  gr_overlap_threshold,
  gr_subnetwork_num,
  gset_list,
  adj_method,
  enrichment_threshold,
  list_active_snw_genes
) {
  ## Per-iteration seed.  Matches the original pathfindR, which runs the search
  ## (and its Monte-Carlo calibration) with `-seedForRandom = i` on iteration i,
  ## and with 1234 for a single (i = NULL) run.
  iter_seed <- ifelse(is.null(i), 1234, i)

  ## Build this iteration's score context with that seed.  The z-scores are
  ## seed-independent; the means/stds (and hence the calibrated scores driving
  ## the score-quantile filter) depend on the seed, which is what gives each
  ## iteration a different — and collectively more diverse — set of subnetworks.
  score_context <- build_score_context(
    network, experiment_df,
    list(p_for_nonsignificant = 0.5, seed = iter_seed)
  )

  snws <- get_active_subnetworks(
    significant_genes = experiment_df$gene,
    network = network,
    score_context = score_context,
    score_quan_thr = score_quan_thr,
    sig_gene_thr = sig_gene_thr,
    search_method = search_method,
    seed_for_stochastic_methods = iter_seed,
    verbose = verbose,
    start_with_all_positives = start_with_all_positives,
    gene_init_prob = ifelse(!is.null(i), gene_init_prob[i], gene_init_prob),
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
    gr_subnetwork_num = gr_subnetwork_num
  )
  enrichment_res <- enrichment_analyses(
    snws = snws,
    sig_genes_vec = experiment_df$gene,
    pin_name_path = pin_path,
    genes_by_term = gset_list$genes_by_term,
    term_descriptions = gset_list$term_descriptions,
    adj_method = adj_method,
    enrichment_threshold = enrichment_threshold,
    list_active_snw_genes = list_active_snw_genes
  )
  return(enrichment_res)
}
