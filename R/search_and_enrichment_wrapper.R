#' Wrapper for Active Subnetwork Search + Enrichment over Single/Multiple Iteration(s)
#'
#' @param input_processed processed input data frame
#' @param pin_path path/to/PIN/file
#' @param gset_list list for gene sets.
#' @param disable_parallel boolean to indicate whether to disable parallel runs
#'  via \code{foreach} (default = FALSE)
#' @inheritParams run_pathfindR
#' @inheritParams active_snw_search
#' @inheritParams enrichment_analyses
#' @param iterations number of iterations for active subnetwork search and
#'  enrichment analyses (Default = 10)
#' @param n_processes optional argument for specifying the number of processes
#'  used by foreach. If not specified, the function determines this
#'  automatically (Default == NULL. Gets set to 1 for Genetic Algorithm)
#'
#' @return Data frame of combined pathfindR enrichment results
active_snw_enrichment_wrapper <- function(input_processed, pin_path, gset_list, enrichment_threshold,
                                          list_active_snw_genes, adj_method = "bonferroni", search_method = "GR", disable_parallel = FALSE,
                                          use_all_positives = FALSE, iterations = 10, n_processes = NULL, score_quan_thr = 0.8,
                                          sig_gene_thr = 0.02, saTemp0 = 1, saTemp1 = 0.01, saIter = 10000, gaPop = 400,
                                          gaIter = 200, gaThread = 5, gaCrossover = 1, gaMut = 0, grMaxDepth = 1, grSearchDepth = 1,
                                          grOverlap = 0.5, grSubNum = 1000, silent_option = TRUE) {
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
      stop("`n_processes` should be > 1")
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

  if (!is.logical(use_all_positives)) {
    stop("`use_all_positives` should be either TRUE or FALSE")
  }

  if (!is.logical(silent_option)) {
    stop("`silent_option` should be either TRUE or FALSE")
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

  geneInitProbs <- 0.1
  dirs <- c()
  if (iterations > 1) {
    geneInitProbs <- seq(from = 0.01, to = 0.2, length.out = iterations)

    for (i in base::seq_len(iterations)) {
      dir_i <- file.path("active_snw_searches", paste0("Iteration_", i))
      dir.create(dir_i, recursive = TRUE, showWarnings = FALSE)
      dirs <- c(dirs, dir_i)
    }
  }

  if (iterations == 1) {
    combined_res <- single_iter_wrapper(i = NULL, dirs, input_processed, pin_path,
                                        score_quan_thr, sig_gene_thr, search_method, silent_option, use_all_positives,
                                        geneInitProbs, saTemp0, saTemp1, saIter, gaPop, gaIter, gaThread, gaCrossover,
                                        gaMut, grMaxDepth, grSearchDepth, grOverlap, grSubNum, gset_list, adj_method,
                                        enrichment_threshold, list_active_snw_genes)
  } else {
    if (!disable_parallel) {
      cl <- parallel::makeCluster(n_processes, setup_strategy = "sequential")
      doParallel::registerDoParallel(cl)
      `%dopar%` <- foreach::`%dopar%`
      combined_res <- foreach::foreach(i = 1:iterations, .combine = rbind,
                                       .packages = "pathfindR") %dopar% {
                                         single_iter_wrapper(i, dirs, input_processed, pin_path, score_quan_thr,
                                                             sig_gene_thr, search_method, silent_option, use_all_positives,
                                                             geneInitProbs, saTemp0, saTemp1, saIter, gaPop, gaIter, gaThread,
                                                             gaCrossover, gaMut, grMaxDepth, grSearchDepth, grOverlap, grSubNum,
                                                             gset_list, adj_method, enrichment_threshold, list_active_snw_genes)
                                       }
      parallel::stopCluster(cl)
    } else {
      combined_res <- c()
      for (i in 1:iterations) {
        current_res <- single_iter_wrapper(i, dirs, score_quan_thr, sig_gene_thr,
                                           search_method, silent_option, use_all_positives, geneInitProbs,
                                           saTemp0, saTemp1, saIter, gaPop, gaIter, gaThread, gaCrossover,
                                           gaMut, grMaxDepth, grSearchDepth, grOverlap, grSubNum, gset_list,
                                           adj_method, enrichment_threshold, list_active_snw_genes)
        combined_res <- rbind(combined_res, current_res)
      }
    }
  }
  return(combined_res)
}

#' Active Subnetwork Search + Enrichment Analysis Wrapper for a Single Iteration
#'
#' @param i current iteration index (default = \code{NULL})
#' @param dirs vector of directories for parallel runs
#' @inheritParams active_snw_search
#' @inheritParams enrichment_analyses
#' @inheritParams active_snw_enrichment_wrapper
#'
#' @return Data frame of enrichment results using active subnetwork search results
single_iter_wrapper <- function(i = NULL, dirs, input_processed, pin_path, score_quan_thr,
                                sig_gene_thr, search_method, silent_option, use_all_positives, geneInitProbs,
                                saTemp0, saTemp1, saIter, gaPop, gaIter, gaThread, gaCrossover, gaMut, grMaxDepth,
                                grSearchDepth, grOverlap, grSubNum, gset_list, adj_method, enrichment_threshold,
                                list_active_snw_genes) {
  snws_file <- "active_snws"
  dir_for_parallel_run <- NULL
  if (!is.null(i)) {
    snws_file <- paste0("active_snws_", i)
    dir_for_parallel_run <- dirs[i]
  }
  snws <- active_snw_search(input_for_search = input_processed, pin_name_path = pin_path,
                            snws_file = snws_file, dir_for_parallel_run = dir_for_parallel_run, score_quan_thr = score_quan_thr,
                            sig_gene_thr = sig_gene_thr, search_method = search_method, seedForRandom = ifelse(is.null(i),
                                                                                                               1234, i), silent_option = silent_option, use_all_positives = use_all_positives,
                            geneInitProbs = ifelse(!is.null(i), geneInitProbs[i], geneInitProbs), saTemp0 = saTemp0,
                            saTemp1 = saTemp1, saIter = saIter, gaPop = gaPop, gaIter = gaIter, gaThread = gaThread,
                            gaCrossover = gaCrossover, gaMut = gaMut, grMaxDepth = grMaxDepth, grSearchDepth = grSearchDepth,
                            grOverlap = grOverlap, grSubNum = grSubNum)
  enrichment_res <- enrichment_analyses(snws = snws, sig_genes_vec = input_processed$GENE,
                                        pin_name_path = pin_path, genes_by_term = gset_list$genes_by_term, term_descriptions = gset_list$term_descriptions,
                                        adj_method = adj_method, enrichment_threshold = enrichment_threshold, list_active_snw_genes = list_active_snw_genes)
  return(enrichment_res)
}
