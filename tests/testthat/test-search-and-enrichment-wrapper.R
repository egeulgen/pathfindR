set.seed(42)

test_that("`active_snw_enrichment_wrapper()` -- wraps and combines each output as expected", {
  mock_input <- example_pathfindR_input[1:10, c(1, 3)]
  colnames(mock_input) <- c("GENE", "P_VALUE")

  mock_pin_path <- "/path/to/pin"

  mock_output <- data.frame(A = c(1, 2), B = c(3, 4))

  with_mocked_bindings(
    {
      expect_is(res_1_iter <- pathfindR:::active_snw_enrichment_wrapper(
        input_processed = mock_input, pin_path = mock_pin_path,
        gset_list = list(), enrichment_threshold = 0.05, list_active_snw_genes = FALSE,
        iterations = 1
      ), "data.frame")
      expect_identical(res_1_iter, mock_output)

      expect_is(res_3iters <- pathfindR:::active_snw_enrichment_wrapper(
        input_processed = mock_input, pin_path = mock_pin_path,
        gset_list = list(), enrichment_threshold = 0.05, list_active_snw_genes = FALSE,
        iterations = 3, disable_parallel = TRUE
      ), "data.frame")
      expect_identical(res_3iters, rbind(mock_output, mock_output, mock_output))
    },
    single_iter_wrapper = function(...) mock_output,
    build_network = function(...) list(),
    .package = "pathfindR"
  )
})

test_that("`active_snw_enrichment_wrapper()` -- if using GA method, warning is raised and iterations set to 1", {
  mock_input <- example_pathfindR_input[1:10, c(1, 3)]
  colnames(mock_input) <- c("GENE", "P_VALUE")

  mock_pin_path <- "/path/to/pin"

  mock_output <- data.frame(A = c(1, 2), B = c(3, 4))

  with_mocked_bindings(
    {
      expect_warning(
        res <- pathfindR:::active_snw_enrichment_wrapper(
          input_processed = mock_input, pin_path = mock_pin_path, search_method = "GA",
          gset_list = list(), enrichment_threshold = 0.05, list_active_snw_genes = FALSE,
          iterations = 123
        )
      )
      expect_identical(res, mock_output)
    },
    single_iter_wrapper = function(...) mock_output,
    build_network = function(...) list(),
    .package = "pathfindR"
  )
})

test_that("`active_snw_enrichment_wrapper()` -- produces same result when run sequentially and in parallel", {
  ## Given
  # Input data - experiment
  n_input_genes <- 10
  input_genes <- paste0("GENE", seq_len(n_input_genes))
  input_p_vals <- runif(n_input_genes, min = 1e-10, max = 0.001)
  toy_input_df <- data.frame(GENE = input_genes, P_VALUE = input_p_vals)

  # Input data - interaction network
  pool <- paste0("GENE", 1:25)
  n_edges <- 50
  toy_pin_df <- data.frame(
    InteractorA = sample(pool, n_edges, replace = TRUE),
    pp = "pp",
    InteractorB = sample(pool, n_edges, replace = TRUE),
    stringsAsFactors = FALSE
  )
  # remove self-loops
  toy_pin_df <- subset(toy_pin_df, InteractorA != InteractorB)
  # remove duplicate edges
  toy_pin_df <- toy_pin_df[
    !duplicated(
      t(apply(toy_pin_df[c("InteractorA", "InteractorB")], 1, sort))
    ),
  ]
  sif_file <- tempfile(fileext = ".sif")
  write.table(
    toy_pin_df,
    sif_file,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )
  on.exit(unlink(sif_file), add = TRUE)

  # input data - gene sets
  mock_genes_by_term <- list(A = c("GENE5", "GENE3", "GENE7", "GENE1", "GENE4"), B = c("GENE6", "GENE9", "GENE17", "GENE3"))
  mock_term_descriptions <- c(A = "gene set A", B = "genes set B")
  mock_gset_list <- list(genes_by_term = mock_genes_by_term, term_descriptions = mock_term_descriptions)

  ## When
  parallel_res <- active_snw_enrichment_wrapper(
    input_processed = toy_input_df,
    pin_path = sif_file,
    gset_list = list(genes_by_term = mock_genes_by_term, term_descriptions = mock_term_descriptions),
    enrichment_threshold = 0.05,
    list_active_snw_genes = FALSE,
    iterations = 2,
    disable_parallel = FALSE
  )
  sequential_res <- active_snw_enrichment_wrapper(
    input_processed = toy_input_df,
    pin_path = sif_file,
    gset_list = mock_gset_list,
    enrichment_threshold = 0.05,
    list_active_snw_genes = FALSE,
    iterations = 2,
    disable_parallel = TRUE
  )

  ## Then
  expect_identical(parallel_res, sequential_res)
})

test_that("`single_iter_wrapper()` produces same results when same seed is set", {
  set.seed(42)
  ## Given
  pool <- paste0("GENE", 1:250)
  n_input_genes <- 75
  input_genes <- paste0("GENE", seq_len(n_input_genes))

  prop_alt <- 0.3
  is_alt <- rbinom(n_input_genes, 1, prop_alt)
  signal_genes <- input_genes[is_alt == 1]

  # near-zero p-values for signal genes -> strong, clearly separable signal
  input_p_vals <- ifelse(
    is_alt == 1,
    runif(n_input_genes, min = 0, max = 1e-4),
    runif(n_input_genes)
  )
  toy_input_df <- data.frame(GENE = input_genes, P_VALUE = input_p_vals)

  # PIN: a dense clique among the signal genes guarantees a discoverable
  # connected module, plus random background edges as noise.
  seed_edges <- as.data.frame(t(combn(signal_genes, 2)))
  names(seed_edges) <- c("InteractorA", "InteractorB")
  seed_edges$pp <- "pp"
  seed_edges <- seed_edges[c("InteractorA", "pp", "InteractorB")]

  n_background_edges <- 800
  background_edges <- data.frame(
    InteractorA = sample(pool, n_background_edges, replace = TRUE),
    pp = "pp",
    InteractorB = sample(pool, n_background_edges, replace = TRUE),
    stringsAsFactors = FALSE
  )

  toy_pin_df <- rbind(seed_edges, background_edges)
  toy_pin_df <- subset(toy_pin_df, InteractorA != InteractorB)
  toy_pin_df <- toy_pin_df[
    !duplicated(t(apply(toy_pin_df[c("InteractorA", "InteractorB")], 1, sort))),
  ]

  sif_file <- tempfile(fileext = ".sif")
  write.table(toy_pin_df, sif_file,
    sep = "\t",
    row.names = FALSE, col.names = FALSE, quote = FALSE
  )
  on.exit(unlink(sif_file), add = TRUE)

  mock_gset_list <- list(
    genes_by_term = list(
      A = c("GENE5", "GENE3", "GENE7", "GENE1", "GENE4"),
      B = c("GENE6", "GENE9", "GENE17", "GENE3")
    ),
    term_descriptions = c(A = "gene set A", B = "genes set B")
  )

  network <- build_network(sif_file)
  toy_experiment_df <- data.frame(
    gene = toupper(toy_input_df$GENE),
    pvalue = toy_input_df$P_VALUE
  )

  n_iter <- 3
  gene_init_prob_vec <- rep(0.4, n_iter)

  run_once <- function(i, search_method) {
    pathfindR:::single_iter_wrapper(
      i = i,
      pin_path = sif_file,
      network = network,
      experiment_df = toy_experiment_df,
      score_quan_thr = -1,
      sig_gene_thr = 0,
      search_method = search_method,
      verbose = FALSE,
      start_with_all_positives = FALSE,
      gene_init_prob = gene_init_prob_vec,
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
      gset_list = mock_gset_list,
      adj_method = "bonferroni",
      enrichment_threshold = 0.5,
      list_active_snw_genes = FALSE
    )
  }

  ## When / Then
  for (search_method in c("GR", "SA")) {
    iter_idx <- c(1, 2, 1)
    results <- lapply(iter_idx, run_once, search_method = search_method)

    expect_false(
      identical(results[[1]], results[[2]]),
      info = sprintf("[%s] different seeds unexpectedly gave identical results", search_method)
    )
    expect_identical(
      results[[1]], results[[3]],
      info = sprintf("[%s] same seed did not reproduce identical results", search_method)
    )
  }
})

test_that("`active_snw_enrichment_wrapper()` -- argument checks work", {
  valid_mets <- c("GR", "SA", "GA")
  mock_pin_path <- "/path/to/pin"

  expect_error(active_snw_enrichment_wrapper(
    input_processed = mock_input,
    pin_path = mock_pin_path, gset_list = list(), enrichment_threshold = 0.05, list_active_snw_genes = FALSE,
    search_method = "INVALID"
  ), paste0("`search_method` should be one of ", paste(dQuote(valid_mets),
    collapse = ", "
  )))

  expect_error(active_snw_enrichment_wrapper(
    input_processed = mock_input,
    pin_path = mock_pin_path, gset_list = list(), enrichment_threshold = 0.05, list_active_snw_genes = FALSE,
    start_with_all_positives = "INVALID"
  ), "`start_with_all_positives` should be either TRUE or FALSE")

  expect_error(active_snw_enrichment_wrapper(
    input_processed = mock_input,
    pin_path = mock_pin_path, gset_list = list(), enrichment_threshold = 0.05, list_active_snw_genes = FALSE,
    verbose = "INVALID"
  ), "`verbose` should be either TRUE or FALSE")

  expect_error(active_snw_enrichment_wrapper(
    input_processed = mock_input,
    pin_path = mock_pin_path, gset_list = list(), enrichment_threshold = 0.05, list_active_snw_genes = FALSE,
    disable_parallel = "INVALID"
  ), "`disable_parallel` should be either TRUE or FALSE")

  expect_error(active_snw_enrichment_wrapper(
    input_processed = mock_input,
    pin_path = mock_pin_path, gset_list = list(), enrichment_threshold = 0.05, list_active_snw_genes = FALSE,
    iterations = "INVALID"
  ), "`iterations` should be a positive integer")

  expect_error(active_snw_enrichment_wrapper(
    input_processed = mock_input,
    pin_path = mock_pin_path, gset_list = list(), enrichment_threshold = 0.05, list_active_snw_genes = FALSE,
    iterations = 0
  ), "`iterations` should be >= 1")

  expect_error(active_snw_enrichment_wrapper(
    input_processed = mock_input,
    pin_path = mock_pin_path, gset_list = list(), enrichment_threshold = 0.05, list_active_snw_genes = FALSE,
    n_processes = "INVALID"
  ), "`n_processes` should be either NULL or a positive integer")

  expect_error(active_snw_enrichment_wrapper(
    input_processed = mock_input,
    pin_path = mock_pin_path, gset_list = list(), enrichment_threshold = 0.05, list_active_snw_genes = FALSE,
    n_processes = 0
  ), "`n_processes` should be >= 1")
})
