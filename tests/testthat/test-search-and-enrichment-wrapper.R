set.seed(123)

test_that("`active_snw_enrichment_wrapper()` -- works as expected", {
  input_df <- example_pathfindR_input[, c(1, 3)]
  colnames(input_df) <- c("GENE", "P_VALUE")

  org_dir <- getwd()
  test_directory <- file.path(tempdir(check = TRUE), "snw_wrapper_test")
  dir.create(test_directory)
  setwd(test_directory)
  on.exit(setwd(org_dir))
  on.exit(unlink(test_directory), add = TRUE)

  with_mocked_bindings(
    {
      expect_is(active_snw_enrichment_wrapper(
        input_processed = input_df, pin_path = "Biogrid",
        gset_list = list(), enrichment_threshold = 0.05, list_active_snw_genes = FALSE,
        iterations = 1
      ), "data.frame")

      expect_is(active_snw_enrichment_wrapper(
        input_processed = input_df, pin_path = "Biogrid",
        gset_list = list(), enrichment_threshold = 0.05, list_active_snw_genes = FALSE,
        iterations = 2, disable_parallel = TRUE
      ), "data.frame")

      expect_warning(active_snw_enrichment_wrapper(
        input_processed = input_df,
        pin_path = "Biogrid", gset_list = list(), enrichment_threshold = 0.05,
        list_active_snw_genes = FALSE, search_method = "GA", iterations = 2
      ))
    },
    single_iter_wrapper = function(...) example_pathfindR_output,
    .package = "pathfindR"
  )

  skip_on_cran()
  expect_is(
    active_snw_enrichment_wrapper(
      input_processed = input_df[1:10, ], pin_path = "Biogrid",
      gset_list = list(genes_by_term = kegg_genes[1:2], term_descriptions = kegg_descriptions[names(kegg_genes[1:2])]),
      enrichment_threshold = 0.05, list_active_snw_genes = FALSE, iterations = 2
    ),
    "NULL"
  )
})

test_that("`active_snw_enrichment_wrapper()` -- argument checks work", {
  valid_mets <- c("GR", "SA", "GA")
  expect_error(active_snw_enrichment_wrapper(
    input_processed = input_processed,
    pin_path = pin_path, gset_list = list(), enrichment_threshold = 0.05, list_active_snw_genes = FALSE,
    search_method = "INVALID"
  ), paste0("`search_method` should be one of ", paste(dQuote(valid_mets),
    collapse = ", "
  )))

  expect_error(active_snw_enrichment_wrapper(
    input_processed = input_processed,
    pin_path = pin_path, gset_list = list(), enrichment_threshold = 0.05, list_active_snw_genes = FALSE,
    use_all_positives = "INVALID"
  ), "`use_all_positives` should be either TRUE or FALSE")

  expect_error(active_snw_enrichment_wrapper(
    input_processed = input_processed,
    pin_path = pin_path, gset_list = list(), enrichment_threshold = 0.05, list_active_snw_genes = FALSE,
    silent_option = "INVALID"
  ), "`silent_option` should be either TRUE or FALSE")

  expect_error(active_snw_enrichment_wrapper(
    input_processed = input_processed,
    pin_path = pin_path, gset_list = list(), enrichment_threshold = 0.05, list_active_snw_genes = FALSE,
    disable_parallel = "INVALID"
  ), "`disable_parallel` should be either TRUE or FALSE")

  expect_error(active_snw_enrichment_wrapper(
    input_processed = input_processed,
    pin_path = pin_path, gset_list = list(), enrichment_threshold = 0.05, list_active_snw_genes = FALSE,
    iterations = "INVALID"
  ), "`iterations` should be a positive integer")

  expect_error(active_snw_enrichment_wrapper(
    input_processed = input_processed,
    pin_path = pin_path, gset_list = list(), enrichment_threshold = 0.05, list_active_snw_genes = FALSE,
    iterations = 0
  ), "`iterations` should be >= 1")

  expect_error(active_snw_enrichment_wrapper(
    input_processed = input_processed,
    pin_path = pin_path, gset_list = list(), enrichment_threshold = 0.05, list_active_snw_genes = FALSE,
    n_processes = "INVALID"
  ), "`n_processes` should be either NULL or a positive integer")

  expect_error(active_snw_enrichment_wrapper(
    input_processed = input_processed,
    pin_path = pin_path, gset_list = list(), enrichment_threshold = 0.05, list_active_snw_genes = FALSE,
    n_processes = 0
  ), "`n_processes` should be > 1")
})
