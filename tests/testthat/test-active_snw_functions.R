##################################################
## Project: pathfindR
## Script purpose: Testthat testing script for
## active subnetwork search functions
## Date: Oct 20, 2019
## Author: Ege Ulgen
##################################################

# filterActiveSnws --------------------------------------------------------
sample_path <- system.file("extdata/resultActiveSubnetworkSearch.txt",
                           package = "pathfindR")
def_smw_len <- 1000

test_that("`filterActiveSnws()` returns list object", {
  tmp_filtered <- filterActiveSnws(active_snw_path = sample_path,
                                   sig_genes_vec = RA_input$Gene.symbol)
  expect_is(tmp_filtered, "list")
  expect_is(tmp_filtered[[1]], "character")

  # empty file case
  empty_path <- file.path(tempdir(), "empty.txt")
  file.create(empty_path)
  expect_null(suppressWarnings(filterActiveSnws(active_snw_path = empty_path,
                                                sig_genes_vec = RA_input$Gene.symbol)))
})

test_that("`sig_genes_vec` in `filterActiveSnws()` works", {
  tmp_filtered1 <- filterActiveSnws(active_snw_path = sample_path,
                                    sig_genes_vec = RA_input$Gene.symbol)
  tmp_filtered2 <- filterActiveSnws(active_snw_path = sample_path,
                                    sig_genes_vec = RA_input$Gene.symbol[1:100])
  expect_true(length(tmp_filtered1) > length(tmp_filtered2))
})

test_that("`score_quan_thr` in `filterActiveSnws()` works", {
  tmp_filtered <- filterActiveSnws(active_snw_path = sample_path,
                                   sig_genes_vec = RA_input$Gene.symbol,
                                   score_quan_thr = -1,
                                   sig_gene_thr = 0)
  expect_length(tmp_filtered, def_smw_len)

  for (q_thr in seq(.1, 1, by = .1)) {
    tmp_filtered <- filterActiveSnws(active_snw_path = sample_path,
                                     sig_genes_vec = RA_input$Gene.symbol,
                                     score_quan_thr = q_thr,
                                     sig_gene_thr = 0)
    exp_len <- def_smw_len * (1 - q_thr)
    expect_length(tmp_filtered, as.integer(exp_len + .5))
  }
})

test_that("`sig_gene_thr` in `filterActiveSnws()` works", {
  tmp_filtered1 <- filterActiveSnws(active_snw_path = sample_path,
                                    sig_genes_vec = RA_input$Gene.symbol,
                                    sig_gene_thr = 10, # default
                                    score_quan_thr = -1)

  tmp_filtered2 <- filterActiveSnws(active_snw_path = sample_path,
                                    sig_genes_vec = RA_input$Gene.symbol,
                                    sig_gene_thr = 20,
                                    score_quan_thr = -1)

  expect_true(length(tmp_filtered2) < def_smw_len)
  expect_true(length(tmp_filtered1) > length(tmp_filtered2))
})

test_that("`filterActiveSnws()` arg checks work", {
  expect_error(filterActiveSnws(active_snw_path = "this/is/not/a/path"),
               "The active subnetwork file does not exist! Check the `active_snw_path` argument")

  expect_error(filterActiveSnws(active_snw_path = sample_path,
                                sig_genes_vec = 1:3),
               "`sig_genes_vec` must be a character vector")

  expect_error(filterActiveSnws(active_snw_path = sample_path,
                                sig_genes_vec = RA_input$Gene.symbol,
                                score_quan_thr = "NOT_NUMERIC"),
               "`score_quan_thr` must be numeric!")

  expect_error(filterActiveSnws(active_snw_path = sample_path,
                                sig_genes_vec = RA_input$Gene.symbol,
                                score_quan_thr = -2),
               "`score_quan_thr` must be in \\[0, 1\\] or -1 \\(if not filtering\\)")

  expect_error(filterActiveSnws(active_snw_path = sample_path,
                                sig_genes_vec = RA_input$Gene.symbol,
                                score_quan_thr = 2),
               "`score_quan_thr` must be in \\[0, 1\\] or -1 \\(if not filtering\\)")
})

# active_snw_search -------------------------------------------------------
pin_path <- return_pin_path()
input_df1 <- suppressMessages(input_processing(RA_input[1:100, ],
                                               p_val_threshold = 0.05,
                                               pin_path))
input_df2 <- suppressMessages(input_processing(RA_input[1:3, ],
                                               p_val_threshold = 0.05,
                                               pin_path))

test_that("`active_snw_search()` returns list object", {
  # Expect > 0 active snws
  expect_message(snw_list <- active_snw_search(input_for_search = input_df1,
                                               pin_path = pin_path),
                 "Found [1-9]\\d* active subnetworks")
  expect_is(snw_list, "list")
  expect_is(snw_list[[1]], "character")
  unlink("active_snw_search", recursive = TRUE)

  # Expect no active snws
  expect_message(snw_list <- active_snw_search(input_for_search = input_df2,
                                               pin_path = pin_path),
                 "Found 0 active subnetworks")
  expect_identical(snw_list, list())
  unlink("active_snw_search", recursive = TRUE)

  # dir_for_parallel_run works?
  dir.create("dummy_dir")
  expect_message(snw_list <- active_snw_search(input_for_search = input_df1,
                                               pin_path = pin_path,
                                               dir_for_parallel_run = "dummy_dir"),
                 "Found [1-9]\\d* active subnetworks")
  expect_true(file.exists("dummy_dir/active_snw_search/active_snws.txt"))
  unlink("dummy_dir", recursive = TRUE)
})

test_that("All search methods for `active_snw_search()` work", {
  ## GR
  expect_message(snw_list <- active_snw_search(input_for_search = input_df1,
                                               pin_path = pin_path,
                                               search_method = "GR"),
                 "Found [1-9]\\d* active subnetworks")
  expect_is(snw_list, "list")
  expect_is(snw_list[[1]], "character")
  unlink("active_snw_search", recursive = TRUE)

  skip("will test SA and GA if we can create a suitable (faster and non-empty) test case")
  ## SA
  expect_message(snw_list <- active_snw_search(input_for_search = input_df1[1:100, ],
                                               pin_path = pin_path,
                                               search_method = "SA"),
                 "Found [1-9]\\d* active subnetworks")
  expect_is(snw_list, "list")
  expect_is(snw_list[[1]], "character")
  unlink("active_snw_search", recursive = TRUE)

  ## GA
  expect_message(snw_list <- active_snw_search(input_for_search = input_df1,
                                               pin_path = pin_path,
                                               search_method = "GA"),
                 "Found [1-9]\\d* active subnetworks")
  expect_is(snw_list, "list")
  expect_is(snw_list[[1]], "character")
  unlink("active_snw_search", recursive = TRUE)
})

test_that("`active_snw_search()` error messages work", {
  expect_error(active_snw_search(input_for_search = input_df2,
                                 pin_path = pin_path,
                                 search_method = "WRONG"),
    "`search_method` must be one of \"GR\", \"SA\", \"GA\"")

  expect_error(active_snw_search(input_for_search = input_df2,
                                 pin_path = pin_path,
                                 use_all_positives = "WRONG"),
    "the argument `use_all_positives` must be either TRUE or FALSE")

  expect_error(active_snw_search(input_for_search = input_df2,
                                 pin_path = pin_path,
                                 silent_option = "WRONG"),
               "the argument `silent_option` must be either TRUE or FALSE")
})
