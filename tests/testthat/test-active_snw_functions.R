##################################################
## Project: pathfindR
## Script purpose: Testthat testing script for
## active snw search functions
## Date: Aug 30, 2019
## Author: Ege Ulgen
##################################################

# filterActiveSnws --------------------------------------------------------
sample_path <- system.file("extdata/resultActiveSubnetworkSearch.txt",
                           package = "pathfindR")

test_that("Filter function returns list object", {

  tmp_filtered <- filterActiveSnws(sample_path, RA_input$Gene.symbol)
  expect_is(tmp_filtered, "list")
  expect_is(tmp_filtered[[1]], "character")

  # empty file case
  expect_null(suppressWarnings(filterActiveSnws("", RA_input$Gene.symbol)))
})

test_that("Filter function arguments work OK", {

  # supplied num. of genes
  tmp_filtered1 <- filterActiveSnws(sample_path, RA_input$Gene.symbol)
  tmp_filtered2 <- filterActiveSnws(sample_path, RA_input$Gene.symbol[1:100])
  expect_true(length(tmp_filtered1) > length(tmp_filtered2))


  # score_quan_thr
  tmp_filtered <- filterActiveSnws(active_snw_path = sample_path,
                                    signif_genes = RA_input$Gene.symbol,
                                    score_quan_thr = 0.8) # default
  expect_true(length(tmp_filtered1) <= 200)

  tmp_filtered <- filterActiveSnws(active_snw_path = sample_path,
                                    signif_genes = RA_input$Gene.symbol,
                                    score_quan_thr = 0.9)
  expect_true(length(tmp_filtered) <= 100)

  # sig_gene_thr
  tmp_filtered <- filterActiveSnws(active_snw_path = sample_path,
                                   signif_genes = RA_input$Gene.symbol,
                                   sig_gene_thr = 0)
  expect_true(length(tmp_filtered) == 200)

  tmp_filtered1 <- filterActiveSnws(active_snw_path = sample_path,
                                   signif_genes = RA_input$Gene.symbol,
                                   sig_gene_thr = 10) # default

  tmp_filtered2 <- filterActiveSnws(active_snw_path = sample_path,
                                    signif_genes = RA_input$Gene.symbol,
                                    sig_gene_thr = 20)

  expect_true(length(tmp_filtered1) <= 200)
  expect_true(length(tmp_filtered1) > length(tmp_filtered2))
})

# active_snw_search -------------------------------------------------------
pin_path <- return_pin_path()
input_df1 <- suppressMessages(input_processing(RA_input[1:50, ],
                                               p_val_threshold = 0.05,
                                               pin_path))
input_df2 <- suppressMessages(input_processing(RA_input[1:3, ],
                                               p_val_threshold = 0.05,
                                               pin_path))

test_that("Active snw search function returns list object", {
  # Expect > 0 active snws
  expect_message(snw_list <- active_snw_search(input_df1, pin_path),
                 "Found \\d+ active subnetworks")
  expect_is(snw_list, "list")
  expect_is(snw_list[[1]], "character")

  # Expect no active snws
  expect_message(snw_list <- active_snw_search(input_df2, pin_path),
                 "Found 0 active subnetworks")
  expect_is(snw_list, "list")
  expect_equal(length(snw_list), 0)

  # dir_for_parallel_run works?
  expect_message(snw_list <- active_snw_search(input_df2,
                                               pin_path,
                                               dir_for_parallel_run = "."),
                 "Found 0 active subnetworks")

})

test_that("Active snw search function error messages work", {
  expect_error(active_snw_search(input_df2,
                                 pin_path,
                                 search_method = "WRONG"),
               "`search_method` must be one of \"GR\", \"SA\", \"GA\"")

  expect_error(active_snw_search(input_df2,
                                 pin_path,
                                 use_all_positives = "WRONG"),
               "the argument `use_all_positives` must be either TRUE or FALSE")

  expect_error(active_snw_search(input_df2,
                                 pin_path,
                                 silent_option = "WRONG"),
               "the argument `silent_option` must be either TRUE or FALSE")
})



