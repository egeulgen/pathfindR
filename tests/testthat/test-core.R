##################################################
## Package: pathfindR
## Script purpose: Unit testing script for
## core function
## Date: Apr 27, 2023
## Author: Ege Ulgen
##################################################

test_that("`run_pathfindR()` works as expected", {
  skip_on_cran()

  out_dir <- file.path(tempdir(check = TRUE), "pathfindR_results")
  ## GR
  expect_is(
    run_pathfindR(
      input = example_pathfindR_input[1:50, ],
      iterations = 1,
      score_quan_thr = 0.8
    ),
    "data.frame"
  )
  expect_true(dir.exists(out_dir))

  expect_is(
    run_pathfindR(
      input = example_pathfindR_input[1:50, ],
      iterations = 2,
      n_processes = 5,
      gene_sets = "BioCarta",
      pin_name_path = "GeneMania",
      score_quan_thr = 0.8,
      plot_enrichment_chart = FALSE,
      output_dir = out_dir
    ),
    "data.frame"
  )
  expect_true(dir.exists(paste0(out_dir, "(1)")))

  ## GA - n_processes <- 1 and n_processes <- iterations (iterations < n_processes)
  expected_warns <- c(
    "Did not find any enriched terms!",
    "`iterations` is set to 1 because `search_method = \"GA\""
  )
  expect_warning(
    run_pathfindR(
      input = example_pathfindR_input[3:4, ],
      search_method = "GA",
      iterations = 2,
      score_quan_thr = 0.8,
      pin_name_path = "KEGG",
      output_dir = file.path(tempdir(), "GA_example")
    ),
    paste0(paste(expected_warns, collapse = "|")),
    all = TRUE, perl = TRUE
  )

  skip("will test SA and GA if we can create a suitable (faster and non-empty) test case")

  ## SA
  expect_is(
    run_pathfindR(
      input = example_pathfindR_input[1:50, ],
      iterations = 1,
      gene_sets = "GO-BP",
      pin_name_path = "GeneMania",
      search_method = "SA",
      plot_enrichment_chart = FALSE,
      output_dir = out_dir
    ),
    "data.frame"
  )

  ## GA
  expect_is(
    run_pathfindR(
      input = example_pathfindR_input[1:50, ],
      gene_sets = "GO-BP",
      pin_name_path = "GeneMania",
      search_method = "GA",
      plot_enrichment_chart = FALSE,
      output_dir = out_dir
    ),
    "data.frame"
  )
})

test_that("Expect warning with empty result from `run_pathfindR()`", {
  expect_warning(
    res <- run_pathfindR(
      input = example_pathfindR_input[1:2, ],
      iterations = 1,
      output_dir = file.path(tempdir(check = TRUE), "pathfindR_results")
    ),
    "Did not find any enriched terms!"
  )
  expect_identical(res, data.frame())
})

test_that("`run_pathfindR()` argument checks work", {
  expect_error(
    run_pathfindR(example_pathfindR_input, plot_enrichment_chart = "INVALID"),
    "`plot_enrichment_chart` should be either TRUE or FALSE"
  )

  expect_error(
    run_pathfindR(example_pathfindR_input, list_active_snw_genes = "INVALID"),
    "`list_active_snw_genes` should be either TRUE or FALSE"
  )
})
