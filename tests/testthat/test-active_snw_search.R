##################################################
## Package: pathfindR
## Script purpose: Unit testing script for
## functions related to active subnetwork searh
## Date: May 6, 2023
## Author: Ege Ulgen
##################################################

# active_snw_search -------------------------------------------------------
input_df1 <- suppressMessages(input_processing(example_pathfindR_input[1:100, ],
  p_val_threshold = 0.05,
  pin_name_path = "Biogrid"
))
input_df2 <- suppressMessages(input_processing(example_pathfindR_input[1:2, ],
  p_val_threshold = 0.05,
  pin_name_path = "Biogrid"
))

test_that("`active_snw_search()` returns list object", {
  skip_on_cran()
  # Expect > 0 active snws
  expect_message(
    snw_list <- active_snw_search(input_for_search = input_df1),
    "Found [1-9]\\d* active subnetworks"
  )
  expect_is(snw_list, "list")
  expect_is(snw_list[[1]], "character")
  expect_true(length(snw_list) > 0)
  unlink("active_snw_search", recursive = TRUE)

  # dir_for_parallel_run works
  dummy_dir <- file.path(tempdir(check = TRUE), "dummy_dir")
  dir.create(dummy_dir)
  expect_message(
    snw_list <- active_snw_search(
      input_for_search = input_df1,
      dir_for_parallel_run = dummy_dir
    ),
    "Found [1-9]\\d* active subnetworks"
  )
  expect_true(file.exists(file.path(dummy_dir, "active_snw_search/active_snws.txt")))


  # Expect no active snws
  expect_message(
    snw_list <- active_snw_search(
      input_for_search = input_df2,
      sig_gene_thr = 1
    ),
    "Found 0 active subnetworks"
  )
  expect_identical(snw_list, list())
  unlink("active_snw_search", recursive = TRUE)

  mockery::stub(active_snw_search, "filterActiveSnws", NULL)
  expect_message(
    snw_list <- active_snw_search(input_for_search = input_df2),
    "Found 0 active subnetworks"
  )
  expect_identical(snw_list, list())
  unlink("active_snw_search", recursive = TRUE)
})

test_that("All search methods for `active_snw_search()` work", {
  skip_on_cran()
  ## GR
  expect_message(
    snw_list <- active_snw_search(
      input_for_search = input_df1,
      pin_name_path = "Biogrid",
      search_method = "GR"
    ),
    "Found [1-9]\\d* active subnetworks"
  )
  expect_is(snw_list, "list")
  expect_is(snw_list[[1]], "character")
  unlink("active_snw_search", recursive = TRUE)

  skip("will test SA and GA if we can create a suitable (faster and non-empty) test case")
  ## SA
  expect_message(
    snw_list <- active_snw_search(
      input_for_search = input_df1,
      pin_name_path = "Biogrid",
      search_method = "SA"
    ),
    "Found [1-9]\\d* active subnetworks"
  )
  expect_is(snw_list, "list")
  expect_is(snw_list[[1]], "character")
  unlink("active_snw_search", recursive = TRUE)

  ## GA
  expect_message(
    snw_list <- active_snw_search(
      input_for_search = input_df1,
      pin_name_path = "Biogrid",
      search_method = "GA"
    ),
    "Found [1-9]\\d* active subnetworks"
  )
  expect_is(snw_list, "list")
  expect_is(snw_list[[1]], "character")
  unlink("active_snw_search", recursive = TRUE)
})


test_that("`active_snw_search()` results are reproducible", {
  skip_on_cran()
  snw_list1 <- active_snw_search(input_for_search = input_df1, seedForRandom = 123)
  snw_list2 <- active_snw_search(input_for_search = input_df1, seedForRandom = 123)
  snw_list3 <- active_snw_search(input_for_search = input_df1, seedForRandom = 456)
  expect_identical(snw_list1, snw_list2)
  expect_false(identical(snw_list1, snw_list3))
})




test_that("`active_snw_search()` argument checks work", {
  # input_for_search
  expect_error(
    snw_list <- active_snw_search(input_for_search = list()),
    "`input_for_search` should be data frame"
  )

  invalid_input <- input_df2[, 3:4]
  cnames <- c("GENE", "P_VALUE")
  expect_error(
    snw_list <- active_snw_search(input_for_search = invalid_input),
    paste0(
      "`input_for_search` should contain the columns ",
      paste(dQuote(cnames), collapse = ",")
    )
  )

  # snws_file
  expect_error(
    snw_list <- active_snw_search(
      input_for_search = input_df2,
      snws_file = "[/]"
    ),
    "`snws_file` may be containing forbidden characters. Please change and try again"
  )

  # search_method
  valid_mets <- c("GR", "SA", "GA")
  expect_error(
    active_snw_search(
      input_for_search = input_df2,
      search_method = "INVALID"
    ),
    paste0(
      "`search_method` should be one of ",
      paste(dQuote(valid_mets), collapse = ", ")
    )
  )

  # silent_option
  expect_error(
    active_snw_search(
      input_for_search = input_df2,
      silent_option = "WRONG"
    ),
    "`silent_option` should be either TRUE or FALSE"
  )

  expect_error(
    active_snw_search(
      input_for_search = input_df2,
      use_all_positives = "INVALID"
    ),
    "`use_all_positives` should be either TRUE or FALSE"
  )
})



# filterActiveSnws --------------------------------------------------------
sample_path <- system.file("extdata/resultActiveSubnetworkSearch.txt",
  package = "pathfindR"
)
example_snws_len <- 20

test_that("`filterActiveSnws()` returns list object", {
  skip_on_cran()
  tmp_filtered <- filterActiveSnws(
    active_snw_path = sample_path,
    sig_genes_vec = example_pathfindR_input$Gene.symbol
  )
  expect_is(tmp_filtered, "list")
  expect_length(tmp_filtered, 2)
  expect_is(tmp_filtered$subnetworks, "list")
  expect_is(tmp_filtered$scores, "numeric")

  expect_is(tmp_filtered$subnetworks[[1]], "character")
  expect_true(length(tmp_filtered$subnetworks) <= example_snws_len)

  # empty file case
  empty_path <- tempfile("empty", fileext = ".txt")
  file.create(empty_path)
  expect_null(suppressWarnings(filterActiveSnws(
    active_snw_path = empty_path,
    sig_genes_vec = example_pathfindR_input$Gene.symbol
  )))
})

test_that("`score_quan_thr` in `filterActiveSnws()` works", {
  skip_on_cran()
  tmp_filtered <- filterActiveSnws(
    active_snw_path = sample_path,
    sig_genes_vec = example_pathfindR_input$Gene.symbol,
    score_quan_thr = -1,
    sig_gene_thr = 0
  )
  expect_length(tmp_filtered$subnetworks, example_snws_len)

  for (q_thr in seq(.1, 1, by = .1)) {
    tmp_filtered <- filterActiveSnws(
      active_snw_path = sample_path,
      sig_genes_vec = example_pathfindR_input$Gene.symbol,
      score_quan_thr = q_thr,
      sig_gene_thr = 0
    )
    exp_len <- example_snws_len * (1 - q_thr)
    expect_length(tmp_filtered$subnetworks, as.integer(exp_len + .5))
  }
})

test_that("`sig_gene_thr` in `filterActiveSnws()` works", {
  skip_on_cran()
  tmp_filtered1 <- filterActiveSnws(
    active_snw_path = sample_path,
    sig_genes_vec = example_pathfindR_input$Gene.symbol,
    sig_gene_thr = 0.02, # default
    score_quan_thr = -1
  )

  tmp_filtered2 <- filterActiveSnws(
    active_snw_path = sample_path,
    sig_genes_vec = example_pathfindR_input$Gene.symbol,
    sig_gene_thr = 0.1,
    score_quan_thr = -1
  )

  expect_true(length(tmp_filtered2$subnetworks) < example_snws_len)
  expect_true(length(tmp_filtered1$subnetworks) > length(tmp_filtered2$subnetworks))
})

test_that("`filterActiveSnws()` argument checks work", {
  expect_error(
    filterActiveSnws(active_snw_path = "this/is/not/a/valid/path"),
    "The active subnetwork file does not exist! Check the `active_snw_path` argument"
  )

  expect_error(
    filterActiveSnws(
      active_snw_path = sample_path,
      sig_genes_vec = list()
    ),
    "`sig_genes_vec` should be a vector"
  )

  expect_error(
    filterActiveSnws(
      active_snw_path = sample_path,
      sig_genes_vec = example_pathfindR_input$Gene.symbol,
      score_quan_thr = "INVALID"
    ),
    "`score_quan_thr` should be numeric"
  )
  expect_error(
    filterActiveSnws(
      active_snw_path = sample_path,
      sig_genes_vec = example_pathfindR_input$Gene.symbol,
      score_quan_thr = -2
    ),
    "`score_quan_thr` should be in \\[0, 1\\] or -1 \\(if not filtering\\)"
  )
  expect_error(
    filterActiveSnws(
      active_snw_path = sample_path,
      sig_genes_vec = example_pathfindR_input$Gene.symbol,
      score_quan_thr = 2
    ),
    "`score_quan_thr` should be in \\[0, 1\\] or -1 \\(if not filtering\\)"
  )

  expect_error(
    filterActiveSnws(
      active_snw_path = sample_path,
      sig_genes_vec = example_pathfindR_input$Gene.symbol,
      sig_gene_thr = "INVALID"
    ),
    "`sig_gene_thr` should be numeric"
  )
  expect_error(
    filterActiveSnws(
      active_snw_path = sample_path,
      sig_genes_vec = example_pathfindR_input$Gene.symbol,
      sig_gene_thr = -1
    ),
    "`sig_gene_thr` should be in \\[0, 1\\]"
  )
})

# visualize_active_subnetworks --------------------------------------------
test_that("`visualize_active_subnetworks()` returns list of ggraph objects", {
  # empty file case
  empty_path <- tempfile("empty", fileext = ".txt")
  file.create(empty_path)
  expect_null(visualize_active_subnetworks(
    active_snw_path = empty_path,
    genes_df = example_pathfindR_input[1:5, ]
  ))

  skip_on_cran()
  # default
  input_df <- example_pathfindR_input[1:10, ]
  g_list <- visualize_active_subnetworks(sample_path, input_df)
  expect_is(g_list, "list")
  expect_is(g_list[[1]], "ggraph")
  expect_true(length(g_list) <= example_snws_len)

  # set `num_snws` to larger than actual number
  g_list <- visualize_active_subnetworks(sample_path, input_df, num_snws = 21)
  expect_is(g_list, "list")
  expect_is(g_list[[1]], "ggraph")
  expect_length(g_list, 3)
})
