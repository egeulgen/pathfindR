# set up input data
input_data_frame <- example_pathfindR_input[1:10, c(1, 3)]
colnames(input_data_frame) <- c("GENE", "P_VALUE")

example_snws_len <- 1000
example_snw_output <- system.file("extdata", "resultActiveSubnetworkSearch.txt",
  package = "pathfindR"
)
mock_file_path <- function(...) {
  args <- list(...)
  if (args[[1]] == "active_snw_search") {
    return(example_snw_output)
  }
  return(file.path(...))
}

test_that("`get_active_subnetworks()` -- returns a list object", {
  mockery::stub(get_active_subnetworks, "dir.exists", TRUE)
  mockery::stub(get_active_subnetworks, "file.exists", TRUE)
  mockery::stub(get_active_subnetworks, "normalizePath", NULL)
  mockery::stub(get_active_subnetworks, "system", NULL)
  mockery::stub(get_active_subnetworks, "file.path", mock_file_path)
  mockery::stub(get_active_subnetworks, "file.rename", NULL)

  # Expect > 0 active snws
  expect_message(
    snw_list <- get_active_subnetworks(input_for_search = input_data_frame),
    "Found [1-9]\\d* active subnetworks"
  )
  expect_is(snw_list, "list")
  expect_is(snw_list[[1]], "character")
  expect_true(length(snw_list) > 0)

  # Expect no active snws
  mockery::stub(get_active_subnetworks, "filter_active_subnetworks", NULL)
  expect_message(
    snw_list <- get_active_subnetworks(input_for_search = input_data_frame),
    "Found 0 active subnetworks"
  )
  expect_identical(snw_list, list())
})

test_that("`get_active_subnetworks()` -- `dir_for_parallel_run` arg is used when provided", {
  mockery::stub(get_active_subnetworks, "dir.exists", TRUE)
  mockery::stub(get_active_subnetworks, "file.exists", TRUE)
  mockery::stub(get_active_subnetworks, "normalizePath", NULL)
  mockery::stub(get_active_subnetworks, "system", NULL)
  mockery::stub(get_active_subnetworks, "file.path", mock_file_path)
  mockery::stub(get_active_subnetworks, "file.rename", NULL)

  m <- mockery::mock(NULL, cycle = TRUE)
  mockery::stub(get_active_subnetworks, "setwd", m)
  res <- get_active_subnetworks(input_for_search = input_data_frame, dir_for_parallel_run = tempdir())
  mockery::expect_called(m, 2)
})

test_that("`get_active_subnetworks()` -- argument checks work", {
  # input_for_search
  expect_error(snw_list <- get_active_subnetworks(input_for_search = list()), "`input_for_search` should be data frame")

  invalid_input <- input_data_frame
  colnames(invalid_input) <- c("A", "B")
  expect_error(
    snw_list <- get_active_subnetworks(input_for_search = invalid_input),
    paste0("`input_for_search` should contain the columns ", paste(dQuote(c(
      "GENE",
      "P_VALUE"
    )), collapse = ","))
  )

  # snws_file
  expect_error(snw_list <- get_active_subnetworks(
    input_for_search = input_data_frame,
    snws_file = "[/]"
  ), "`snws_file` may be containing forbidden characters. Please change and try again")

  # search_method
  valid_mets <- c("GR", "SA", "GA")
  expect_error(
    get_active_subnetworks(input_for_search = input_data_frame, search_method = "INVALID"),
    paste0("`search_method` should be one of ", paste(dQuote(valid_mets), collapse = ", "))
  )

  # silent_option
  expect_error(
    get_active_subnetworks(input_for_search = input_data_frame, silent_option = "WRONG"),
    "`silent_option` should be either TRUE or FALSE"
  )

  expect_error(
    get_active_subnetworks(input_for_search = input_data_frame, use_all_positives = "INVALID"),
    "`use_all_positives` should be either TRUE or FALSE"
  )
})

test_that("`get_active_subnetworks()` -- all search methods work", {
  skip_on_cran()
  ## GR
  expect_message(
    snw_list <- get_active_subnetworks(
      input_for_search = input_data_frame,
      pin_name_path = "Biogrid", search_method = "GR", dir_for_parallel_run = tempdir(check = TRUE)
    ),
    "Found [1-9]\\d* active subnetworks"
  )
  expect_is(snw_list, "list")
  expect_is(snw_list[[1]], "character")

  skip("will test SA and GA if we can create a suitable (faster and non-empty) test case")
  ## SA
  expect_message(
    snw_list <- get_active_subnetworks(
      input_for_search = input_data_frame,
      pin_name_path = "Biogrid", search_method = "SA", dir_for_parallel_run = tempdir(check = TRUE)
    ),
    "Found [1-9]\\d* active subnetworks"
  )
  expect_is(snw_list, "list")
  expect_is(snw_list[[1]], "character")

  ## GA
  expect_message(
    snw_list <- get_active_subnetworks(
      input_for_search = input_data_frame,
      pin_name_path = "Biogrid", search_method = "GA", dir_for_parallel_run = tempdir(check = TRUE)
    ),
    "Found [1-9]\\d* active subnetworks"
  )
  expect_is(snw_list, "list")
  expect_is(snw_list[[1]], "character")
})

test_that("`get_active_subnetworks()` -- results are reproducible", {
  skip_on_cran()
  snw_lists <- list()
  seed_vals <- c(123, 123, 456)
  for (idx in 1:3) {
    seed <- seed_vals[idx]
    snw_lists[[idx]] <- get_active_subnetworks(
      input_for_search = input_data_frame,
      seedForRandom = seed, dir_for_parallel_run = tempdir(check = TRUE)
    )
  }
  expect_identical(snw_lists[[1]], snw_lists[[2]])
  expect_false(identical(snw_lists[[1]], snw_lists[[3]]))
})
