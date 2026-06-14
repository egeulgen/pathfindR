set.seed(123)

# set up input data
input_genes <- paste0("GENE", 1:10)
input_p_vals <- runif(10, min = 0.00001, max = 0.05)
input_data_frame <- data.frame(GENE = input_genes, P_VALUE = input_p_vals)

pool <- paste0("GENE", 1:50)
n_edges <- 100
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

test_that("`get_active_subnetworks()` -- returns a list object", {
  # Expect > 0 active snws
  expect_message(
    snw_list <- get_active_subnetworks(input_for_search = input_data_frame, pin_name_path = sif_file),
    "Found [1-9]\\d* active subnetworks"
  )
  expect_is(snw_list, "list")
  expect_is(snw_list[[1]], "character")
  expect_true(length(snw_list) > 0)

  # Expect no active snws
  mockery::stub(get_active_subnetworks, "filter_active_subnetworks", NULL)
  expect_message(
    snw_list <- get_active_subnetworks(input_for_search = input_data_frame, pin_name_path = sif_file),
    "Found 0 active subnetworks"
  )
  expect_identical(snw_list, list())
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

  # search_method
  valid_mets <- c("GR", "SA", "GA")
  expect_error(
    get_active_subnetworks(input_for_search = input_data_frame, search_method = "INVALID"),
    paste0("`search_method` should be one of ", paste(dQuote(valid_mets), collapse = ", "))
  )

  # verbose
  expect_error(
    get_active_subnetworks(input_for_search = input_data_frame, verbose = "WRONG"),
    "`verbose` should be either TRUE or FALSE"
  )

  expect_error(
    get_active_subnetworks(input_for_search = input_data_frame, start_with_all_positives = "INVALID"),
    "`start_with_all_positives` should be either TRUE or FALSE"
  )
})

test_that("`get_active_subnetworks()` -- all search methods work", {
  ## GR
  expect_message(
    snw_list <- get_active_subnetworks(
      input_for_search = input_data_frame,
      pin_name_path = sif_file, search_method = "GR"
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
      pin_name_path = sif_file, search_method = "SA"
    ),
    "Found [1-9]\\d* active subnetworks"
  )
  expect_is(snw_list, "list")
  expect_is(snw_list[[1]], "character")

  ## GA
  expect_message(
    snw_list <- get_active_subnetworks(
      input_for_search = input_data_frame,
      pin_name_path = sif_file, search_method = "GA"
    ),
    "Found [1-9]\\d* active subnetworks"
  )
  expect_is(snw_list, "list")
  expect_is(snw_list[[1]], "character")
})

test_that("`get_active_subnetworks()` -- results are reproducible", {
  larger_pool <- paste0("GENE", 1:1000)
  toy_pin_df <- data.frame(
    InteractorA = sample(pool, 500, replace = TRUE),
    pp = "pp",
    InteractorB = sample(pool, 500, replace = TRUE),
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

  snw_lists <- list()
  seed_vals <- c(123, 456, 123)
  for (idx in 1:3) {
    seed <- seed_vals[idx]
    snw_lists[[idx]] <- get_active_subnetworks(
      input_for_search = input_data_frame,
      pin_name_path = sif_file,
      seed_for_stochastic_methods = seed
    )
  }
  expect_false(identical(snw_lists[[1]], snw_lists[[2]]))
  expect_identical(snw_lists[[1]], snw_lists[[3]])
})
