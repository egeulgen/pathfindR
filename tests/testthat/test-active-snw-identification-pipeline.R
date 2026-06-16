set.seed(123)

# set up input data
n_input_genes <- 10
input_genes <- paste0("GENE", seq_len(n_input_genes))
input_p_vals <- runif(n_input_genes, min = 1e-10, max = 0.001)
input_data_frame <- data.frame(GENE = input_genes, P_VALUE = input_p_vals)

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

network <- build_network(sif_file)
mock_params <- list(
  p_for_nonsignificant = 0.5,
  seed                 = 1234L
)
score_context <- build_score_context(network, input_data_frame, mock_params)

test_that("`get_active_subnetworks()` -- returns a list object", {
  # Expect > 0 active snws
  expect_message(
    snw_list <- get_active_subnetworks(significant_genes = input_genes, network, score_context),
    "Found [1-9]\\d* active subnetworks"
  )
  expect_is(snw_list, "list")
  expect_is(snw_list[[1]], "character")
  expect_true(length(snw_list) > 0)

  # Expect no active snws
  mockery::stub(get_active_subnetworks, "filter_active_subnetworks", NULL)
  expect_message(
    snw_list <- get_active_subnetworks(significant_genes = input_genes, network, score_context),
    "Found 0 active subnetworks"
  )
  expect_identical(snw_list, list())
})

test_that("`get_active_subnetworks()` -- argument checks work", {
  # search_method
  valid_mets <- c("GR", "SA", "GA")
  expect_error(
    get_active_subnetworks(significant_genes = input_genes, network, score_context, search_method = "INVALID"),
    paste0("`search_method` should be one of ", paste(dQuote(valid_mets), collapse = ", "))
  )

  # verbose
  expect_error(
    get_active_subnetworks(significant_genes = input_genes, network, score_context, verbose = "WRONG"),
    "`verbose` should be either TRUE or FALSE"
  )

  expect_error(
    get_active_subnetworks(significant_genes = input_genes, network, score_context, start_with_all_positives = "INVALID"),
    "`start_with_all_positives` should be either TRUE or FALSE"
  )
})

test_that("`get_active_subnetworks()` -- all search methods work", {
  ## GR
  expect_message(
    snw_list <- get_active_subnetworks(
      significant_genes = input_genes, network, score_context,
      search_method = "GR"
    ),
    "Found [1-9]\\d* active subnetworks"
  )
  expect_is(snw_list, "list")
  expect_is(snw_list[[1]], "character")

  ## SA
  expect_message(
    snw_list <- get_active_subnetworks(
      significant_genes = input_genes, network, score_context,
      search_method = "SA", sig_gene_thr = 0, score_quan_thr = -1 # needed not to filter on toy results
    ),
    "Found [1-9]\\d* active subnetworks"
  )
  expect_is(snw_list, "list")
  expect_is(snw_list[[1]], "character")

  ## GA
  expect_message(
    snw_list <- get_active_subnetworks(
      significant_genes = input_genes, network, score_context,
      search_method = "GA"
    ),
    "Found [1-9]\\d* active subnetworks"
  )
  expect_is(snw_list, "list")
  expect_is(snw_list[[1]], "character")
})
