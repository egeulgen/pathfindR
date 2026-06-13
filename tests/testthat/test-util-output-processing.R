example_gene_data <- example_pathfindR_input[1:10, ]
colnames(example_gene_data) <- c("GENE", "CHANGE", "P_VALUE")
tmp_res <- example_pathfindR_output[1:5, -c(7, 8)]

test_that("`annotate_term_genes()` -- adds input genes for each term", {
  expect_is(
    annotated_result <- annotate_term_genes(result_df = tmp_res, input_processed = example_gene_data),
    "data.frame"
  )
  expect_true("Up_regulated" %in% colnames(annotated_result) & "Down_regulated" %in%
    colnames(annotated_result))
  expect_true(nrow(annotated_result) == nrow(tmp_res))
})

test_that("annotate_term_genes() -- argument checks work", {
  expect_error(
    annotate_term_genes(result_df = list(), input_processed = example_gene_data),
    "`result_df` should be a data frame"
  )
  expect_error(
    annotate_term_genes(result_df = tmp_res[, -1], input_processed = example_gene_data),
    "`result_df` should contain an \"ID\" column"
  )

  expect_error(
    annotate_term_genes(result_df = tmp_res, input_processed = list()),
    "`input_processed` should be a data frame"
  )
  expect_error(annotate_term_genes(result_df = tmp_res, input_processed = example_gene_data[
    ,
    -1
  ]), "`input_processed` should contain the columns \"GENE\" and \"CHANGE\"")


  expect_error(annotate_term_genes(
    result_df = tmp_res, input_processed = example_gene_data,
    genes_by_term = "INVALID"
  ), "`genes_by_term` should be a list of term gene sets")
  expect_error(annotate_term_genes(
    result_df = tmp_res, input_processed = example_gene_data,
    genes_by_term = list(1)
  ), "`genes_by_term` should be a named list \\(names are gene set IDs\\)")
})

test_that("`create_HTML_report()` -- works a expected", {
  mock_render <- mockery::mock(NULL, cycle = TRUE)
  mockery::stub(create_HTML_report, "rmarkdown::render", mock_render)

  create_HTML_report(
    input = data.frame(), input_processed = data.frame(), final_res = data.frame(),
    dir_for_report = "/path/to/report/dir"
  )
  mockery::expect_called(mock_render, 3)
})
