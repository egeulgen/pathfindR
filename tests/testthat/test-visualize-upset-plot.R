test_that("`UpSet_plot()` -- produces a ggplot object", {
  # Top 10 (default)
  expect_is(p <- UpSet_plot(example_pathfindR_output), "ggplot")

  # Top 3
  expect_is(p <- UpSet_plot(example_pathfindR_output, num_terms = 3), "ggplot")

  # All terms
  expect_is(
    p <- UpSet_plot(example_pathfindR_output[1:15, ], num_terms = NULL),
    "ggplot"
  )

  # No genes in 'Down_regulated'
  res_df <- example_pathfindR_output
  res_df$Down_regulated <- ""
  expect_is(p <- UpSet_plot(res_df, num_terms = 3), "ggplot")

  # No genes in 'Up_regulated'
  res_df <- example_pathfindR_output
  res_df$Up_regulated <- ""
  expect_is(p <- UpSet_plot(res_df, num_terms = 3), "ggplot")

  # use_description = TRUE
  expect_is(
    p <- UpSet_plot(example_pathfindR_output, use_description = TRUE),
    "ggplot"
  )

  # Other visualization types
  expect_is(p <- UpSet_plot(example_pathfindR_output[1:3, ], example_pathfindR_input[1:10, ]), "ggplot")
  expect_is(p <- UpSet_plot(example_pathfindR_output[1:3, ], example_pathfindR_input[1:10, ], method = "boxplot"), "ggplot")
  expect_is(
    p <- UpSet_plot(example_pathfindR_output[1:3, ], method = "barplot"),
    "ggplot"
  )
})

test_that("`UpSet_plot()` -- argument checks work", {
  expect_error(
    UpSet_plot(result_df = example_pathfindR_output, use_description = "INVALID"),
    "`use_description` must either be TRUE or FALSE!"
  )

  expect_error(UpSet_plot(result_df = "INVALID"), "`result_df` should be a data frame")

  wrong_df <- example_pathfindR_output[, -c(1, 2)]
  ID_column <- "ID"
  nec_cols <- c(ID_column, "lowest_p", "Up_regulated", "Down_regulated")
  expect_error(UpSet_plot(wrong_df, use_description = FALSE), paste0(
    "`result_df` should have the following columns: ",
    paste(dQuote(nec_cols), collapse = ", ")
  ))

  ID_column <- "Term_Description"
  nec_cols <- c(ID_column, "lowest_p", "Up_regulated", "Down_regulated")
  expect_error(UpSet_plot(wrong_df, use_description = TRUE), paste0(
    "`result_df` should have the following columns: ",
    paste(dQuote(nec_cols), collapse = ", ")
  ))

  expect_error(UpSet_plot(result_df = example_pathfindR_output, genes_df = "INVALID"))

  expect_error(
    UpSet_plot(result_df = example_pathfindR_output, num_terms = "INVALID"),
    "`num_terms` should be numeric or NULL"
  )

  expect_error(
    UpSet_plot(result_df = example_pathfindR_output, num_terms = -1),
    "`num_terms` should be > 0 or NULL"
  )

  valid_opts <- c("heatmap", "boxplot", "barplot")
  expect_error(
    UpSet_plot(result_df = example_pathfindR_output, method = "INVALID"),
    paste("`method` should be one of`", paste(dQuote(valid_opts), collapse = ", "))
  )

  expect_error(
    UpSet_plot(result_df = example_pathfindR_output, method = "boxplot"),
    "For `method = boxplot`, you must provide `genes_df`"
  )

  expect_error(UpSet_plot(example_pathfindR_output, low = ""))
  expect_error(UpSet_plot(example_pathfindR_output, mid = ""))
  expect_error(UpSet_plot(example_pathfindR_output, high = ""))
})
