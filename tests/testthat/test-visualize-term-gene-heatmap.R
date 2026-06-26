test_that("`term_gene_heatmap()` -- produces a ggplot object using the correct data", {
  # Top 10 (default)
  expect_is(p <- term_gene_heatmap(example_pathfindR_output), "ggplot")
  expect_equal(length(unique(p$data$Enriched_Term)), 10)
  expect_true(all(p$data$Enriched_Term %in% example_pathfindR_output$ID))

  # Top 3
  expect_is(
    p <- term_gene_heatmap(example_pathfindR_output, num_terms = 3),
    "ggplot"
  )
  expect_equal(length(unique(p$data$Enriched_Term)), 3)

  # No genes in 'Down_regulated'
  res_df <- example_pathfindR_output[1:3, ]
  res_df$Down_regulated <- ""
  expect_is(p <- term_gene_heatmap(res_df), "ggplot")
  expect_equal(length(unique(p$data$Enriched_Term)), 3)

  # No genes in 'Up_regulated'
  res_df <- example_pathfindR_output[1:3, ]
  res_df$Up_regulated <- ""
  expect_is(p <- term_gene_heatmap(res_df), "ggplot")
  expect_equal(length(unique(p$data$Enriched_Term)), 3)

  # All terms
  expect_is(
    p <- term_gene_heatmap(example_pathfindR_output[1:15, ], num_terms = NULL),
    "ggplot"
  )
  expect_equal(length(unique(p$data$Enriched_Term)), 15)

  # Top 1000, expect to plot top nrow(output)
  expect_is(
    p <- term_gene_heatmap(example_pathfindR_output[1:15, ], num_terms = 1000),
    "ggplot"
  )
  expect_equal(length(unique(p$data$Enriched_Term)), 15)

  # use_description = TRUE
  expect_is(
    p <- term_gene_heatmap(example_pathfindR_output, use_description = TRUE),
    "ggplot"
  )
  expect_equal(length(unique(p$data$Enriched_Term)), 10)
  expect_true(all(p$data$Enriched_Term %in% example_pathfindR_output$Term_Description))

  # genes_df supplied
  expect_is(
    p <- term_gene_heatmap(example_pathfindR_output[1:3, ], example_pathfindR_input),
    "ggplot"
  )

  # genes_df supplied - wihout change column
  expect_is(p <- term_gene_heatmap(example_pathfindR_output[1:3, ], example_pathfindR_input[
    ,
    -2
  ]), "ggplot")

  # sort by lowest_p instead
  expect_is(p <- term_gene_heatmap(example_pathfindR_output[1:3, ], example_pathfindR_input,
    sort_terms_by_p = TRUE
  ), "ggplot")
})

test_that("`term_gene_graph()` -- argument checks work", {
  expect_error(
    term_gene_heatmap(result_df = example_pathfindR_output, use_description = "INVALID"),
    "`use_description` must either be TRUE or FALSE!"
  )

  expect_error(term_gene_heatmap(result_df = "INVALID"), "`result_df` should be a data frame")

  wrong_df <- example_pathfindR_output[, -c(1, 2)]
  ID_column <- "ID"
  nec_cols <- c(ID_column, "lowest_p", "Up_regulated", "Down_regulated")
  expect_error(term_gene_heatmap(wrong_df, use_description = FALSE), paste0(
    "`result_df` should have the following columns: ",
    paste(dQuote(nec_cols), collapse = ", ")
  ))

  ID_column <- "Term_Description"
  nec_cols <- c(ID_column, "lowest_p", "Up_regulated", "Down_regulated")
  expect_error(term_gene_heatmap(wrong_df, use_description = TRUE), paste0(
    "`result_df` should have the following columns: ",
    paste(dQuote(nec_cols), collapse = ", ")
  ))

  expect_error(term_gene_heatmap(result_df = example_pathfindR_output, genes_df = "INVALID"))

  expect_error(
    term_gene_heatmap(result_df = example_pathfindR_output, num_terms = "INVALID"),
    "`num_terms` should be numeric or NULL"
  )

  expect_error(
    term_gene_heatmap(result_df = example_pathfindR_output, num_terms = -1),
    "`num_terms` should be > 0 or NULL"
  )

  expect_error(term_gene_heatmap(example_pathfindR_output, low = ""))
  expect_error(term_gene_heatmap(example_pathfindR_output, mid = ""))
  expect_error(term_gene_heatmap(example_pathfindR_output, high = ""))
})
