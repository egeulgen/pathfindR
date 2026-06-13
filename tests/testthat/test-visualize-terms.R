single_result <- example_pathfindR_output[1, ]
processed_input <- example_pathfindR_input[, c(1, 1, 2, 3)]
colnames(processed_input) <- c("old_GENE", "GENE", "CHANGE", "P_VALUE")

test_that("`visualize_terms()` -- calls the appropriate function", {
  mock_vis_kegg <- mockery::mock(NULL)
  mockery::stub(visualize_terms, "visualize_KEGG_diagram", mock_vis_kegg)
  expect_silent(visualize_terms(
    result_df = single_result, input_processed = data.frame(),
    is_KEGG_result = TRUE
  ))
  mockery::expect_called(mock_vis_kegg, 1)

  mock_vis_term_inter <- mockery::mock(NULL)
  mockery::stub(visualize_terms, "visualize_term_interactions", mock_vis_term_inter)
  expect_silent(visualize_terms(result_df = single_result, is_KEGG_result = FALSE))
  mockery::expect_called(mock_vis_term_inter, 1)
})

test_that("`visualize_terms()` -- argumment checks work", {
  expect_error(visualize_terms(result_df = "INVALID"), "`result_df` should be a data frame")

  # is_KEGG_result = TRUE
  nec_cols <- "ID"
  expect_error(visualize_terms(single_result[, -1], is_KEGG_result = TRUE), paste0(
    "`result_df` should contain the following columns: ",
    paste(dQuote(nec_cols), collapse = ", ")
  ))
  # is_KEGG_result = FALSE
  nec_cols <- c("Term_Description", "Up_regulated", "Down_regulated")
  expect_error(visualize_terms(single_result[, -2], is_KEGG_result = FALSE), paste0(
    "`result_df` should contain the following columns: ",
    paste(dQuote(nec_cols), collapse = ", ")
  ))

  expect_error(visualize_terms(result_df = single_result, is_KEGG_result = TRUE), "`input_processed` should be specified when `is_KEGG_result = TRUE`")

  expect_error(
    visualize_terms(result_df = single_result, is_KEGG_result = "INVALID"),
    "the argument `is_KEGG_result` should be either TRUE or FALSE"
  )
})

test_that("`visualize_term_interactions()` -- creates expected list of ggraph objects", {
  skip_on_cran()
  expect_is(res <- visualize_term_interactions(single_result, pin_name_path = "Biogrid"), "list")
  expect_is(res[[1]], "ggraph")

  tmp_res <- rbind(single_result, single_result)
  tmp_res$Term_Description[2] <- "SKIP"
  tmp_res$Up_regulated[2] <- "Gene1"
  tmp_res$Down_regulated[2] <- ""
  expect_message(
    res <- visualize_term_interactions(tmp_res, pin_name_path = "KEGG"),
    paste0("< 2 genes, skipping visualization of ", tmp_res$Term_Description[2])
  )

  # Non-empty non_Signif_Snw_Genes
  tmp_res <- single_result
  tmp_res$non_Signif_Snw_Genes <- example_pathfindR_output$Up_regulated[2]
  expect_is(res <- visualize_term_interactions(tmp_res, pin_name_path = "Biogrid"), "list")
  expect_is(res[[1]], "ggraph")
})

test_that("`visualize_KEGG_diagram()` -- creates expected list of ggraph objects", {
  skip_on_cran()
  skip_if_not_installed("org.Hs.eg.db")

  expect_is(res <- visualize_KEGG_diagram(kegg_pw_ids = single_result$ID, input_processed = processed_input), "list")
  expect_is(res[[1]], "ggraph")

  constant_input <- processed_input
  constant_input$CHANGE <- 1e+06
  expect_is(visualize_KEGG_diagram(kegg_pw_ids = single_result$ID, input_processed = constant_input), "list")
  expect_is(res[[1]], "ggraph")
})

test_that("`visualize_KEGG_diagram()` -- skips pathway if non-existent", {
  skip_on_cran()
  skip_if_not_installed("org.Hs.eg.db")
  temp_res <- example_pathfindR_output[1:2, ]
  temp_res$ID[2] <- "hsa12345"

  expect_is(res <- visualize_KEGG_diagram(kegg_pw_ids = temp_res$ID, input_processed = processed_input), "list")
  expect_is(res[[1]], "ggraph")
  expect_length(expect_is, 1)
})

test_that("`visualize_KEGG_diagram()` -- argument checks work", {
  expect_error(
    visualize_KEGG_diagram(kegg_pw_ids = list(), input_processed = processed_input),
    "`kegg_pw_ids` should be a vector of KEGG IDs"
  )
  expect_error(
    visualize_KEGG_diagram(kegg_pw_ids = c("X", "Y", "Z"), input_processed = processed_input),
    "`kegg_pw_ids` should be a vector of valid hsa KEGG IDs"
  )

  expect_error(
    visualize_KEGG_diagram(kegg_pw_ids = "abc12345", input_processed = list()),
    "`input_processed` should be a data frame"
  )
  expect_error(visualize_KEGG_diagram(kegg_pw_ids = "abc12345", input_processed = processed_input[
    ,
    -2
  ]), paste0(
    "`input_processed` should contain the following columns: ",
    paste(dQuote(c("GENE", "CHANGE")), collapse = ", ")
  ))
})
