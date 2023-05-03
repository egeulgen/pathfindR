##################################################
## Package: pathfindR
## Script purpose: Unit testing script for
## various visualization functions
## Date: Apr 27, 2023
## Author: Ege Ulgen
##################################################

single_result <- example_pathfindR_output[1, ]
processed_input <- suppressMessages(
  input_processing(example_pathfindR_input, 0.05, "Biogrid")
)

# visualize_terms ---------------------------------------------------------
test_that("`visualize_terms()` creates expected png", {
  skip_on_cran()

  ## non-hsa-KEGG (visualize_term_interactions)
  fname <- gsub("\\s", "_", single_result$Term_Description)
  expected_out_file <- file.path(
    "term_visualizations", paste0(fname, ".png")
  )
  suppressMessages(
    visualize_terms(
      result_df = single_result,
      input_processed = processed_input,
      hsa_KEGG = FALSE,
      pin_name_path = "Biogrid"
    )
  )
  expect_true(file.exists(expected_out_file))
  unlink("term_visualizations", recursive = TRUE)

  suppressMessages(
    visualize_terms(
      result_df = single_result,
      input_processed = processed_input,
      hsa_KEGG = FALSE,
      pin_name_path = "GeneMania"
    )
  )
  expect_true(file.exists(expected_out_file))
  unlink("term_visualizations", recursive = TRUE)

  ## hsa KEGG
  expected_out_file <- file.path(
    "term_visualizations",
    paste0(single_result$ID, "_pathfindR.png")
  )
  suppressMessages(
    visualize_terms(
      result_df = single_result,
      input_processed = processed_input,
      hsa_KEGG = TRUE
    )
  )
  expect_true(file.exists(expected_out_file))
  unlink("term_visualizations", recursive = TRUE)
})


test_that("`visualize_terms()` argumment checks work", {
  expect_error(
    visualize_terms(result_df = "INVALID"),
    "`result_df` should be a data frame"
  )

  # hsa_KEGG = TRUE
  nec_cols <- "ID"
  expect_error(
    visualize_terms(single_result[, -1],
      hsa_KEGG = TRUE
    ),
    paste0(
      "`result_df` should contain the following columns: ",
      paste(dQuote(nec_cols), collapse = ", ")
    )
  )
  # hsa_KEGG = FALSE
  nec_cols <- c("Term_Description", "Up_regulated", "Down_regulated")
  expect_error(
    visualize_terms(single_result[, -2],
      hsa_KEGG = FALSE
    ),
    paste0(
      "`result_df` should contain the following columns: ",
      paste(dQuote(nec_cols), collapse = ", ")
    )
  )

  expect_error(
    visualize_terms(
      result_df = single_result,
      hsa_KEGG = TRUE
    ),
    "`input_processed` should be specified when `hsa_KEGG = TRUE`"
  )

  expect_error(
    visualize_terms(
      result_df = single_result,
      hsa_KEGG = "INVALID"
    ),
    "the argument `hsa_KEGG` should be either TRUE or FALSE"
  )
})

# visualize_term_interactions ---------------------------------------------
test_that("`visualize_term_interactions()` creates expected png file(s)", {
  skip_on_cran()
  expected_out_file <- file.path(
    "term_visualizations",
    paste0(gsub("\\s", "_", single_result$Term_Description), ".png")
  )
  expect_null(
    visualize_term_interactions(single_result, pin_name_path = "Biogrid")
  )
  expect_true(file.exists(expected_out_file))
  unlink("term_visualizations", recursive = TRUE)

  tmp_res <- rbind(single_result, single_result)
  tmp_res$Term_Description[2] <- "SKIP"
  tmp_res$Up_regulated[2] <- "Gene1"
  tmp_res$Down_regulated[2] <- ""
  expect_message(
    visualize_term_interactions(tmp_res, pin_name_path = "KEGG"),
    paste0(
      "< 2 genes, skipping visualization of ",
      tmp_res$Term_Description[2]
    )
  )
  unlink("term_visualizations", recursive = TRUE)

  # Non-empty non_Signif_Snw_Genes
  tmp_res <- single_result
  tmp_res$non_Signif_Snw_Genes <- example_pathfindR_output$Up_regulated[2]
  expect_null(
    visualize_term_interactions(tmp_res, pin_name_path = "Biogrid")
  )
  expect_true(file.exists(expected_out_file))
  unlink("term_visualizations", recursive = TRUE)
})

# visualize_hsa_KEGG ------------------------------------------------------
test_that("`visualize_hsa_KEGG()` creates expected png file(s)", {
  skip_on_cran()

  expected_out_file <- file.path(
    "term_visualizations", paste0(single_result$ID, "_pathfindR.png")
  )

  ###### Continuous change values
  ### scale_vals = TRUE
  expect_null(
    visualize_hsa_KEGG(
      hsa_kegg_ids = single_result$ID,
      input_processed = processed_input
    )
  )
  expect_true(file.exists(expected_out_file))
  unlink("term_visualizations", recursive = TRUE)

  expect_null(
    visualize_hsa_KEGG(
      hsa_kegg_ids = single_result$ID,
      input_processed = processed_input,
      node_cols = c("red", "green", "blue")
    )
  )
  expect_true(file.exists(expected_out_file))
  unlink("term_visualizations", recursive = TRUE)

  ### scale_vals = FALSE
  expect_null(
    visualize_hsa_KEGG(
      hsa_kegg_ids = single_result$ID,
      input_processed = processed_input,
      scale_vals = FALSE
    )
  )
  expect_true(file.exists(expected_out_file))
  unlink("term_visualizations", recursive = TRUE)

  ###### Constant change values (if no change values supplied in input)
  processed_input_no_change <- processed_input
  processed_input_no_change$CHANGE <- 1e6
  expect_null(
    visualize_hsa_KEGG(
      hsa_kegg_ids = single_result$ID,
      input_processed = processed_input_no_change
    )
  )
  expect_true(file.exists(expected_out_file))
  unlink("term_visualizations", recursive = TRUE)

  expect_message(
    visualize_hsa_KEGG(
      hsa_kegg_ids = single_result$ID,
      input_processed = processed_input_no_change,
      node_cols = "red"
    ),
    "all `change_vec` values are 1e6, using the first color in `node_cols`"
  )
  expect_true(file.exists(expected_out_file))
  unlink("term_visualizations", recursive = TRUE)

  ###### skips NULL
  temp_res <- example_pathfindR_output[1:2, ]
  temp_res$ID[2] <- "hsa12345"
  expect_null(
    visualize_hsa_KEGG(
      hsa_kegg_ids = temp_res$ID,
      input_processed = processed_input
    )
  )
  expect_true(
    file.exists(
      file.path("term_visualizations", paste0(temp_res$ID[1], "_pathfindR.png"))
    )
  )
  unlink("term_visualizations", recursive = TRUE)
})

temp_ids <- example_pathfindR_output$ID[1:2]
test_that("KEGML download error is handled properly", {
  skip_on_cran()
  temp_ids[2] <- "hsa12345"
  expected_out_files <- file.path(
    "term_visualizations",
    paste0(temp_ids, "_pathfindR.png")
  )
  visualize_hsa_KEGG(
    hsa_kegg_ids = temp_ids,
    input_processed = processed_input
  )
  expect_equal(file.exists(expected_out_files), c(TRUE, FALSE))
  unlink("term_visualizations", recursive = TRUE)
})

test_that("`visualize_hsa_KEGG()` argument checks work", {
  expect_error(
    visualize_hsa_KEGG(
      hsa_kegg_ids = list(),
      input_processed = processed_input
    ),
    "`hsa_kegg_ids` should be a vector of hsa KEGG IDs"
  )
  expect_error(
    visualize_hsa_KEGG(
      hsa_kegg_ids = c("X", "Y", "Z"),
      input_processed = processed_input
    ),
    "`hsa_kegg_ids` should be a vector of valid hsa KEGG IDs"
  )

  expect_error(
    visualize_hsa_KEGG(
      hsa_kegg_ids = single_result$ID,
      input_processed = list()
    ),
    "`input_processed` should be a data frame"
  )
  expect_error(
    visualize_hsa_KEGG(
      hsa_kegg_ids = single_result$ID,
      input_processed = processed_input[, -2]
    ),
    paste0(
      "`input_processed` should contain the following columns: ",
      paste(dQuote(c("GENE", "CHANGE")), collapse = ", ")
    )
  )

  expect_error(
    visualize_hsa_KEGG(
      hsa_kegg_ids = single_result$ID,
      input_processed = processed_input,
      scale_vals = "INVALID"
    ),
    "`scale_vals` should be logical"
  )

  expect_error(
    visualize_hsa_KEGG(
      hsa_kegg_ids = single_result$ID,
      input_processed = processed_input,
      node_cols = list()
    ),
    "`node_cols` should be a vector of colors"
  )
  expect_error(
    visualize_hsa_KEGG(
      hsa_kegg_ids = single_result$ID,
      input_processed = processed_input,
      node_cols = rep("red", 4)
    ),
    "the length of `node_cols` should be 3"
  )
  expect_error(
    visualize_hsa_KEGG(
      hsa_kegg_ids = single_result$ID,
      input_processed = processed_input,
      node_cols = c("red", "#FFFFFF", "INVALID")
    ),
    "`node_cols` should be a vector of valid colors"
  )
})

# color_kegg_pathway ------------------------------------------------------
test_that("`color_kegg_pathway()` exceptions are handled properly", {
  expect_null(
    suppressWarnings(
      color_kegg_pathway(pw_id = "hsa03040", change_vec = NULL)
    )
  )
  expect_message(
    color_kegg_pathway(pw_id = "hsa11111", change_vec = c())
  )

  expect_message(
    tmp <- obtain_KEGGML_URL("INVAID", tempfile())
  )
  mockery::stub(obtain_KEGGML_URL, "KEGGgraph::getKGMLurl", function(...) stop())
  expect_message(
    tmp <- obtain_KEGGML_URL("hsa11111", tempfile())
  )

  expect_message(
    tmp2 <- obtain_colored_url("INVALID", "INVALID", "white", "white")
  )
  expect_message(
    tmp3 <- suppressWarnings(download_kegg_png("INVALID", tempfile()))
  )

  expect_true(is.na(tmp))
  expect_true(is.na(tmp2))
  expect_true(is.na(tmp3))
})

# enrichment_chart --------------------------------------------------------
test_that("enrichment_chart produces a ggplot object with correct labels", {
  skip_on_cran()
  # default - top 10
  expect_is(g <- enrichment_chart(example_pathfindR_output), "ggplot")
  expect_equal(ggplot2::quo_name(g$mapping$x), "Fold_Enrichment")
  expect_equal(ggplot2::quo_name(g$mapping$y), "Term_Description")
  expect_equal(g$labels$size, "# genes")
  expect_equal(g$labels$colour, expression(-log[10](p)))
  expect_equal(g$labels$x, "Fold Enrichment")
  expect_equal(g$labels$y, "Term_Description")

  # plot_by_cluster
  expect_is(
    g <- enrichment_chart(example_pathfindR_output_clustered, plot_by_cluster = TRUE),
    "ggplot"
  )
  expect_equal(ggplot2::quo_name(g$mapping$x), "Fold_Enrichment")
  expect_equal(ggplot2::quo_name(g$mapping$y), "Term_Description")
  expect_equal(g$labels$size, "# genes")
  expect_equal(g$labels$colour, expression(-log[10](p)))
  expect_equal(g$labels$x, "Fold Enrichment")
  expect_equal(g$labels$y, "Term_Description")

  # chang top_terms
  expect_is(
    g <- enrichment_chart(example_pathfindR_output, top_terms = NULL),
    "ggplot"
  )
  expect_equal(ggplot2::quo_name(g$mapping$x), "Fold_Enrichment")
  expect_equal(ggplot2::quo_name(g$mapping$y), "Term_Description")
  expect_equal(g$labels$size, "# genes")
  expect_equal(g$labels$colour, expression(-log[10](p)))
  expect_equal(g$labels$x, "Fold Enrichment")
  expect_equal(g$labels$y, "Term_Description")

  expect_is(
    g <- enrichment_chart(example_pathfindR_output, top_terms = 1e3),
    "ggplot"
  )
  expect_equal(ggplot2::quo_name(g$mapping$x), "Fold_Enrichment")
  expect_equal(ggplot2::quo_name(g$mapping$y), "Term_Description")
  expect_equal(g$labels$size, "# genes")
  expect_equal(g$labels$colour, expression(-log[10](p)))
  expect_equal(g$labels$x, "Fold Enrichment")
  expect_equal(g$labels$y, "Term_Description")

  # change num_bubbles
  expect_is(
    g <- enrichment_chart(example_pathfindR_output_clustered,
      num_bubbles = 30
    ),
    "ggplot"
  )
  expect_equal(ggplot2::quo_name(g$mapping$x), "Fold_Enrichment")
  expect_equal(ggplot2::quo_name(g$mapping$y), "Term_Description")
  expect_equal(g$labels$size, "# genes")
  expect_equal(g$labels$colour, expression(-log[10](p)))
  expect_equal(g$labels$x, "Fold Enrichment")
  expect_equal(g$labels$y, "Term_Description")

  # change even_breaks
  expect_is(
    g <- enrichment_chart(example_pathfindR_output_clustered, even_breaks = FALSE),
    "ggplot"
  )
  expect_equal(ggplot2::quo_name(g$mapping$x), "Fold_Enrichment")
  expect_equal(ggplot2::quo_name(g$mapping$y), "Term_Description")
  expect_equal(g$labels$size, "# genes")
  expect_equal(g$labels$colour, expression(-log[10](p)))
  expect_equal(g$labels$x, "Fold Enrichment")
  expect_equal(g$labels$y, "Term_Description")
})

test_that("`enrichment_chart()` argument checks work", {
  necessary <- c(
    "Term_Description", "Fold_Enrichment", "lowest_p",
    "Up_regulated", "Down_regulated"
  )
  expect_error(
    enrichment_chart(example_pathfindR_output[, -2]),
    paste0(
      "The input data frame must have the columns:\n",
      paste(necessary, collapse = ", ")
    )
  )

  expect_error(
    enrichment_chart(example_pathfindR_output, plot_by_cluster = "INVALID"),
    "`plot_by_cluster` must be either TRUE or FALSE"
  )

  expect_message(
    enrichment_chart(example_pathfindR_output, plot_by_cluster = TRUE),
    "For plotting by cluster, there must a column named `Cluster` in the input data frame!"
  )

  expect_error(
    enrichment_chart(example_pathfindR_output, top_terms = "INVALID"),
    "`top_terms` must be either numeric or NULL"
  )

  expect_error(
    enrichment_chart(example_pathfindR_output, top_terms = 0),
    "`top_terms` must be > 1"
  )
})

# term_gene_graph ---------------------------------------------------------
test_that("`term_gene_graph()` produces a ggplot object using the correct data", {
  skip_on_cran()
  # Top 10 (default)
  expect_is(p <- term_gene_graph(example_pathfindR_output), "ggplot")
  expect_equal(sum(p$data$type == "term"), 10)

  # Top 3
  expect_is(p <- term_gene_graph(example_pathfindR_output, num_terms = 3), "ggplot")
  expect_equal(sum(p$data$type == "term"), 3)

  # All terms
  expect_is(p <- term_gene_graph(example_pathfindR_output[1:15, ], num_terms = NULL), "ggplot")
  expect_equal(sum(p$data$type == "term"), 15)

  # Top 1000, expect to plot top nrow(output)
  expect_is(p <- term_gene_graph(example_pathfindR_output[1:15, ], num_terms = 1e3), "ggplot")
  expect_equal(sum(p$data$type == "term"), 15)

  # use_description = TRUE
  expect_is(p <- term_gene_graph(example_pathfindR_output, use_description = TRUE), "ggplot")
  expect_equal(sum(p$data$type == "term"), 10)

  # node_size = "p_val"
  expect_is(p <- term_gene_graph(example_pathfindR_output, node_size = "p_val"), "ggplot")
  expect_equal(sum(p$data$type == "term"), 10)
})

test_that("`term_gene_graph()` argument checks work", {
  expect_error(
    term_gene_graph(example_pathfindR_output, num_terms = "INVALID"),
    "`num_terms` must either be numeric or NULL!"
  )

  expect_error(
    term_gene_graph(example_pathfindR_output, use_description = "INVALID"),
    "`use_description` must either be TRUE or FALSE!"
  )

  val_node_size <- c("num_genes", "p_val")
  expect_error(
    term_gene_graph(example_pathfindR_output, node_size = "INVALID"),
    paste0(
      "`node_size` should be one of ",
      paste(dQuote(val_node_size), collapse = ", ")
    )
  )

  expect_error(
    term_gene_graph(result_df = "INVALID"),
    "`result_df` should be a data frame"
  )

  wrong_df <- example_pathfindR_output[, -c(1, 2)]
  ID_column <- "ID"
  necessary_cols <- c(ID_column, "lowest_p", "Up_regulated", "Down_regulated")
  expect_error(
    term_gene_graph(wrong_df, use_description = FALSE),
    paste(c(
      "All of", paste(necessary_cols, collapse = ", "),
      "must be present in `results_df`!"
    ), collapse = " ")
  )

  ID_column <- "Term_Description"
  necessary_cols <- c(ID_column, "lowest_p", "Up_regulated", "Down_regulated")
  expect_error(
    term_gene_graph(wrong_df, use_description = TRUE),
    paste(c(
      "All of", paste(necessary_cols, collapse = ", "),
      "must be present in `results_df`!"
    ), collapse = " ")
  )
})

# term_gene_heatmap -------------------------------------------------------
test_that("`term_gene_heatmap()` produces a ggplot object using the correct data", {
  skip_on_cran()
  # Top 10 (default)
  expect_is(
    p <- term_gene_heatmap(example_pathfindR_output),
    "ggplot"
  )
  expect_equal(length(unique(p$data$Enriched_Term)), 10)
  expect_true(all(p$data$Enriched_Term %in% example_pathfindR_output$ID))

  # Top 3
  expect_is(
    p <- term_gene_heatmap(example_pathfindR_output, num_terms = 3),
    "ggplot"
  )
  expect_equal(length(unique(p$data$Enriched_Term)), 3)

  # No genes in "Down_regulated"
  res_df <- example_pathfindR_output[1:3, ]
  res_df$Down_regulated <- ""
  expect_is(p <- term_gene_heatmap(res_df), "ggplot")
  expect_equal(length(unique(p$data$Enriched_Term)), 3)

  # No genes in "Up_regulated"
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
    p <- term_gene_heatmap(example_pathfindR_output[1:15, ], num_terms = 1e3),
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
  expect_is(
    p <- term_gene_heatmap(example_pathfindR_output[1:3, ], example_pathfindR_input[, -2]),
    "ggplot"
  )

  # sort by lowest_p instead
  expect_is(
    p <- term_gene_heatmap(example_pathfindR_output[1:3, ], example_pathfindR_input, sort_terms_by_p = TRUE),
    "ggplot"
  )
})

test_that("`term_gene_graph()` argument checks work", {
  expect_error(
    term_gene_heatmap(result_df = example_pathfindR_output, use_description = "INVALID"),
    "`use_description` must either be TRUE or FALSE!"
  )

  expect_error(
    term_gene_heatmap(result_df = "INVALID"),
    "`result_df` should be a data frame"
  )

  wrong_df <- example_pathfindR_output[, -c(1, 2)]
  ID_column <- "ID"
  nec_cols <- c(ID_column, "lowest_p", "Up_regulated", "Down_regulated")
  expect_error(
    term_gene_heatmap(wrong_df, use_description = FALSE),
    paste0(
      "`result_df` should have the following columns: ",
      paste(dQuote(nec_cols), collapse = ", ")
    )
  )

  ID_column <- "Term_Description"
  nec_cols <- c(ID_column, "lowest_p", "Up_regulated", "Down_regulated")
  expect_error(
    term_gene_heatmap(wrong_df, use_description = TRUE),
    paste0(
      "`result_df` should have the following columns: ",
      paste(dQuote(nec_cols), collapse = ", ")
    )
  )

  expect_error(
    term_gene_heatmap(result_df = example_pathfindR_output, genes_df = "INVALID")
  )

  expect_error(
    term_gene_heatmap(result_df = example_pathfindR_output, num_terms = "INVALID"),
    "`num_terms` should be numeric or NULL"
  )

  expect_error(
    term_gene_heatmap(result_df = example_pathfindR_output, num_terms = -1),
    "`num_terms` should be > 0 or NULL"
  )
})

# UpSet_plot --------------------------------------------------------------
test_that("`UpSet_plot()` produces a ggplot object", {
  skip_on_cran()
  # Top 10 (default)
  expect_is(p <- UpSet_plot(example_pathfindR_output), "ggplot")

  # Top 3
  expect_is(p <- UpSet_plot(example_pathfindR_output, num_terms = 3), "ggplot")

  # All terms
  expect_is(p <- UpSet_plot(example_pathfindR_output[1:15, ], num_terms = NULL), "ggplot")

  # No genes in "Down_regulated"
  res_df <- example_pathfindR_output
  res_df$Down_regulated <- ""
  expect_is(p <- UpSet_plot(res_df, num_terms = 3), "ggplot")

  # No genes in "Up_regulated"
  res_df <- example_pathfindR_output
  res_df$Up_regulated <- ""
  expect_is(p <- UpSet_plot(res_df, num_terms = 3), "ggplot")

  # use_description = TRUE
  expect_is(
    p <- UpSet_plot(example_pathfindR_output, use_description = TRUE),
    "ggplot"
  )

  # Other visualization types
  expect_is(
    p <- UpSet_plot(example_pathfindR_output[1:3, ], example_pathfindR_input[1:10, ]),
    "ggplot"
  )
  expect_is(
    p <- UpSet_plot(example_pathfindR_output[1:3, ], example_pathfindR_input[1:10, ], method = "boxplot"),
    "ggplot"
  )
  expect_is(
    p <- UpSet_plot(example_pathfindR_output[1:3, ], method = "barplot"),
    "ggplot"
  )
})

test_that("`UpSet_plot()` argument checks work", {
  expect_error(
    UpSet_plot(result_df = example_pathfindR_output, use_description = "INVALID"),
    "`use_description` must either be TRUE or FALSE!"
  )

  expect_error(
    UpSet_plot(result_df = "INVALID"),
    "`result_df` should be a data frame"
  )

  wrong_df <- example_pathfindR_output[, -c(1, 2)]
  ID_column <- "ID"
  nec_cols <- c(ID_column, "lowest_p", "Up_regulated", "Down_regulated")
  expect_error(
    UpSet_plot(wrong_df, use_description = FALSE),
    paste0(
      "`result_df` should have the following columns: ",
      paste(dQuote(nec_cols), collapse = ", ")
    )
  )

  ID_column <- "Term_Description"
  nec_cols <- c(ID_column, "lowest_p", "Up_regulated", "Down_regulated")
  expect_error(
    UpSet_plot(wrong_df, use_description = TRUE),
    paste0(
      "`result_df` should have the following columns: ",
      paste(dQuote(nec_cols), collapse = ", ")
    )
  )

  expect_error(UpSet_plot(
    result_df = example_pathfindR_output,
    genes_df = "INVALID"
  ))

  expect_error(
    UpSet_plot(
      result_df = example_pathfindR_output,
      num_terms = "INVALID"
    ),
    "`num_terms` should be numeric or NULL"
  )

  expect_error(
    UpSet_plot(
      result_df = example_pathfindR_output,
      num_terms = -1
    ),
    "`num_terms` should be > 0 or NULL"
  )

  valid_opts <- c("heatmap", "boxplot", "barplot")
  expect_error(
    UpSet_plot(
      result_df = example_pathfindR_output,
      method = "INVALID"
    ),
    paste("`method` should be one of`", paste(dQuote(valid_opts), collapse = ", "))
  )

  expect_error(
    UpSet_plot(
      result_df = example_pathfindR_output,
      method = "boxplot"
    ),
    "For `method = boxplot`, you must provide `genes_df`"
  )
})
