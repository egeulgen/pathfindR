## Tests for functions related to various visualization functions - Apr 2024

single_result <- example_pathfindR_output[1, ]
processed_input <- example_pathfindR_input[, c(1, 1, 2, 3)]
colnames(processed_input) <- c("old_GENE", "GENE", "CHANGE", "P_VALUE")

test_that("`visualize_terms()` -- calls the appropriate function", {
    mock_vis_kegg <- mockery::mock(NULL)
    mockery::stub(visualize_terms, "visualize_KEGG_diagram", mock_vis_kegg)
    expect_silent(visualize_terms(result_df = single_result, input_processed = data.frame(),
        is_KEGG_result = TRUE))
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
    expect_error(visualize_terms(single_result[, -1], is_KEGG_result = TRUE), paste0("`result_df` should contain the following columns: ",
        paste(dQuote(nec_cols), collapse = ", ")))
    # is_KEGG_result = FALSE
    nec_cols <- c("Term_Description", "Up_regulated", "Down_regulated")
    expect_error(visualize_terms(single_result[, -2], is_KEGG_result = FALSE), paste0("`result_df` should contain the following columns: ",
        paste(dQuote(nec_cols), collapse = ", ")))

    expect_error(visualize_terms(result_df = single_result, is_KEGG_result = TRUE), "`input_processed` should be specified when `is_KEGG_result = TRUE`")

    expect_error(visualize_terms(result_df = single_result, is_KEGG_result = "INVALID"),
        "the argument `is_KEGG_result` should be either TRUE or FALSE")
})

test_that("`visualize_term_interactions()` -- creates expected list of ggraph objects", {
    skip_on_cran()
    expect_is(res <- visualize_term_interactions(single_result, pin_name_path = "Biogrid"), "list")
    expect_is(res[[1]], "ggraph")

    tmp_res <- rbind(single_result, single_result)
    tmp_res$Term_Description[2] <- "SKIP"
    tmp_res$Up_regulated[2] <- "Gene1"
    tmp_res$Down_regulated[2] <- ""
    expect_message(res <- visualize_term_interactions(tmp_res, pin_name_path = "KEGG"),
        paste0("< 2 genes, skipping visualization of ", tmp_res$Term_Description[2]))

    # Non-empty non_Signif_Snw_Genes
    tmp_res <- single_result
    tmp_res$non_Signif_Snw_Genes <- example_pathfindR_output$Up_regulated[2]
    expect_is(res <- visualize_term_interactions(tmp_res, pin_name_path = "Biogrid"), "list")
    expect_is(res[[1]], "ggraph")
})

test_that("`visualize_KEGG_diagram()` -- creates expected list of ggraph objects", {
    skip_on_cran()

    expect_is(res <- visualize_KEGG_diagram(kegg_pw_ids = single_result$ID, input_processed = processed_input), "list")
    expect_is(res[[1]], "ggraph")

    constant_input <- processed_input
    constant_input$CHANGE <- 1e+06
    expect_is(visualize_KEGG_diagram(kegg_pw_ids = single_result$ID, input_processed = constant_input), "list")
    expect_is(res[[1]], "ggraph")
})

test_that("`visualize_KEGG_diagram()` -- skips pathway if non-existent", {
    skip_on_cran()
    temp_res <- example_pathfindR_output[1:2, ]
    temp_res$ID[2] <- "hsa12345"

    expect_is(res <- visualize_KEGG_diagram(kegg_pw_ids = temp_res$ID, input_processed = processed_input), "list")
    expect_is(res[[1]], "ggraph")
    expect_length(expect_is, 1)
})

test_that("`visualize_KEGG_diagram()` -- argument checks work", {
    expect_error(visualize_KEGG_diagram(kegg_pw_ids = list(), input_processed = processed_input),
        "`kegg_pw_ids` should be a vector of KEGG IDs")
    expect_error(visualize_KEGG_diagram(kegg_pw_ids = c("X", "Y", "Z"), input_processed = processed_input),
        "`kegg_pw_ids` should be a vector of valid hsa KEGG IDs")

    expect_error(visualize_KEGG_diagram(kegg_pw_ids = "abc12345", input_processed = list()),
        "`input_processed` should be a data frame")
    expect_error(visualize_KEGG_diagram(kegg_pw_ids = "abc12345", input_processed = processed_input[,
        -2]), paste0("`input_processed` should contain the following columns: ",
        paste(dQuote(c("GENE", "CHANGE")), collapse = ", ")))
})

test_that("`color_kegg_pathway()` -- works as expected", {
  skip_on_cran()

  pw_id <- "hsa00010"
  change_vec <- c(-2, 4, 6)
  names(change_vec) <- c("hsa:2821", "hsa:226", "hsa:229")

  expect_is(result <- color_kegg_pathway(pw_id, change_vec), "ggraph")

  names(change_vec) <- rep("missing", 3)
  expect_is(result <- color_kegg_pathway(pw_id, change_vec), "NULL")
})

test_that("`color_kegg_pathway()` -- exceptions are handled properly", {
    change_vec <- c(-2, 4, 6)
    names(change_vec) <- c("hsa:2821", "hsa:226", "hsa:229")

    expect_error(color_kegg_pathway(pw_id = "hsa03040", change_vec = change_vec,
        scale_vals = "INVALID"), "`scale_vals` should be logical")
    expect_error(color_kegg_pathway(pw_id = "hsa03040", change_vec = change_vec,
        node_cols = list()), "`node_cols` should be a vector of colors")
    expect_error(color_kegg_pathway(pw_id = "hsa03040", change_vec = change_vec,
        node_cols = rep("red", 4)), "the length of `node_cols` should be 3")
    expect_error(color_kegg_pathway(pw_id = "hsa03040", change_vec = change_vec,
        node_cols = c("red", "#FFFFFF", "INVALID")), "`node_cols` should be a vector of valid colors")

    skip_on_cran()

    constant_vec <- rep(1e+06, 3)
    names(constant_vec) <- c("hsa:2821", "hsa:226", "hsa:229")

    expect_silent(color_kegg_pathway(pw_id = "hsa03040", change_vec = change_vec,
        node_cols = c("red", "blue", "green")))
    expect_message(color_kegg_pathway(pw_id = "hsa03040", change_vec = constant_vec,
        node_cols = c("red", "blue", "green")))

    expect_null(suppressWarnings(color_kegg_pathway(pw_id = "hsa03040", change_vec = NULL)))
    expect_message(color_kegg_pathway(pw_id = "hsa11111", change_vec = c()))
})

test_that("`enrichment_chart()` -- produces a ggplot object with correct labels",
    {
        # default - top 10
        expect_is(g <- enrichment_chart(example_pathfindR_output), "ggplot")
        expect_equal(ggplot2::quo_name(g$mapping$x), "Fold_Enrichment")
        expect_equal(ggplot2::quo_name(g$mapping$y), "Term_Description")
        
        labels <- ggplot2::get_labs(g)
        expect_equal(labels$size, "# genes")
        expect_equal(labels$colour, expression(-log[10](p)))
        expect_equal(labels$x, "Fold Enrichment")
        expect_equal(labels$y, "Term_Description")

        # plot_by_cluster
        expect_is(g <- enrichment_chart(example_pathfindR_output_clustered, plot_by_cluster = TRUE),
            "ggplot")
        expect_equal(ggplot2::quo_name(g$mapping$x), "Fold_Enrichment")
        expect_equal(ggplot2::quo_name(g$mapping$y), "Term_Description")
        
        labels <- ggplot2::get_labs(g)
        expect_equal(labels$size, "# genes")
        expect_equal(labels$colour, expression(-log[10](p)))
        expect_equal(labels$x, "Fold Enrichment")
        expect_equal(labels$y, "Term_Description")

        # chang top_terms
        expect_is(g <- enrichment_chart(example_pathfindR_output, top_terms = NULL),
            "ggplot")
        expect_equal(ggplot2::quo_name(g$mapping$x), "Fold_Enrichment")
        expect_equal(ggplot2::quo_name(g$mapping$y), "Term_Description")
        
        labels <- ggplot2::get_labs(g)
        expect_equal(labels$size, "# genes")
        expect_equal(labels$colour, expression(-log[10](p)))
        expect_equal(labels$x, "Fold Enrichment")
        expect_equal(labels$y, "Term_Description")

        expect_is(g <- enrichment_chart(example_pathfindR_output, top_terms = 1000),
            "ggplot")
        expect_equal(ggplot2::quo_name(g$mapping$x), "Fold_Enrichment")
        expect_equal(ggplot2::quo_name(g$mapping$y), "Term_Description")
        
        labels <- ggplot2::get_labs(g)
        expect_equal(labels$size, "# genes")
        expect_equal(labels$colour, expression(-log[10](p)))
        expect_equal(labels$x, "Fold Enrichment")
        expect_equal(labels$y, "Term_Description")

        # change num_bubbles
        expect_is(g <- enrichment_chart(example_pathfindR_output_clustered, num_bubbles = 30),
            "ggplot")
        expect_equal(ggplot2::quo_name(g$mapping$x), "Fold_Enrichment")
        expect_equal(ggplot2::quo_name(g$mapping$y), "Term_Description")
        
        labels <- ggplot2::get_labs(g)
        expect_equal(labels$size, "# genes")
        expect_equal(labels$colour, expression(-log[10](p)))
        expect_equal(labels$x, "Fold Enrichment")
        expect_equal(labels$y, "Term_Description")

        # change even_breaks
        expect_is(g <- enrichment_chart(example_pathfindR_output_clustered, even_breaks = FALSE),
            "ggplot")
        expect_equal(ggplot2::quo_name(g$mapping$x), "Fold_Enrichment")
        expect_equal(ggplot2::quo_name(g$mapping$y), "Term_Description")
        
        labels <- ggplot2::get_labs(g)
        expect_equal(labels$size, "# genes")
        expect_equal(labels$colour, expression(-log[10](p)))
        expect_equal(labels$x, "Fold Enrichment")
        expect_equal(labels$y, "Term_Description")
    })

test_that("`enrichment_chart()` -- argument checks work", {
    necessary <- c("Term_Description", "Fold_Enrichment", "lowest_p", "Up_regulated",
        "Down_regulated")
    expect_error(enrichment_chart(example_pathfindR_output[, -2]), paste0("The input data frame must have the columns:\n",
        paste(necessary, collapse = ", ")))

    expect_error(enrichment_chart(example_pathfindR_output, plot_by_cluster = "INVALID"),
        "`plot_by_cluster` must be either TRUE or FALSE")

    expect_message(enrichment_chart(example_pathfindR_output, plot_by_cluster = TRUE),
        "For plotting by cluster, there must a column named `Cluster` in the input data frame!")

    expect_error(enrichment_chart(example_pathfindR_output, top_terms = "INVALID"),
        "`top_terms` must be either numeric or NULL")

    expect_error(enrichment_chart(example_pathfindR_output, top_terms = 0), "`top_terms` must be > 1")
})

test_that("`term_gene_graph()` -- produces a ggplot object using the correct data",
    {
        # Top 10 (default)
        expect_is(p <- term_gene_graph(example_pathfindR_output), "ggplot")
        expect_equal(sum(p$data$type == "term"), 10)

        # Top 3
        expect_is(p <- term_gene_graph(example_pathfindR_output, num_terms = 3),
            "ggplot")
        expect_equal(sum(p$data$type == "term"), 3)

        # All terms
        expect_is(p <- term_gene_graph(example_pathfindR_output[1:15, ], num_terms = NULL),
            "ggplot")
        expect_equal(sum(p$data$type == "term"), 15)

        # Top 1000, expect to plot top nrow(output)
        expect_is(p <- term_gene_graph(example_pathfindR_output[1:15, ], num_terms = 1000),
            "ggplot")
        expect_equal(sum(p$data$type == "term"), 15)

        # use_description = TRUE
        expect_is(p <- term_gene_graph(example_pathfindR_output, use_description = TRUE),
            "ggplot")
        expect_equal(sum(p$data$type == "term"), 10)

        # node_size = 'p_val'
        expect_is(p <- term_gene_graph(example_pathfindR_output, node_size = "p_val"),
            "ggplot")
        expect_equal(sum(p$data$type == "term"), 10)
    })

test_that("`term_gene_graph()` -- argument checks work", {
    expect_error(term_gene_graph(example_pathfindR_output, num_terms = "INVALID"),
        "`num_terms` must either be numeric or NULL!")

    expect_error(term_gene_graph(example_pathfindR_output, use_description = "INVALID"),
        "`use_description` must either be TRUE or FALSE!")

    val_node_size <- c("num_genes", "p_val")
    expect_error(term_gene_graph(example_pathfindR_output, node_size = "INVALID"),
        paste0("`node_size` should be one of ", paste(dQuote(val_node_size), collapse = ", ")))

    expect_error(term_gene_graph(result_df = "INVALID"), "`result_df` should be a data frame")

    wrong_df <- example_pathfindR_output[, -c(1, 2)]
    ID_column <- "ID"
    necessary_cols <- c(ID_column, "lowest_p", "Up_regulated", "Down_regulated")
    expect_error(term_gene_graph(wrong_df, use_description = FALSE), paste(c("All of",
        paste(necessary_cols, collapse = ", "), "must be present in `results_df`!"),
        collapse = " "))

    ID_column <- "Term_Description"
    necessary_cols <- c(ID_column, "lowest_p", "Up_regulated", "Down_regulated")
    expect_error(term_gene_graph(wrong_df, use_description = TRUE), paste(c("All of",
        paste(necessary_cols, collapse = ", "), "must be present in `results_df`!"),
        collapse = " "))

    expect_error(term_gene_graph(example_pathfindR_output, node_colors = list()))
    expect_error(term_gene_graph(example_pathfindR_output, node_colors = c(1, 2, 3)))
    expect_error(term_gene_graph(example_pathfindR_output, node_colors = c("red", "blue")))
})

test_that("`term_gene_heatmap()` -- produces a ggplot object using the correct data",
    {
        skip_on_cran()
        # Top 10 (default)
        expect_is(p <- term_gene_heatmap(example_pathfindR_output), "ggplot")
        expect_equal(length(unique(p$data$Enriched_Term)), 10)
        expect_true(all(p$data$Enriched_Term %in% example_pathfindR_output$ID))

        # Top 3
        expect_is(p <- term_gene_heatmap(example_pathfindR_output, num_terms = 3),
            "ggplot")
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
        expect_is(p <- term_gene_heatmap(example_pathfindR_output[1:15, ], num_terms = NULL),
            "ggplot")
        expect_equal(length(unique(p$data$Enriched_Term)), 15)

        # Top 1000, expect to plot top nrow(output)
        expect_is(p <- term_gene_heatmap(example_pathfindR_output[1:15, ], num_terms = 1000),
            "ggplot")
        expect_equal(length(unique(p$data$Enriched_Term)), 15)

        # use_description = TRUE
        expect_is(p <- term_gene_heatmap(example_pathfindR_output, use_description = TRUE),
            "ggplot")
        expect_equal(length(unique(p$data$Enriched_Term)), 10)
        expect_true(all(p$data$Enriched_Term %in% example_pathfindR_output$Term_Description))

        # genes_df supplied
        expect_is(p <- term_gene_heatmap(example_pathfindR_output[1:3, ], example_pathfindR_input),
            "ggplot")

        # genes_df supplied - wihout change column
        expect_is(p <- term_gene_heatmap(example_pathfindR_output[1:3, ], example_pathfindR_input[,
            -2]), "ggplot")

        # sort by lowest_p instead
        expect_is(p <- term_gene_heatmap(example_pathfindR_output[1:3, ], example_pathfindR_input,
            sort_terms_by_p = TRUE), "ggplot")
    })

test_that("`term_gene_graph()` -- argument checks work", {
    expect_error(term_gene_heatmap(result_df = example_pathfindR_output, use_description = "INVALID"),
        "`use_description` must either be TRUE or FALSE!")

    expect_error(term_gene_heatmap(result_df = "INVALID"), "`result_df` should be a data frame")

    wrong_df <- example_pathfindR_output[, -c(1, 2)]
    ID_column <- "ID"
    nec_cols <- c(ID_column, "lowest_p", "Up_regulated", "Down_regulated")
    expect_error(term_gene_heatmap(wrong_df, use_description = FALSE), paste0("`result_df` should have the following columns: ",
        paste(dQuote(nec_cols), collapse = ", ")))

    ID_column <- "Term_Description"
    nec_cols <- c(ID_column, "lowest_p", "Up_regulated", "Down_regulated")
    expect_error(term_gene_heatmap(wrong_df, use_description = TRUE), paste0("`result_df` should have the following columns: ",
        paste(dQuote(nec_cols), collapse = ", ")))

    expect_error(term_gene_heatmap(result_df = example_pathfindR_output, genes_df = "INVALID"))

    expect_error(term_gene_heatmap(result_df = example_pathfindR_output, num_terms = "INVALID"),
        "`num_terms` should be numeric or NULL")

    expect_error(term_gene_heatmap(result_df = example_pathfindR_output, num_terms = -1),
        "`num_terms` should be > 0 or NULL")

    expect_error(term_gene_heatmap(example_pathfindR_output, low = ""))
    expect_error(term_gene_heatmap(example_pathfindR_output, mid = ""))
    expect_error(term_gene_heatmap(example_pathfindR_output, high = ""))
})

test_that("`UpSet_plot()` -- produces a ggplot object", {
    skip_on_cran()
    # Top 10 (default)
    expect_is(p <- UpSet_plot(example_pathfindR_output), "ggplot")

    # Top 3
    expect_is(p <- UpSet_plot(example_pathfindR_output, num_terms = 3), "ggplot")

    # All terms
    expect_is(p <- UpSet_plot(example_pathfindR_output[1:15, ], num_terms = NULL),
        "ggplot")

    # No genes in 'Down_regulated'
    res_df <- example_pathfindR_output
    res_df$Down_regulated <- ""
    expect_is(p <- UpSet_plot(res_df, num_terms = 3), "ggplot")

    # No genes in 'Up_regulated'
    res_df <- example_pathfindR_output
    res_df$Up_regulated <- ""
    expect_is(p <- UpSet_plot(res_df, num_terms = 3), "ggplot")

    # use_description = TRUE
    expect_is(p <- UpSet_plot(example_pathfindR_output, use_description = TRUE),
        "ggplot")

    # Other visualization types
    expect_is(p <- UpSet_plot(example_pathfindR_output[1:3, ], example_pathfindR_input[1:10,
        ]), "ggplot")
    expect_is(p <- UpSet_plot(example_pathfindR_output[1:3, ], example_pathfindR_input[1:10,
        ], method = "boxplot"), "ggplot")
    expect_is(p <- UpSet_plot(example_pathfindR_output[1:3, ], method = "barplot"),
        "ggplot")
})

test_that("`UpSet_plot()` -- argument checks work", {
    expect_error(UpSet_plot(result_df = example_pathfindR_output, use_description = "INVALID"),
        "`use_description` must either be TRUE or FALSE!")

    expect_error(UpSet_plot(result_df = "INVALID"), "`result_df` should be a data frame")

    wrong_df <- example_pathfindR_output[, -c(1, 2)]
    ID_column <- "ID"
    nec_cols <- c(ID_column, "lowest_p", "Up_regulated", "Down_regulated")
    expect_error(UpSet_plot(wrong_df, use_description = FALSE), paste0("`result_df` should have the following columns: ",
        paste(dQuote(nec_cols), collapse = ", ")))

    ID_column <- "Term_Description"
    nec_cols <- c(ID_column, "lowest_p", "Up_regulated", "Down_regulated")
    expect_error(UpSet_plot(wrong_df, use_description = TRUE), paste0("`result_df` should have the following columns: ",
        paste(dQuote(nec_cols), collapse = ", ")))

    expect_error(UpSet_plot(result_df = example_pathfindR_output, genes_df = "INVALID"))

    expect_error(UpSet_plot(result_df = example_pathfindR_output, num_terms = "INVALID"),
        "`num_terms` should be numeric or NULL")

    expect_error(UpSet_plot(result_df = example_pathfindR_output, num_terms = -1),
        "`num_terms` should be > 0 or NULL")

    valid_opts <- c("heatmap", "boxplot", "barplot")
    expect_error(UpSet_plot(result_df = example_pathfindR_output, method = "INVALID"),
        paste("`method` should be one of`", paste(dQuote(valid_opts), collapse = ", ")))

    expect_error(UpSet_plot(result_df = example_pathfindR_output, method = "boxplot"),
        "For `method = boxplot`, you must provide `genes_df`")

    expect_error(UpSet_plot(example_pathfindR_output, low = ""))
    expect_error(UpSet_plot(example_pathfindR_output, mid = ""))
    expect_error(UpSet_plot(example_pathfindR_output, high = ""))
})

test_that("`isColor()` -- identifies colors correctly", {
  expect_true(isColor("red"))
  expect_true(isColor("green"))
  expect_true(isColor("black"))
  expect_true(isColor("gray60"))
  expect_true(isColor("#E5D7BF"))

  expect_false(isColor(""))
  expect_false(isColor("a"))
  expect_false(isColor(FALSE))
  expect_false(isColor(1))
  expect_false(isColor(c()))
  expect_false(isColor(list()))
})
