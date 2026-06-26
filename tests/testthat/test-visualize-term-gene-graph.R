test_that("`create_term_gene_graph()` -- check arguments", {
  ## Checking Error handling
  expect_error(
    create_term_gene_graph(list(ID = 1, lowest_p = 0.01, Up_regulated = "A", Down_regulated = "B")),
    "`result_df` should be a data.frame!"
  )

  bad_df <- data.frame(
    lowest_p = 0.01, Up_regulated = "A", Down_regulated = "B"
  )
  expect_error(
    create_term_gene_graph(bad_df),
    "All of ID, lowest_p, Up_regulated, Down_regulated must be present in `results_df`!"
  )

  expect_error(
    create_term_gene_graph(bad_df, use_description = TRUE),
    "All of Term_Description, lowest_p, Up_regulated, Down_regulated must be present in `results_df`!"
  )

  expect_error(
    create_term_gene_graph(example_pathfindR_output, genes_df = list(Gene.symbol = "A", logFC = 1)),
    "`genes_df` should be a data.frame!"
  )

  expect_error(
    create_term_gene_graph(example_pathfindR_output, term_fill = "nonexistent"),
    "`term_fill` is not found in the supplied `result_df`!"
  )

  expect_error(
    create_term_gene_graph(example_pathfindR_output, term_size = "nonexistent"),
    "`term_size` should be one of"
  )

  expect_error(
    create_term_gene_graph(example_pathfindR_output, num_terms = "five"),
    "`num_terms` must either be numeric or NULL!"
  )

  expect_error(
    create_term_gene_graph(example_pathfindR_output, use_description = "FALSE"),
    "`use_description` must either be TRUE or FALSE!"
  )

  expect_error(
    create_term_gene_graph(example_pathfindR_output, use_edge_weights = "FALSE"),
    "`use_edge_weights` must either be TRUE or FALSE!"
  )
})

test_that("`create_term_gene_graph()` -- check igraph creation", {
  genes_df <- example_pathfindR_input[1:3, ]
  colnames(genes_df) <- c("X", "Y", "Z")
  expect_error(
    create_term_gene_graph(example_pathfindR_output, genes_df)
  )

  ## Checking functional behavior
  genes_df <- example_pathfindR_input[1:5, ]
  input_terms_df <- example_pathfindR_output[1:10, ]

  expect_is(g <- create_term_gene_graph(input_terms_df), "igraph")
  expect_null(igraph::E(g)$weight)
  expect_null(igraph::V(g)$logFC)
  expect_null(igraph::V(g)$term_fill)
  expect_is(igraph::V(g)$size, "numeric")
  expect_equal(sum(igraph::V(g)$type == "term"), 10)

  expect_is(g <- create_term_gene_graph(input_terms_df, genes_df), "igraph")
  expect_null(igraph::E(g)$weight)
  expect_is(igraph::V(g)$logFC, "numeric")
  expect_null(igraph::V(g)$term_fill)
  expect_is(igraph::V(g)$size, "numeric")
  expect_equal(sum(igraph::V(g)$type == "term"), 10)

  expect_is(g <- create_term_gene_graph(input_terms_df, genes_df, term_fill = "Fold_Enrichment"), "igraph")
  expect_null(igraph::E(g)$weight)
  expect_is(igraph::V(g)$logFC, "numeric")
  expect_is(igraph::V(g)$term_fill, "numeric")
  expect_is(igraph::V(g)$size, "numeric")
  expect_equal(sum(igraph::V(g)$type == "term"), 10)

  expect_is(g <- create_term_gene_graph(input_terms_df, genes_df, term_fill = "Fold_Enrichment", use_edge_weights = TRUE), "igraph")
  expect_is(igraph::E(g)$weight, "numeric")
  expect_is(igraph::V(g)$logFC, "numeric")
  expect_is(igraph::V(g)$term_fill, "numeric")
  expect_is(igraph::V(g)$size, "numeric")
  expect_equal(sum(igraph::V(g)$type == "term"), 10)

  expect_is(g <- create_term_gene_graph(input_terms_df, genes_df, term_fill = "Fold_Enrichment", use_edge_weights = TRUE, term_size = "p_val"), "igraph")
  expect_is(igraph::E(g)$weight, "numeric")
  expect_is(igraph::V(g)$logFC, "numeric")
  expect_is(igraph::V(g)$term_fill, "numeric")
  expect_is(igraph::V(g)$size, "numeric")
  expect_equal(sum(igraph::V(g)$type == "term"), 10)

  expect_is(g <- create_term_gene_graph(input_terms_df, genes_df, term_fill = "Fold_Enrichment", use_edge_weights = TRUE, num_terms = 3), "igraph")
  expect_is(igraph::E(g)$weight, "numeric")
  expect_is(igraph::V(g)$logFC, "numeric")
  expect_is(igraph::V(g)$term_fill, "numeric")
  expect_is(igraph::V(g)$size, "numeric")
  expect_equal(sum(igraph::V(g)$type == "term"), 3)

  ## Corrects `num_terms` to maximum number of rows of input
  expect_is(g <- create_term_gene_graph(input_terms_df, num_terms = 150), "igraph")
  expect_equal(sum(igraph::V(g)$type == "term"), nrow(input_terms_df))
})


test_that("`create_term_gene_plot()` -- check arguments", {
  ## Checking Error handling
  expect_error(
    create_term_gene_plot(list()),
    "`graph` needs to be of class 'igraph'!"
  )

  genes_df <- example_pathfindR_input[1:3, ]
  g <- create_term_gene_graph(example_pathfindR_output, genes_df, term_fill = "Fold_Enrichment")

  expect_error(
    create_term_gene_plot(g, gene_node_fill = c("green", "red")),
    "`gene_node_fill` needs to be of length 3!"
  )

  expect_error(
    create_term_gene_plot(g, gene_node_fill = c("green", "red", "invalid")),
    "Not all elements in `gene_node_fill` are valid colors!"
  )

  expect_error(
    create_term_gene_plot(g, term_node_fill = c("#CCBB44", "#4477AA")),
    "`term_node_fill` needs to be of length 3!"
  )

  expect_error(
    create_term_gene_plot(g, term_node_fill = c("#CCBB44", "invalid", "#4477AA")),
    "Not all elements in `term_node_fill` are valid colors!"
  )

  expect_error(
    create_term_gene_plot(g, gene_node_color = c("green")),
    "`gene_node_color` needs to be of length 2!"
  )

  expect_error(
    create_term_gene_plot(g, gene_node_color = c("green", "red", "blue")),
    "`gene_node_color` needs to be of length 2!"
  )

  expect_error(
    create_term_gene_plot(g, gene_node_color = c("green", "notacolor")),
    "Not all elements in `gene_node_color` are valid colors!"
  )

  expect_error(
    create_term_gene_plot(g, term_node_color = "#INVALID"),
    "`term_node_color` is not a valid color!"
  )
})

test_that("`create_term_gene_plot()` -- Check ggraph creation", {
  input_terms_df <- example_pathfindR_output[1:10, ]
  ## Default functionality
  g0 <- create_term_gene_graph(input_terms_df, term_fill = "Fold_Enrichment")
  expect_is(create_term_gene_plot(g0), "ggraph")

  g0 <- create_term_gene_graph(input_terms_df)
  expect_is(create_term_gene_plot(g0), "ggraph")

  ## `genes_df` is included
  genes_df <- example_pathfindR_input[1:10, ]
  g1 <- create_term_gene_graph(input_terms_df, genes_df)
  expect_is(create_term_gene_plot(g1), "ggraph")

  g2 <- create_term_gene_graph(input_terms_df, genes_df, term_fill = "Fold_Enrichment")
  expect_is(create_term_gene_plot(g2), "ggraph")

  g3 <- create_term_gene_graph(input_terms_df, genes_df, term_fill = "Fold_Enrichment", use_edge_weights = TRUE)
  expect_is(create_term_gene_plot(g3), "ggraph")

  expect_is(create_term_gene_plot(g3, term_fill_label = "Fold Enrichment"), "ggraph")
  expect_is(create_term_gene_plot(g3, term_size_label = "# genes"), "ggraph")

  expect_is(create_term_gene_plot(g0, layout = "stress"), "ggraph")
  expect_is(create_term_gene_plot(g0, layout = "kk"), "ggraph")
  expect_is(create_term_gene_plot(g0, layout = "fr"), "ggraph")
  expect_error(create_term_gene_plot(g0, layout = "INVALID"))
})
