test_that("`enrichment_chart()` -- produces a ggplot object with correct labels", {
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
  expect_is(
    g <- enrichment_chart(example_pathfindR_output_clustered, plot_by_cluster = TRUE),
    "ggplot"
  )
  expect_equal(ggplot2::quo_name(g$mapping$x), "Fold_Enrichment")
  expect_equal(ggplot2::quo_name(g$mapping$y), "Term_Description")

  labels <- ggplot2::get_labs(g)
  expect_equal(labels$size, "# genes")
  expect_equal(labels$colour, expression(-log[10](p)))
  expect_equal(labels$x, "Fold Enrichment")
  expect_equal(labels$y, "Term_Description")

  # chang top_terms
  expect_is(
    g <- enrichment_chart(example_pathfindR_output, top_terms = NULL),
    "ggplot"
  )
  expect_equal(ggplot2::quo_name(g$mapping$x), "Fold_Enrichment")
  expect_equal(ggplot2::quo_name(g$mapping$y), "Term_Description")

  labels <- ggplot2::get_labs(g)
  expect_equal(labels$size, "# genes")
  expect_equal(labels$colour, expression(-log[10](p)))
  expect_equal(labels$x, "Fold Enrichment")
  expect_equal(labels$y, "Term_Description")

  expect_is(
    g <- enrichment_chart(example_pathfindR_output, top_terms = 1000),
    "ggplot"
  )
  expect_equal(ggplot2::quo_name(g$mapping$x), "Fold_Enrichment")
  expect_equal(ggplot2::quo_name(g$mapping$y), "Term_Description")

  labels <- ggplot2::get_labs(g)
  expect_equal(labels$size, "# genes")
  expect_equal(labels$colour, expression(-log[10](p)))
  expect_equal(labels$x, "Fold Enrichment")
  expect_equal(labels$y, "Term_Description")

  # change num_bubbles
  expect_is(
    g <- enrichment_chart(example_pathfindR_output_clustered, num_bubbles = 30),
    "ggplot"
  )
  expect_equal(ggplot2::quo_name(g$mapping$x), "Fold_Enrichment")
  expect_equal(ggplot2::quo_name(g$mapping$y), "Term_Description")

  labels <- ggplot2::get_labs(g)
  expect_equal(labels$size, "# genes")
  expect_equal(labels$colour, expression(-log[10](p)))
  expect_equal(labels$x, "Fold Enrichment")
  expect_equal(labels$y, "Term_Description")

  # change even_breaks
  expect_is(
    g <- enrichment_chart(example_pathfindR_output_clustered, even_breaks = FALSE),
    "ggplot"
  )
  expect_equal(ggplot2::quo_name(g$mapping$x), "Fold_Enrichment")
  expect_equal(ggplot2::quo_name(g$mapping$y), "Term_Description")

  labels <- ggplot2::get_labs(g)
  expect_equal(labels$size, "# genes")
  expect_equal(labels$colour, expression(-log[10](p)))
  expect_equal(labels$x, "Fold Enrichment")
  expect_equal(labels$y, "Term_Description")

  # change even_breaks
  expect_is(
    g <- enrichment_chart(example_pathfindR_output_clustered, order_by = "Fold_Enrichment"),
    "ggplot"
  )
  expect_equal(ggplot2::quo_name(g$mapping$x), "Fold_Enrichment")
  expect_equal(ggplot2::quo_name(g$mapping$y), "Term_Description")

  labels <- ggplot2::get_labs(g)
  expect_equal(labels$size, "# genes")
  expect_equal(labels$colour, "Fold_Enrichment")
  expect_equal(labels$x, "Fold Enrichment")
  expect_equal(labels$y, "Term_Description")
})

test_that("`enrichment_chart()` -- order_by arg-related tests", {
  # Change order_by
  expect_is(
    g <- enrichment_chart(example_pathfindR_output_clustered, order_by = "highest_p"),
    "ggplot"
  )
  expect_equal(ggplot2::quo_name(g$mapping$x), "Fold_Enrichment")
  expect_equal(ggplot2::quo_name(g$mapping$y), "Term_Description")

  labels <- ggplot2::get_labs(g)
  expect_equal(labels$size, "# genes")
  expect_equal(labels$colour, expression(-log[10](p)))
  expect_equal(labels$x, "Fold Enrichment")
  expect_equal(labels$y, "Term_Description")

  # check if order is correct
  input_ordered <- example_pathfindR_output_clustered[order(example_pathfindR_output_clustered[["highest_p"]], decreasing = FALSE), ]
  expect_equal(input_ordered$ID[1:10], g$data$ID)
})

test_that("`enrichment_chart()` -- argument checks work", {
  necessary <- c(
    "Term_Description", "Fold_Enrichment", "lowest_p", "Up_regulated",
    "Down_regulated"
  )
  expect_error(enrichment_chart(example_pathfindR_output[, -2]), paste0(
    "The input data frame must have the columns:\n",
    paste(necessary, collapse = ", ")
  ))

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

  expect_error(enrichment_chart(example_pathfindR_output, top_terms = 0), "`top_terms` must be > 1")
})
