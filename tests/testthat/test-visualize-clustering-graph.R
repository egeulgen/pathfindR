enrichment_res <- example_pathfindR_output[1:5, ]
input_kappa_mat <- create_kappa_matrix(enrichment_res)

test_that("`cluster_graph_vis()` -- graph visualization of clusters works OK", {
  mockery::stub(hierarchical_term_clustering, "graphics::plot", NULL)
  mockery::stub(hierarchical_term_clustering, "stats::rect.hclust", NULL)
  mock_plot.igraph <- mockery::mock(NULL, cycle = TRUE)
  mockery::stub(cluster_graph_vis, "igraph::plot.igraph", mock_plot.igraph)
  ## use_description = FALSE
  for (clustering_func in c(hierarchical_term_clustering, fuzzy_term_clustering)) {
    clu_obj <- clustering_func(input_kappa_mat, enrichment_res)
    expect_silent(cluster_graph_vis(clu_obj, input_kappa_mat, enrichment_res))
  }
  mockery::expect_called(mock_plot.igraph, 2)
})

test_that("`cluster_graph_vis()` -- coloring of 'extra' clusters work", {
  mockery::stub(cluster_graph_vis, "igraph::plot.igraph", NULL)
  ### more than 41 clusters (number of colors available)
  selected_num_terms <- 45
  clu_input_df <- example_pathfindR_output[seq_len(selected_num_terms), ]
  mock_kappa_mat <- matrix(0,
    nrow = selected_num_terms, ncol = selected_num_terms,
    dimnames = list(clu_input_df$ID, clu_input_df$ID)
  )

  # dummy hierarchical result
  hierarchical_clu_obj <- seq_len(selected_num_terms)
  names(hierarchical_clu_obj) <- clu_input_df$ID
  expect_silent(cluster_graph_vis(hierarchical_clu_obj, mock_kappa_mat, clu_input_df))

  # dummy fuzzy result
  fuzzy_clu_obj <- matrix(FALSE,
    nrow = selected_num_terms, ncol = selected_num_terms,
    dimnames = list(clu_input_df$ID, seq_len(selected_num_terms))
  )
  diag(fuzzy_clu_obj) <- TRUE
  expect_silent(cluster_graph_vis(fuzzy_clu_obj, mock_kappa_mat, clu_input_df))
})

test_that("`cluster_graph_vis()` -- check errors are raised appropriately", {
  expect_error(cluster_graph_vis(list(), matrix(), data.frame(ID = 1)), "Invalid class for `clu_obj`!")

  # hierarchical - missing terms in kappa matrix
  clu_obj <- hierarchical_term_clustering(input_kappa_mat, enrichment_res, plot_dend = FALSE)
  expect_error(
    cluster_graph_vis(c(clu_obj, EXTRA = 1L), input_kappa_mat, enrichment_res),
    "Not all terms in `clu_obj` present in `kappa_mat`!"
  )

  # fuzzy - missing terms in kappa matrix
  clu_obj <- fuzzy_term_clustering(input_kappa_mat, enrichment_res)
  expect_error(cluster_graph_vis(
    rbind(clu_obj, EXTRA = rep(FALSE, ncol(clu_obj))),
    input_kappa_mat, enrichment_res
  ), "Not all terms in `clu_obj` present in `kappa_mat`!")
})
