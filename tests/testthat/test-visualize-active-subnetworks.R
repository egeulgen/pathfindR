input_data_frame <- example_pathfindR_input[1:10, c(1, 3)]



test_that("`visualize_active_subnetworks()` -- returns list of ggraph objects", {
  # empty snws case
  empty_snws <- list()
  expect_null(visualize_active_subnetworks(active_snws = empty_snws, genes_df = input_data_frame))

  skip_on_cran()

  example_snws_len <- length(example_unfiltered_snws)

  # default
  g_list <- visualize_active_subnetworks(example_unfiltered_snws, input_data_frame)
  expect_is(g_list, "list")
  expect_is(g_list[[1]], "ggraph")
  expect_true(length(g_list) <= example_snws_len)

  # set `num_snws` to larger than actual number
  g_list <- visualize_active_subnetworks(example_unfiltered_snws, input_data_frame,
    num_snws = 21
  )
  expect_is(g_list, "list")
  expect_is(g_list[[1]], "ggraph")
  expect_length(g_list, 13)
})
