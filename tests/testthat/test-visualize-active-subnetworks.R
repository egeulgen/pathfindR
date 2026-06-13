input_data_frame <- example_pathfindR_input[1:10, c(1, 3)]

example_snw_output <- system.file("extdata", "resultActiveSubnetworkSearch.txt",
  package = "pathfindR"
)
example_snws_len <- 1000

test_that("`visualize_active_subnetworks()` -- returns list of ggraph objects", {
  # empty file case
  empty_path <- tempfile("empty", fileext = ".txt")
  file.create(empty_path)
  expect_null(visualize_active_subnetworks(active_snw_path = empty_path, genes_df = input_data_frame))

  skip_on_cran()
  # default
  g_list <- visualize_active_subnetworks(example_snw_output, input_data_frame)
  expect_is(g_list, "list")
  expect_is(g_list[[1]], "ggraph")
  expect_true(length(g_list) <= example_snws_len)

  # set `num_snws` to larger than actual number
  g_list <- visualize_active_subnetworks(example_snw_output, input_data_frame,
    num_snws = 21
  )
  expect_is(g_list, "list")
  expect_is(g_list[[1]], "ggraph")
  expect_length(g_list, 13)
})
