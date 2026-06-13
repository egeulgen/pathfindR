test_that("`configure_output_dir()` -- works as expected", {
  expected_dir <- file.path(tempdir(), "test_pathfindR_results")
  mockery::stub(configure_output_dir, "file.path", expected_dir)
  expect_equal(configure_output_dir(), expected_dir)

  test_out_dir <- file.path(tempdir(), "TEST")
  for (i in 1:3) {
    actual_dir <- configure_output_dir(test_out_dir)
    dir_to_check <- test_out_dir
    if (i > 1) {
      dir_to_check <- paste0(dir_to_check, "(", i - 1, ")")
    }
    expect_equal(actual_dir, dir_to_check)
    dir.create(actual_dir)
  }
})
