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

test_that("`order_df_by_columnn()` -- identifies if column can be ordered and orders data frame", {
  input_df <- data.frame(
    values = c(5, 1, 10),
    floats = c(0.5, 0.3, 0.2),
    chars = c("second", "first", "third"),
    incomplete = c(1, NA, 4),
    bad = I(list(c(1, 2), c(3, 4), c(5, 6)))
  )

  expect_is(out <- order_df_by_columnn(input_df, "values"), "data.frame")
  expect_equal(out$values, c(1, 5, 10))

  expect_is(out <- order_df_by_columnn(input_df, "floats"), "data.frame")
  expect_equal(out$floats, c(0.2, 0.3, 0.5))

  expect_is(out <- order_df_by_columnn(input_df, "chars"), "data.frame")
  expect_equal(out$chars, c("first", "second", "third"))

  expect_error(
    order_df_by_columnn(input_df, "incomplete"),
    "Column values of `order_by` cannot have NAs!"
  )
  expect_error(
    order_df_by_columnn(input_df, "bad"),
    "`order_by` \\(bad\\) cannot be used to order the `df`:"
  )
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

  expect_error(color_kegg_pathway(
    pw_id = "hsa03040", change_vec = change_vec,
    scale_vals = "INVALID"
  ), "`scale_vals` should be logical")
  expect_error(color_kegg_pathway(
    pw_id = "hsa03040", change_vec = change_vec,
    node_cols = list()
  ), "`node_cols` should be a vector of colors")
  expect_error(color_kegg_pathway(
    pw_id = "hsa03040", change_vec = change_vec,
    node_cols = rep("red", 4)
  ), "the length of `node_cols` should be 3")
  expect_error(color_kegg_pathway(
    pw_id = "hsa03040", change_vec = change_vec,
    node_cols = c("red", "#FFFFFF", "INVALID")
  ), "`node_cols` should be a vector of valid colors")

  skip_on_cran()

  constant_vec <- rep(1e+06, 3)
  names(constant_vec) <- c("hsa:2821", "hsa:226", "hsa:229")

  expect_silent(color_kegg_pathway(
    pw_id = "hsa03040", change_vec = change_vec,
    node_cols = c("red", "blue", "green")
  ))
  expect_message(color_kegg_pathway(
    pw_id = "hsa03040", change_vec = constant_vec,
    node_cols = c("red", "blue", "green")
  ))

  expect_null(suppressWarnings(color_kegg_pathway(pw_id = "hsa03040", change_vec = NULL)))
  expect_message(color_kegg_pathway(pw_id = "hsa11111", change_vec = c()))
})
