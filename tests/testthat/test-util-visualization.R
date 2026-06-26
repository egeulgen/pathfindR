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
  expect_error(
    order_df_by_columnn(input_df, "non-existent"),
    "`order_by` column doesn't exist in `df`"
  )
})

# A minimal stand-in for what ggkegg::pathway() returns: an igraph graph with
# `name`, `type`, `x`, `y` vertex attributes.
fake_kegg_graph <- function(names = c("hsa:2821", "hsa:226", "hsa:229"),
                            types = NULL) {
  n <- length(names)
  if (is.null(types)) types <- rep("gene", n)
  g <- igraph::make_empty_graph(n = n, directed = FALSE)
  if (n >= 2) g <- igraph::add_edges(g, c(1, 2))
  igraph::V(g)$name <- names
  igraph::V(g)$type <- types
  igraph::V(g)$x <- seq_len(n)
  igraph::V(g)$y <- rep(0, n)
  g
}

skip_kegg <- function() {
  skip_on_cran()
  skip_if_offline()
  skip_if_not_installed("ggkegg")
  skip_if_not_installed("ggraph")
}

test_that("`color_kegg_pathway()` -- rejects a non-logical `scale_vals`", {
  change_vec <- c("hsa:2821" = -2, "hsa:226" = 4, "hsa:229" = 6)
  expect_error(
    color_kegg_pathway("hsa03040", change_vec, scale_vals = "INVALID"),
    "`scale_vals` should be logical"
  )
})

test_that("`color_kegg_pathway()` -- rejects a non-atomic `node_cols`", {
  change_vec <- c("hsa:2821" = -2, "hsa:226" = 4, "hsa:229" = 6)
  expect_error(
    color_kegg_pathway("hsa03040", change_vec, node_cols = list()),
    "`node_cols` should be a vector of colors"
  )
})

test_that("`color_kegg_pathway()` -- requires `node_cols` of length 3 when changes are supplied", {
  change_vec <- c("hsa:2821" = -2, "hsa:226" = 4, "hsa:229" = 6)
  expect_error(
    color_kegg_pathway("hsa03040", change_vec, node_cols = rep("red", 4)),
    "the length of `node_cols` should be 3"
  )
})

test_that("`color_kegg_pathway()` -- rejects invalid colors in `node_cols`", {
  change_vec <- c("hsa:2821" = -2, "hsa:226" = 4, "hsa:229" = 6)
  expect_error(
    color_kegg_pathway("hsa03040", change_vec,
      node_cols = c("red", "#FFFFFF", "INVALID")
    ),
    "`node_cols` should be a vector of valid colors"
  )
})

test_that("`color_kegg_pathway()` -- the length-3 check is bypassed when all changes are the 1e6 dummy", {
  skip_if_not_installed("ggkegg")
  const_vec <- c("hsa:2821" = 1e6, "hsa:226" = 1e6)

  with_mocked_bindings(
    {
      # length-1 node_cols is allowed here, only the palette message fires
      expect_message(
        color_kegg_pathway("hsa00010", const_vec, node_cols = "red"),
        "values are 1e6"
      )
      # invalid colors are still rejected even in the all-1e6 case
      expect_error(
        color_kegg_pathway("hsa00010", const_vec, node_cols = "NOTACOLOR"),
        "valid colors"
      )
    },
    pathway = function(...) NULL,
    .package = "ggkegg"
  )
})

test_that("`color_kegg_pathway()` -- returns NULL (with a message) when pathway parsing errors", {
  skip_if_not_installed("ggkegg")
  change_vec <- c("hsa:2821" = -2)

  with_mocked_bindings(
    {
      expect_message(
        res <- color_kegg_pathway("hsaXXXX", change_vec),
        "Cannot parse KEGG pathway"
      )
      expect_null(res)
    },
    pathway = function(...) stop("boom"),
    .package = "ggkegg"
  )
})

test_that("`color_kegg_pathway()` -- returns NULL (with a message) when pathway parsing warns", {
  skip_if_not_installed("ggkegg")
  change_vec <- c("hsa:2821" = -2)

  with_mocked_bindings(
    {
      expect_message(
        res <- color_kegg_pathway("hsaXXXX", change_vec),
        "Cannot parse KEGG pathway"
      )
      expect_null(res)
    },
    pathway = function(...) warning("watch out"),
    .package = "ggkegg"
  )
})

test_that("`color_kegg_pathway()` -- returns NULL silently when the graph object is NULL", {
  skip_if_not_installed("ggkegg")
  change_vec <- c("hsa:2821" = -2, "hsa:226" = 4, "hsa:229" = 6)

  with_mocked_bindings(
    {
      expect_silent(res <- color_kegg_pathway("hsa00010", change_vec))
      expect_null(res)
    },
    pathway = function(...) NULL,
    .package = "ggkegg"
  )
})

test_that("`color_kegg_pathway()` -- returns NULL when no input genes map to the pathway", {
  skip_if_not_installed("ggkegg")
  change_vec <- c("hsa:2821" = -2, "hsa:226" = 4)
  g <- fake_kegg_graph(names = c("hsa:0001", "hsa:0002", "hsa:0003"))

  with_mocked_bindings(
    {
      expect_null(color_kegg_pathway("hsa00010", change_vec))
    },
    pathway = function(...) g,
    .package = "ggkegg"
  )
})

# The ggkegg geoms / overlay are replaced with no-op layers so the assembly runs
# offline; we assert that a ggraph object comes back across the palette and
# scaling branches.

test_that("`color_kegg_pathway()` -- builds a ggraph for matching genes, across palette/scaling branches", {
  skip_if_not_installed("ggkegg")
  skip_if_not_installed("ggraph")
  skip_if_not_installed("ggplot2")

  change_vec <- c("hsa:2821" = -2, "hsa:226" = 4, "hsa:229" = 6)
  g <- fake_kegg_graph()

  with_mocked_bindings(
    {
      # default palette (red/gray/green), scaled
      expect_s3_class(color_kegg_pathway("hsa00010", change_vec), "ggraph")
      # default palette, unscaled
      expect_s3_class(
        color_kegg_pathway("hsa00010", change_vec, scale_vals = FALSE), "ggraph"
      )
      # user-supplied 3-colour palette
      expect_s3_class(
        color_kegg_pathway("hsa00010", change_vec,
          node_cols = c("red", "gray", "blue")
        ),
        "ggraph"
      )
    },
    pathway = function(...) g,
    geom_node_rect = function(...) ggplot2::geom_blank(),
    overlay_raw_map = function(...) ggplot2::geom_blank(),
    .package = "ggkegg"
  )
})

test_that("`color_kegg_pathway()` -- all-1e6 changes use a single colour and emit a message", {
  skip_if_not_installed("ggkegg")
  skip_if_not_installed("ggraph")
  skip_if_not_installed("ggplot2")

  const_vec <- c("hsa:2821" = 1e6, "hsa:226" = 1e6, "hsa:229" = 1e6)
  g <- fake_kegg_graph()

  with_mocked_bindings(
    {
      # node_cols supplied -> message, first colour used
      expect_message(
        p <- color_kegg_pathway("hsa00010", const_vec,
          node_cols = c("red", "gray", "blue")
        ),
        "values are 1e6"
      )
      expect_s3_class(p, "ggraph")

      # node_cols NULL -> default single colour, no palette message
      expect_no_message(
        p2 <- color_kegg_pathway("hsa00010", const_vec)
      )
      expect_s3_class(p2, "ggraph")
    },
    pathway = function(...) g,
    geom_node_rect = function(...) ggplot2::geom_blank(),
    overlay_raw_map = function(...) ggplot2::geom_blank(),
    .package = "ggkegg"
  )
})

test_that("`color_kegg_pathway()` -- integration: colours a real pathway and returns a ggraph", {
  skip_kegg()
  change_vec <- c("hsa:2821" = -2, "hsa:226" = 4, "hsa:229" = 6)
  expect_s3_class(color_kegg_pathway("hsa00010", change_vec), "ggraph")
})

test_that("`color_kegg_pathway()` -- integration: returns NULL when no genes match the real pathway", {
  skip_kegg()
  change_vec <- c(-2, 4, 6)
  names(change_vec) <- rep("missing", 3)
  expect_null(color_kegg_pathway("hsa00010", change_vec))
})
