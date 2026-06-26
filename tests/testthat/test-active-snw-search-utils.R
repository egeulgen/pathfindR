test_that("`.parse_experiment()` -- non-data-frame input is rejected", {
  expect_error(
    .parse_experiment(list(gene = "a", pvalue = 0.1)),
    "must be a data frame"
  )
  expect_error(.parse_experiment(1:10), "must be a data frame")
  expect_error(.parse_experiment(matrix(1:4, 2)), "must be a data frame")
})

test_that("`.parse_experiment()` -- named gene/pvalue columns are used case-insensitively and by position-independence", {
  # gene is column 2, pvalue is column 1 -> must still be picked by name
  df <- data.frame(
    PValue = c(0.01, 0.2),
    GENE = c("tp53", "egfr"),
    stringsAsFactors = FALSE
  )
  out <- .parse_experiment(df)
  expect_equal(out$gene, c("TP53", "EGFR"))
  expect_equal(out$pvalue, c(0.01, 0.2))
})

test_that("`.parse_experiment()` -- falls back to the first two columns when names are absent", {
  df <- data.frame(
    g = c("a", "b"),
    p = c(0.5, 0.6),
    extra = c(1, 2),
    stringsAsFactors = FALSE
  )
  out <- .parse_experiment(df)
  expect_equal(out$gene, c("A", "B"))
  expect_equal(out$pvalue, c(0.5, 0.6))
})

test_that("`.parse_experiment()` -- gene names are upper-cased and column types are correct", {
  df <- data.frame(
    gene = c("brca1", "brca2"),
    pvalue = c("0.001", "0.5"), # character p-values must be coerced to numeric
    stringsAsFactors = FALSE
  )
  out <- .parse_experiment(df)
  expect_type(out$gene, "character")
  expect_type(out$pvalue, "double")
  expect_equal(out$gene, c("BRCA1", "BRCA2"))
  expect_equal(out$pvalue, c(0.001, 0.5))
})

test_that("`.parse_experiment()` -- factor gene columns are coerced to character", {
  df <- data.frame(
    gene   = factor(c("foo", "bar")),
    pvalue = c(0.1, 0.2)
  )
  out <- .parse_experiment(df)
  expect_type(out$gene, "character")
  expect_equal(out$gene, c("FOO", "BAR"))
})

test_that("`.parse_experiment()` -- row count and ordering are preserved", {
  df <- data.frame(
    gene = c("a", "b", "c"),
    pvalue = c(0.3, 0.1, 0.9),
    stringsAsFactors = FALSE
  )
  out <- .parse_experiment(df)
  expect_equal(nrow(out), 3L)
  expect_equal(out$gene, c("A", "B", "C"))
  expect_equal(out$pvalue, c(0.3, 0.1, 0.9))
})

test_that("`.parse_experiment()` -- a single-row data frame works", {
  out <- .parse_experiment(
    data.frame(gene = "only", pvalue = 0.04, stringsAsFactors = FALSE)
  )
  expect_equal(out$gene, "ONLY")
  expect_equal(out$pvalue, 0.04)
})

test_that("`.parse_experiment()` -- fewer than two columns is an error", {
  df <- data.frame(gene = c("a", "b"), stringsAsFactors = FALSE)
  expect_error(.parse_experiment(df), "at least two columns")
})

test_that("`.parse_experiment()` -- only the named columns are kept; trailing columns are dropped", {
  df <- data.frame(
    gene = c("a", "b"),
    pvalue = c(0.1, 0.2),
    junk = c("x", "y"),
    stringsAsFactors = FALSE
  )
  out <- .parse_experiment(df)
  expect_named(out, c("gene", "pvalue"))
})


make_sif <- function(lines, env = parent.frame()) {
  path <- tempfile(fileext = ".sif")
  writeLines(lines, path)
  path
}

test_that("`build_network()` -- a 2-column SIF yields the expected undirected graph", {
  net <- build_network(make_sif(c("a b", "b c", "c a")))
  expect_setequal(net$nodes, c("A", "B", "C"))
  expect_equal(igraph::vcount(net$g), 3L)
  expect_equal(igraph::ecount(net$g), 3L)
  expect_false(igraph::is_directed(net$g))
})

test_that("`build_network()` -- a 3-column SIF uses columns 1 and 3 and ignores the middle column", {
  net <- build_network(make_sif(c("a pp b", "b pp c")))
  expect_setequal(net$nodes, c("A", "B", "C"))
  expect_equal(igraph::ecount(net$g), 2L)
  # the interaction-type token must not become a node
  expect_false("PP" %in% net$nodes)
})

test_that("`build_network()` -- node names are upper-cased", {
  net <- build_network(make_sif("foo bar"))
  expect_setequal(net$nodes, c("FOO", "BAR"))
})

test_that("`build_network()` -- tab-delimited files are parsed", {
  net <- build_network(make_sif(c("a\tb", "b\tc")))
  expect_setequal(net$nodes, c("A", "B", "C"))
  expect_equal(igraph::ecount(net$g), 2L)
})

test_that("`build_network()` -- self-interactions are discarded", {
  net <- build_network(make_sif(c("a a", "a b")))
  expect_setequal(net$nodes, c("A", "B"))
  expect_equal(igraph::ecount(net$g), 1L)
})

test_that("`build_network()` -- the column count is fixed by the first line; mismatched lines are dropped", {
  # first line is 2-col, so the 3-col line ("c d e") is discarded
  net <- build_network(make_sif(c("a b", "c d e", "b c")))
  expect_setequal(net$nodes, c("A", "B", "C"))
  expect_equal(igraph::ecount(net$g), 2L)
  expect_false(any(c("D", "E") %in% net$nodes))
})

test_that("`build_network()` -- parallel / reversed duplicate edges collapse to one graph edge", {
  net <- build_network(make_sif(c("a b", "b a", "a b")))
  expect_setequal(net$nodes, c("A", "B"))
  expect_equal(igraph::ecount(net$g), 1L)
})

test_that("`build_network()` -- a single-edge file works", {
  net <- build_network(make_sif("x y"))
  expect_setequal(net$nodes, c("X", "Y"))
  expect_equal(igraph::ecount(net$g), 1L)
  expect_length(net$nbr_idx, 2L)
})

test_that("`build_network()` -- an empty SIF errors", {
  expect_error(build_network(make_sif(character(0))), "empty")
  expect_error(build_network(make_sif(c("", ""))), "empty")
})

test_that("`build_network()` -- a first line that is not 2 or 3 columns errors", {
  expect_error(build_network(make_sif("a b c d")), "2 or 3 columns")
  expect_error(build_network(make_sif("a")), "2 or 3 columns")
})

test_that("`build_network()` -- a file with only self-interactions errors", {
  expect_error(
    build_network(make_sif(c("a a", "b b"))),
    "no non-self interactions"
  )
})

test_that("`build_network()` -- `name2id` maps each node to its position in nodes", {
  net <- build_network(make_sif(c("a b", "b c", "c d")))
  expect_equal(
    net$name2id[net$nodes],
    stats::setNames(seq_along(net$nodes), net$nodes)
  )
  expect_equal(net$nodes[net$name2id], net$nodes)
})

test_that("`build_network()` -- nbr and nbr_idx agree, and match the igraph adjacency", {
  net <- build_network(make_sif(c("a b", "b c", "c a", "c d")))
  for (nm in net$nodes) {
    i <- net$name2id[[nm]]
    # nbr_idx indexes into nodes to recover the neighbour names
    expect_equal(net$nodes[net$nbr_idx[[i]]], net$nbr[[i]])
    # neighbour set equals the graph's neighbour set (order-independent)
    expect_setequal(net$nbr[[i]], names(igraph::neighbors(net$g, nm)))
  }
})

test_that("`build_network()` -- neighbour lists are symmetric (undirected)", {
  net <- build_network(make_sif(c("a b", "b c", "c a", "c d")))
  for (nm in net$nodes) {
    i <- net$name2id[[nm]]
    for (j in net$nbr_idx[[i]]) {
      expect_true(i %in% net$nbr_idx[[j]])
    }
  }
})

test_that("`build_network()` -- nbr_idx contains no self-loops or duplicates", {
  net <- build_network(make_sif(c("a b", "b c", "c a", "c d")))
  for (i in seq_along(net$nodes)) {
    ids <- net$nbr_idx[[i]]
    expect_false(i %in% ids) # no self-loop
    expect_equal(anyDuplicated(ids), 0L) # no duplicate neighbours
  }
})

test_that("`build_network()` -- CSR offsets and neighbours are well-formed", {
  net <- build_network(make_sif(c("a b", "b c", "c a", "c d")))
  N <- length(net$nodes)

  expect_length(net$csr_offsets, N + 1L)
  expect_equal(net$csr_offsets[1L], 0L)
  expect_true(all(diff(net$csr_offsets) >= 0)) # non-decreasing
  expect_equal(
    diff(net$csr_offsets), # span == degree
    vapply(net$nbr_idx, length, integer(1))
  )
  expect_equal(
    net$csr_offsets[N + 1L], # 2 * |E|
    2L * igraph::ecount(net$g)
  )

  expect_length(net$csr_nbrs, net$csr_offsets[N + 1L])
  # csr_nbrs is the 0-based flattening of nbr_idx in node order
  expect_equal(
    net$csr_nbrs,
    as.integer(unlist(net$nbr_idx, use.names = FALSE)) - 1L
  )
  expect_true(all(net$csr_nbrs >= 0L & net$csr_nbrs < N))
})

test_that("`build_network()` -- each CSR slice reproduces that node's 0-based neighbour ids", {
  net <- build_network(make_sif(c("a b", "b c", "c a", "c d")))
  for (i in seq_along(net$nodes)) {
    lo <- net$csr_offsets[i] + 1L
    hi <- net$csr_offsets[i + 1L]
    slice <- if (hi >= lo) net$csr_nbrs[lo:hi] else integer(0)
    expect_equal(slice, net$nbr_idx[[i]] - 1L)
  }
})

test_that("`build_network()` -- a node with no neighbours yields empty nbr / nbr_idx entries", {
  testthat::local_mocked_bindings(
    java_node_order = function(src, tgt) c(unique(c(src, tgt)), "ISO")
  )

  net <- build_network(make_sif(c("a b", "b c")))
  expect_true("ISO" %in% net$nodes)
  i <- net$name2id[["ISO"]]

  # the isolated node takes the empty-neighbour branch
  expect_identical(net$nbr[[i]], character(0))
  expect_identical(net$nbr_idx[[i]], integer(0))

  # it is still a real, degree-0 vertex in the graph
  expect_true("ISO" %in% igraph::V(net$g)$name)
  expect_equal(unname(igraph::degree(net$g, "ISO")), 0)

  # and contributes an empty span to the CSR adjacency
  expect_equal(net$csr_offsets[i + 1L] - net$csr_offsets[i], 0L)

  # the real nodes are unaffected
  expect_true(all(c("A", "B", "C") %in% net$nodes))
  expect_equal(igraph::ecount(net$g), 2L)
})

test_that("`build_network()` -- the returned object exposes the documented elements", {
  net <- build_network(make_sif(c("a b", "b c")))
  expect_named(
    net,
    c("g", "nodes", "nbr", "nbr_idx", "name2id", "csr_offsets", "csr_nbrs"),
    ignore.order = TRUE
  )
  expect_s3_class(net$g, "igraph")
})
