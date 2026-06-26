# A graph with two connected components (A-B-C and D-E) plus an isolated node Z.
make_network <- function() {
  g <- igraph::graph_from_literal(A - B, B - C, D - E, Z)
  list(g = g, nodes = c("A", "B", "C", "D", "E", "Z"))
}

# A score context with trivial calibration (mean 0, sd 1 -> normalized == raw).
make_sc <- function() {
  z <- c(A = 1, B = 2, C = 3, D = 4, E = 5, Z = 6)
  list(
    z     = z,
    means = rep(0, length(z)),
    stds  = rep(1, length(z)),
    nodes = names(z)
  )
}

# Normalize a list-of-components to a canonical form for set comparison.
norm_comps <- function(x) {
  x <- lapply(x, sort)
  x[order(vapply(x, function(v) v[[1]], character(1)))]
}

test_that("`build_score_context()` -- keeps the smallest p-value, defaults non-significant nodes, calibrates in node order", {
  captured <- new.env(parent = emptyenv())
  testthat::local_mocked_bindings(
    get_java_zscores = function(p) p, # identity: inspect pmap directly via z
    get_java_mc_calibration = function(z, trials, seed) {
      captured$z <- z
      captured$trials <- trials
      captured$seed <- seed
      list(means = seq_along(z) * 10, stds = seq_along(z))
    }
  )

  network <- list(nodes = c("A", "B", "C"))
  experiment <- data.frame(
    gene = c("a", "b", "a"), # 'a' appears twice; lower-case -> upper-cased
    pvalue = c(0.01, 0.5, 0.001),
    stringsAsFactors = FALSE
  )
  params <- list(seed = 1234, p_for_nonsignificant = 0.7)

  sc <- build_score_context(network, experiment, params)

  # z is named over network$nodes, in node order
  expect_equal(names(sc$z), c("A", "B", "C"))
  expect_equal(unname(sc$z["A"]), 0.001) # smallest p-value for duplicated gene
  expect_equal(unname(sc$z["B"]), 0.5)
  expect_equal(unname(sc$z["C"]), 0.7) # absent from experiment -> default

  # MC replica is driven with 2000 trials, the run seed, and the node-ordered vector
  expect_equal(captured$trials, 2000)
  expect_equal(captured$seed, 1234)
  expect_equal(captured$z, c(0.001, 0.5, 0.7))

  # means / stds are passed straight through; nodes echoed
  expect_equal(sc$means, c(10, 20, 30))
  expect_equal(sc$stds, c(1, 2, 3))
  expect_equal(sc$nodes, c("A", "B", "C"))
})

test_that("`build_score_context()` -- clamps p-values into the safe (0, 1) range", {
  testthat::local_mocked_bindings(
    get_java_zscores = function(p) p,
    get_java_mc_calibration = function(z, trials, seed) {
      list(means = numeric(length(z)), stds = rep(1, length(z)))
    }
  )
  network <- list(nodes = c("A", "B"))
  experiment <- data.frame(
    gene = c("A", "B"), pvalue = c(0, 1), stringsAsFactors = FALSE
  )
  sc <- build_score_context(
    network, experiment, list(seed = 1, p_for_nonsignificant = 0.5)
  )
  expect_equal(unname(sc$z["A"]), 1e-13) # 0 -> MIN_SIG
  expect_equal(unname(sc$z["B"]), 1 - 1e-13) # 1 -> MAX_SIG
})

test_that("`build_score_context()` -- drops genes absent from the network and defaults their nodes", {
  testthat::local_mocked_bindings(
    get_java_zscores = function(p) p,
    get_java_mc_calibration = function(z, trials, seed) {
      list(means = numeric(length(z)), stds = rep(1, length(z)))
    }
  )
  network <- list(nodes = c("A", "B"))
  experiment <- data.frame(
    gene = c("A", "ZZZ"), pvalue = c(0.02, 0.04), stringsAsFactors = FALSE
  )
  sc <- build_score_context(
    network, experiment, list(seed = 1, p_for_nonsignificant = 0.9)
  )
  expect_equal(unname(sc$z["A"]), 0.02)
  expect_equal(unname(sc$z["B"]), 0.9) # no p-value -> default
  expect_equal(length(sc$z), 2L) # foreign gene did not leak in
})


test_that("`.score_subnetwork()` -- a single-node subnetwork always scores 0", {
  sc <- list(means = c(99, 99), stds = c(99, 99))
  expect_equal(.score_subnetwork(sc, 1L, 12.34, TRUE), 0)
  expect_equal(.score_subnetwork(sc, 1L, 12.34, FALSE), 0)
})

test_that("`.score_subnetwork()` -- the raw score is zsum / sqrt(n) when not normalizing", {
  sc <- list(means = c(NA, 1, 2, 3), stds = c(NA, 0.5, 1, 2))
  expect_equal(.score_subnetwork(sc, 4L, 8, FALSE), 8 / sqrt(4))
  expect_equal(.score_subnetwork(sc, 2L, 5, FALSE), 5 / sqrt(2))
})

test_that("`.score_subnetwork()` -- normalizing subtracts the MC mean and divides by the MC sd at size n", {
  sc <- list(means = c(NA, 1, 2, 3), stds = c(NA, 0.5, 1, 2))
  raw <- 8 / sqrt(4)
  expect_equal(.score_subnetwork(sc, 4L, 8, TRUE), (raw - sc$means[4]) / sc$stds[4])
})

test_that("`.score_subnetwork()` -- only a literal TRUE triggers normalization (isTRUE semantics)", {
  sc <- list(means = c(NA, 5, 5), stds = c(NA, 1, 1))
  raw <- 6 / sqrt(2)
  expect_equal(.score_subnetwork(sc, 2L, 6, NA), raw)
  expect_equal(.score_subnetwork(sc, 2L, 6, FALSE), raw)
})


test_that("`.make_subnetwork()` -- an empty node set is the zero subnetwork", {
  out <- .make_subnetwork(make_sc(), character(0))
  expect_equal(out$nodes, character(0))
  expect_equal(out$zsum, 0)
  expect_equal(out$score, 0)
})

test_that("`.make_subnetwork()` -- zsum sums member z-scores and score is the normalized score", {
  sc <- list(z = c(A = 1, B = 2, C = 3), means = c(NA, 1, 2), stds = c(NA, 2, 3))
  out <- .make_subnetwork(sc, c("A", "B"))
  expect_equal(out$nodes, c("A", "B"))
  expect_equal(out$zsum, 3)
  expect_equal(out$score, .score_subnetwork(sc, 2L, 3, TRUE))
})

test_that("`.make_subnetwork()` -- a single-node subnetwork carries its z-sum but scores 0", {
  sc <- list(z = c(A = 7), means = 0, stds = 1)
  out <- .make_subnetwork(sc, "A")
  expect_equal(out$zsum, 7)
  expect_equal(out$score, 0)
})

test_that("`.find_components_named()` -- no on-nodes yields no components", {
  expect_equal(.find_components_named(make_network(), character(0)), list())
})

test_that("`.find_components_named()` -- connected on-nodes form one component", {
  res <- .find_components_named(make_network(), c("A", "B", "C"))
  expect_equal(length(res), 1L)
  expect_setequal(res[[1]], c("A", "B", "C"))
})

test_that("`.find_components_named()` -- disconnected on-nodes split into separate components", {
  res <- norm_comps(.find_components_named(make_network(), c("A", "B", "C", "D", "E")))
  expect_equal(res, list(c("A", "B", "C"), c("D", "E")))
})

test_that("`.find_components_named()` -- on-nodes with no induced edges are singleton components", {
  # A-C is not an edge (path is A-B-C), so with only A and C on they are separate
  res <- norm_comps(.find_components_named(make_network(), c("A", "C")))
  expect_equal(res, list("A", "C"))
})

test_that("`.find_components_named()` -- an isolated node is its own component", {
  expect_equal(.find_components_named(make_network(), "Z"), list("Z"))
})


test_that("`.find_subnetworks()` -- returns one scored subnetwork per component", {
  net <- make_network()
  sc <- make_sc()
  subs <- .find_subnetworks(net, sc, c("A", "B", "C", "D", "E"))
  expect_equal(length(subs), 2L)

  sizes <- vapply(subs, function(s) length(s$nodes), integer(1))
  abc <- subs[[which(sizes == 3L)]]
  de <- subs[[which(sizes == 2L)]]

  expect_setequal(abc$nodes, c("A", "B", "C"))
  expect_equal(abc$zsum, sum(sc$z[c("A", "B", "C")]))
  expect_equal(abc$score, .score_subnetwork(sc, 3L, abc$zsum, TRUE))

  expect_setequal(de$nodes, c("D", "E"))
  expect_equal(de$zsum, sum(sc$z[c("D", "E")]))
  expect_equal(de$score, .score_subnetwork(sc, 2L, de$zsum, TRUE))
})

test_that("`.find_subnetworks()` -- no on-nodes yields no subnetworks", {
  expect_equal(.find_subnetworks(make_network(), make_sc(), character(0)), list())
})


test_that("`.sort_subnetworks_desc()` -- a list of 0 or 1 subnetworks is returned unchanged", {
  expect_equal(.sort_subnetworks_desc(list()), list())
  one <- list(list(nodes = "A", score = 3))
  expect_equal(.sort_subnetworks_desc(one), one)
})

test_that("`.sort_subnetworks_desc()` -- subnetworks are ordered by score, highest first", {
  subs <- list(
    list(nodes = "low", score = 1),
    list(nodes = "high", score = 5),
    list(nodes = "mid", score = 3)
  )
  out <- .sort_subnetworks_desc(subs)
  expect_equal(vapply(out, function(s) s$score, numeric(1)), c(5, 3, 1))
  expect_equal(out[[1]]$nodes, "high") # objects preserved, just reordered
})
