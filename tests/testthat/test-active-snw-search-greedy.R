# Triangle A-B-C; nbr_idx is 1-based, aligned to `nodes`
make_network <- function(nbr_idx = list(c(2L, 3L), c(1L, 3L), c(1L, 2L))) {
  list(nodes = c("A", "B", "C"), nbr_idx = nbr_idx)
}

# z stored out of node order to verify re-alignment to network$nodes
make_sc <- function() {
  list(
    z     = c(C = 3, A = 1, B = 2),
    means = c(0, 0, 0),
    stds  = c(1, 1, 1),
    nodes = c("A", "B", "C")
  )
}

# Defaults use doubles so the integer coercions are exercised
make_params <- function(...) {
  defaults <- list(
    gr_max_depth         = 1,
    gr_search_depth      = 2,
    gr_overlap_threshold = 0.5,
    gr_subnetwork_num    = 1000
  )
  utils::modifyList(defaults, list(...))
}

# ---- orchestration (mocked C++) ---------------------------------------------

test_that("`.greedy_search()` -- an empty network returns an empty list without entering C++", {
  with_mocked_bindings(
    {
      net <- list(nodes = character(0), nbr_idx = list())
      expect_equal(.greedy_search(net, make_sc(), make_params()), list())
    },
    run_greedy_search = function(...) stop("should not be called")
  )
})

test_that("`.greedy_search()` -- errors if network$nbr_idx is missing", {
  with_mocked_bindings(
    {
      net <- list(nodes = c("A", "B"), nbr_idx = NULL)
      expect_error(
        .greedy_search(net, make_sc(), make_params()),
        "nbr_idx is missing"
      )
    },
    run_greedy_search = function(...) stop("should not be called")
  )
})

test_that("`.greedy_search()` -- marshals z (in node order), means/stds, nbr_idx, names and coerced params into the C++ call", {
  captured <- new.env(parent = emptyenv())
  net <- make_network()
  sc <- make_sc()
  params <- make_params(gr_max_depth = 2, gr_search_depth = 3, gr_subnetwork_num = 50)

  with_mocked_bindings(
    {
      .greedy_search(net, sc, params)
    },
    run_greedy_search = function(nbr_idx, z_vec, sc_means, sc_stds, node_names,
                                 max_depth, search_depth, n_nodes,
                                 overlap_threshold, subnetwork_num) {
      captured$nbr_idx <- nbr_idx
      captured$z_vec <- z_vec
      captured$sc_means <- sc_means
      captured$sc_stds <- sc_stds
      captured$node_names <- node_names
      captured$max_depth <- max_depth
      captured$search_depth <- search_depth
      captured$n_nodes <- n_nodes
      captured$overlap_threshold <- overlap_threshold
      captured$subnetwork_num <- subnetwork_num
      list()
    }
  )

  expect_identical(captured$nbr_idx, net$nbr_idx) # forwarded verbatim
  expect_equal(captured$z_vec, c(1, 2, 3)) # realigned to A, B, C
  expect_equal(captured$sc_means, c(0, 0, 0))
  expect_equal(captured$sc_stds, c(1, 1, 1))
  expect_identical(captured$node_names, c("A", "B", "C"))
  expect_identical(captured$n_nodes, 3L)

  expect_identical(captured$max_depth, 2L) # double -> integer
  expect_identical(captured$search_depth, 3L) # double -> integer
  expect_identical(captured$subnetwork_num, 50L) # double -> integer
  expect_identical(captured$overlap_threshold, 0.5)
})

test_that("`.greedy_search()` -- returns the C++ candidates unchanged", {
  canned <- list(
    list(nodes = c("A", "B", "C"), score = 5),
    list(nodes = "A", score = 0)
  )
  with_mocked_bindings(
    {
      expect_identical(.greedy_search(make_network(), make_sc(), make_params()), canned)
    },
    run_greedy_search = function(...) canned
  )
})

test_that("`.greedy_search()` -- emits a progress message only when verbose", {
  with_mocked_bindings(
    {
      expect_message(
        .greedy_search(make_network(), make_sc(), make_params(), verbose = TRUE),
        "Running greedy search"
      )
      expect_no_message(
        .greedy_search(make_network(), make_sc(), make_params(), verbose = FALSE)
      )
    },
    run_greedy_search = function(...) list()
  )
})

# ---- integration: real C++ on a toy graph -----------------------------------

# Build a real network from a SIF and a positive-z (mean 0, sd 1) score context
build_toy <- function(sif_lines, neg = NULL, env = parent.frame()) {
  sif <- tempfile(fileext = ".sif")
  writeLines(sif_lines, sif)
  network <- build_network(sif)
  nm <- network$nodes
  z <- stats::setNames(rep(3, length(nm)), nm)
  if (!is.null(neg)) z[names(neg)] <- neg
  sc <- list(z = z, means = rep(0, length(nm)), stds = rep(1, length(nm)), nodes = nm)
  list(network = network, sc = sc)
}

toy_params <- function() {
  list(
    gr_max_depth         = 1L,
    gr_search_depth      = 1L,
    gr_overlap_threshold = 0.5,
    gr_subnetwork_num    = 1000
  )
}

test_that("`.greedy_search()` -- integration: greedily recovers the full positive module (real C++)", {
  toy <- build_toy(c("A B", "B C", "C A"))
  res <- .greedy_search(toy$network, toy$sc, toy_params())

  expect_true(length(res) >= 1L)

  node_sets <- lapply(res, function(s) sort(s$nodes))
  expect_true(any(vapply(node_sets, function(ns) setequal(ns, c("A", "B", "C")), logical(1))))

  top <- res[[which.max(vapply(res, function(s) s$score, numeric(1)))]]
  expect_setequal(top$nodes, c("A", "B", "C"))
  expect_true(top$score > 0)
})

test_that("`.greedy_search()` -- integration: is deterministic for identical inputs (real C++)", {
  toy <- build_toy(c("A B", "B C", "C A"))
  r1 <- .greedy_search(toy$network, toy$sc, toy_params())
  r2 <- .greedy_search(toy$network, toy$sc, toy_params())

  nodes_of <- function(r) lapply(r, function(s) sort(s$nodes))
  expect_equal(nodes_of(r1), nodes_of(r2))
  expect_equal(
    vapply(r1, function(s) s$score, numeric(1)),
    vapply(r2, function(s) s$score, numeric(1))
  )
})

test_that("`.greedy_search()` -- integration: excludes a strongly negative neighbour from the module (real C++)", {
  # D hangs off C but has a strongly negative z, so adding it never improves the score
  toy <- build_toy(c("A B", "B C", "C A", "C D"), neg = c(D = -10))
  res <- .greedy_search(toy$network, toy$sc, toy_params())

  with_A <- Filter(function(s) "A" %in% s$nodes, res)
  expect_true(length(with_A) >= 1L)
  for (s in with_A) {
    expect_setequal(s$nodes, c("A", "B", "C")) # D is not absorbed into the module
  }
})
