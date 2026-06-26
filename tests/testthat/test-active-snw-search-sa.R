# Graph: A-B-C is one path component; D is isolated
make_network <- function() {
  g <- igraph::graph_from_literal(A - B, B - C, D)
  list(
    nodes       = c("A", "B", "C", "D"),
    g           = g,
    csr_offsets = c(0L, 1L, 3L, 4L, 4L), # plausible CSR; only forwarded to C++
    csr_nbrs    = c(1L, 0L, 2L, 1L)
  )
}

# z deliberately stored out of node order to test re-alignment to network$nodes
make_sc <- function() {
  list(
    z     = c(D = 5, A = 2, C = 3, B = 1),
    means = rep(0, 4),
    stds  = rep(1, 4),
    nodes = c("A", "B", "C", "D")
  )
}

make_params <- function(...) {
  defaults <- list(
    gene_init_prob           = 0.2,
    start_with_all_positives = TRUE,
    sa_initial_temp          = 1,
    sa_final_temp            = 0.01,
    sa_iterations            = 100, # double on purpose -> must be coerced to int
    seed                     = 1234 # double on purpose -> must be coerced to int
  )
  utils::modifyList(defaults, list(...))
}

test_that("`.simulated_annealing()` -- an empty network returns an empty list without entering C++", {
  with_mocked_bindings(
    {
      net <- make_network()
      net$nodes <- character(0)
      expect_equal(.simulated_annealing(net, make_sc(), make_params()), list())
    },
    run_simulated_annealing = function(...) stop("should not be called")
  )
})

test_that("`.simulated_annealing()` -- marshals z (in node order), means/stds, CSR and params into the C++ call", {
  captured <- new.env(parent = emptyenv())
  net <- make_network()
  sc <- make_sc()
  params <- make_params()

  with_mocked_bindings(
    {
      .simulated_annealing(net, sc, params)
    },
    run_simulated_annealing = function(csr_offsets, csr_nbrs, z, means, stds,
                                       n_nodes, gene_init_prob,
                                       start_with_all_positives,
                                       sa_initial_temp, sa_final_temp,
                                       sa_iterations, seed) {
      captured$csr_offsets <- csr_offsets
      captured$csr_nbrs <- csr_nbrs
      captured$z <- z
      captured$means <- means
      captured$stds <- stds
      captured$n_nodes <- n_nodes
      captured$gene_init_prob <- gene_init_prob
      captured$start_with_all_positives <- start_with_all_positives
      captured$sa_initial_temp <- sa_initial_temp
      captured$sa_final_temp <- sa_final_temp
      captured$sa_iterations <- sa_iterations
      captured$seed <- seed
      rep(FALSE, n_nodes)
    }
  )

  # z is realigned to network$nodes order (A, B, C, D)
  expect_equal(captured$z, c(2, 1, 3, 5))
  expect_equal(captured$means, c(0, 0, 0, 0))
  expect_equal(captured$stds, c(1, 1, 1, 1))
  expect_identical(captured$n_nodes, 4L)

  # CSR forwarded verbatim
  expect_identical(captured$csr_offsets, net$csr_offsets)
  expect_identical(captured$csr_nbrs, net$csr_nbrs)

  # scalar params, with the documented type coercions
  expect_identical(captured$gene_init_prob, 0.2)
  expect_identical(captured$start_with_all_positives, TRUE)
  expect_identical(captured$sa_initial_temp, 1)
  expect_identical(captured$sa_final_temp, 0.01)
  expect_identical(captured$sa_iterations, 100L) # double -> integer
  expect_identical(captured$seed, 1234L) # double -> integer
})

# Helper: run the search with a given `start_with_all_positives` value and
# return whatever was actually forwarded to the C++ routine.
forwarded_swap_flag <- function(value, present = TRUE) {
  captured <- new.env(parent = emptyenv())
  params <- make_params()
  if (present) {
    params$start_with_all_positives <- value
  } else {
    params$start_with_all_positives <- NULL
  }
  with_mocked_bindings(
    {
      .simulated_annealing(make_network(), make_sc(), params)
    },
    run_simulated_annealing = function(..., start_with_all_positives) {
      captured$flag <- start_with_all_positives
      rep(FALSE, 4)
    }
  )
  captured$flag
}

test_that("`.simulated_annealing()` -- start_with_all_positives = TRUE is forwarded as TRUE", {
  expect_identical(forwarded_swap_flag(TRUE), TRUE)
})

test_that("`.simulated_annealing()` -- start_with_all_positives = FALSE is forwarded as FALSE", {
  expect_identical(forwarded_swap_flag(FALSE), FALSE)
})

test_that("`.simulated_annealing()` -- non-TRUE start_with_all_positives values collapse to FALSE (isTRUE semantics)", {
  # isTRUE() only accepts a length-1 logical TRUE; everything else is FALSE.
  expect_identical(forwarded_swap_flag(1), FALSE) # numeric 1
  expect_identical(forwarded_swap_flag("TRUE"), FALSE) # character
  expect_identical(forwarded_swap_flag(NA), FALSE) # logical NA
  expect_identical(forwarded_swap_flag(c(TRUE, TRUE)), FALSE) # length > 1
})

test_that("`.simulated_annealing()` -- a missing start_with_all_positives defaults to FALSE", {
  # params$start_with_all_positives is NULL -> isTRUE(NULL) is FALSE, no error.
  expect_identical(forwarded_swap_flag(NULL, present = FALSE), FALSE)
})

test_that("`.simulated_annealing()` -- the flag is always a single logical, never the raw input", {
  for (val in list(TRUE, FALSE, 1, 0, "TRUE", NA, c(TRUE, FALSE))) {
    flag <- forwarded_swap_flag(val)
    expect_type(flag, "logical")
    expect_length(flag, 1L)
    expect_false(is.na(flag))
  }
})

test_that("`.simulated_annealing()` -- the returned mask selects on-nodes whose subnetworks are returned", {
  with_mocked_bindings(
    {
      result <- .simulated_annealing(make_network(), make_sc(), make_params())
      expect_equal(length(result), 1L)
      expect_setequal(result[[1]]$nodes, c("A", "B", "C"))
      expect_true(result[[1]]$score > 0)
    },
    # integer mask (A,B,C on, D off) -> tests as.logical() coercion too
    run_simulated_annealing = function(...) c(1L, 1L, 1L, 0L)
  )
})

test_that("`.simulated_annealing()` -- multiple components come back sorted by score descending", {
  with_mocked_bindings(
    {
      result <- .simulated_annealing(make_network(), make_sc(), make_params())
      expect_equal(length(result), 2L)
      scores <- vapply(result, function(s) s$score, numeric(1))
      expect_true(scores[1] >= scores[2])
      expect_setequal(result[[1]]$nodes, c("A", "B", "C")) # higher-scoring first
      expect_equal(result[[2]]$nodes, "D") # isolated singleton, score 0
    },
    run_simulated_annealing = function(...) rep(TRUE, 4) # everything on
  )
})

test_that("`.simulated_annealing()` -- an all-off mask yields no subnetworks", {
  with_mocked_bindings(
    {
      result <- .simulated_annealing(make_network(), make_sc(), make_params())
      expect_equal(result, list())
    },
    run_simulated_annealing = function(...) rep(FALSE, 4)
  )
})

test_that("`.simulated_annealing()` -- progress messages are emitted only when verbose", {
  with_mocked_bindings(
    {
      expect_message(
        .simulated_annealing(make_network(), make_sc(), make_params(), verbose = TRUE),
        "Running simulated annealing"
      )
      expect_no_message(
        .simulated_annealing(make_network(), make_sc(), make_params(), verbose = FALSE)
      )
    },
    run_simulated_annealing = function(...) c(1L, 1L, 1L, 0L)
  )
})

# ---- integration: real C++ run on a toy case --------------------------------
#
# These tests do NOT mock run_simulated_annealing -- they exercise the real C++
# loop. They additionally need `build_network()` (and its java_* helpers).
# They are deterministic by construction: a connected triangle whose
# nodes all have positive z-scores has a single global optimum (the whole
# triangle), and with the temperature frozen at a tiny value no worsening move
# is ever accepted, so once the solution is at the optimum it cannot leave it.
# That makes the outcome independent of which nodes the RNG happens to toggle.

# Build a real triangle network (A-B-C) and a positive-z score context for it.
make_toy <- function(env = parent.frame()) {
  sif <- tempfile(fileext = ".sif")
  writeLines(c("A B", "B C", "C A"), sif)
  network <- build_network(sif)
  nm <- network$nodes
  sc <- list(
    z     = stats::setNames(rep(3, length(nm)), nm), # all positive
    means = rep(0, length(nm)), # mean 0, sd 1 -> score == raw
    stds  = rep(1, length(nm)),
    nodes = nm
  )
  list(network = network, sc = sc)
}

# Frozen temperature => greedy-with-no-uphill: equal initial and final temps
# means a constant temperature, and a tiny value makes exp(-delta/T) underflow
# to 0 for any worsening move, so they are always rejected.
frozen_params <- function(...) {
  defaults <- list(
    start_with_all_positives = TRUE,
    gene_init_prob           = 0,
    sa_initial_temp          = 1e-12,
    sa_final_temp            = 1e-12,
    sa_iterations            = 200L,
    seed                     = 1234L
  )
  utils::modifyList(defaults, list(...))
}

test_that("`.simulated_annealing()` -- start_with_all_positives=TRUE retains the full optimal module (real C++)", {
  toy <- make_toy()
  res <- .simulated_annealing(
    toy$network, toy$sc,
    frozen_params(start_with_all_positives = TRUE)
  )

  expect_equal(length(res), 1L)
  expect_setequal(res[[1]]$nodes, c("A", "B", "C"))
  expect_true(res[[1]]$score > 0)
})

test_that("`.simulated_annealing()` -- start_with_all_positives=FALSE with gene_init_prob=1 starts all-on and reaches the same module (real C++)", {
  toy <- make_toy()
  res <- .simulated_annealing(
    toy$network, toy$sc,
    frozen_params(start_with_all_positives = FALSE, gene_init_prob = 1)
  )

  expect_equal(length(res), 1L)
  expect_setequal(res[[1]]$nodes, c("A", "B", "C"))
  expect_true(res[[1]]$score > 0)
})

test_that("`.simulated_annealing()` -- TRUE and all-on FALSE converge to the same frozen optimum (real C++)", {
  toy <- make_toy()
  res_true <- .simulated_annealing(
    toy$network, toy$sc, frozen_params(start_with_all_positives = TRUE)
  )
  res_false <- .simulated_annealing(
    toy$network, toy$sc,
    frozen_params(start_with_all_positives = FALSE, gene_init_prob = 1)
  )

  expect_setequal(res_true[[1]]$nodes, res_false[[1]]$nodes)
  expect_equal(res_true[[1]]$score, res_false[[1]]$score)
})

test_that("`.simulated_annealing()` -- a stochastic FALSE initialization is reproducible for a fixed seed (real C++)", {
  toy <- make_toy()
  params <- list(
    start_with_all_positives = FALSE,
    gene_init_prob           = 0.5,
    sa_initial_temp          = 1.0,
    sa_final_temp            = 0.01,
    sa_iterations            = 500L,
    seed                     = 1234L
  )

  r1 <- .simulated_annealing(toy$network, toy$sc, params)
  r2 <- .simulated_annealing(toy$network, toy$sc, params)

  nodes_of <- function(r) lapply(r, function(s) sort(s$nodes))
  expect_equal(nodes_of(r1), nodes_of(r2))
  expect_equal(
    vapply(r1, function(s) s$score, numeric(1)),
    vapply(r2, function(s) s$score, numeric(1))
  )
})
