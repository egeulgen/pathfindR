# Build a GA individual literal (sort/compare only read $scores; $rep is a tag)
gi <- function(scores, id = NULL) list(rep = id, scores = scores)

test_that("`.ga_score()` -- is the top component score, or 0 when there are none", {
  expect_equal(.ga_score(gi(numeric(0))), 0)
  expect_equal(.ga_score(gi(c(5, 3, 1))), 5)
})


test_that("`.ga_compare()` -- the first strict rank difference decides", {
  expect_equal(.ga_compare(gi(c(5, 3)), gi(c(4, 9))), 1L) # 5 > 4
  expect_equal(.ga_compare(gi(c(4, 9)), gi(c(5, 3))), -1L) # 4 < 5
  expect_equal(.ga_compare(gi(c(5, 3)), gi(c(5, 2))), 1L) # tie, then 3 > 2
})

test_that("`.ga_compare()` -- when one score vector is a prefix of the other, more subnetworks wins", {
  expect_equal(.ga_compare(gi(c(5, 1)), gi(c(5))), 1L)
  expect_equal(.ga_compare(gi(c(5)), gi(c(5, 1))), -1L)
})

test_that("`.ga_compare()` -- equivalent individuals compare equal", {
  expect_equal(.ga_compare(gi(c(5, 3)), gi(c(5, 3))), 0L)
  expect_equal(.ga_compare(gi(numeric(0)), gi(numeric(0))), 0L)
})


test_that("`.ga_sort_desc()` -- a population of 0 or 1 is returned unchanged", {
  expect_equal(.ga_sort_desc(list()), list())
  one <- list(gi(c(3), "only"))
  expect_equal(.ga_sort_desc(one), one)
})

test_that("`.ga_sort_desc()` -- an all-empty-score population is returned unchanged", {
  pop <- list(gi(numeric(0), "a"), gi(numeric(0), "b"))
  expect_equal(.ga_sort_desc(pop), pop)
})

test_that("`.ga_sort_desc()` -- orders best-first, matching .ga_compare (length as final tie-break)", {
  pop <- list(
    gi(c(5, 1), "p"),
    gi(c(5), "q"),
    gi(c(6), "r"),
    gi(numeric(0), "s")
  )
  ord <- vapply(.ga_sort_desc(pop), function(x) x$rep, character(1))
  expect_equal(ord, c("r", "p", "q", "s"))
})


test_that("`.ga_pick()` -- selects the bucket the scaled uniform draw lands in", {
  weights <- c(3, 2, 1)
  total <- 6

  pick_with <- function(u) {
    local_mocked_bindings(runif = function(n, ...) u, .package = "stats")
    .ga_pick(weights, total)
  }

  expect_equal(pick_with(0), 1L) # rand 0   -> bucket 1
  expect_equal(pick_with(0.6), 2L) # rand 3.6 -> bucket 2
  expect_equal(pick_with(0.99), 3L) # rand 5.94-> bucket 3
})

test_that("`.ga_pick()` -- falls back to the last index when nothing is selected", {
  local_mocked_bindings(runif = function(n, ...) 1.5, .package = "stats") # rand > total
  expect_equal(.ga_pick(c(3, 2, 1), 6), 3L)
})


# Helper: for crossover without mutation, each position's (c1, c2) pair must be
# a permutation of (p1, p2) -- true for UNIFORM, SINGLEPOINT and MULTIPOINT.
expect_position_preserving <- function(p1, p2, c1, c2) {
  expect_equal(length(c1), length(p1))
  expect_equal(length(c2), length(p1))
  for (i in seq_along(p1)) {
    expect_setequal(c(c1[i], c2[i]), c(p1[i], p2[i]))
  }
}

test_that("`.ga_crossover_mutation()` -- no crossover yields empty children when the rate gate fails", {
  p1 <- c(TRUE, FALSE, TRUE, TRUE)
  p2 <- c(FALSE, FALSE, TRUE, FALSE)
  params <- list(ga_crossover_rate = 0, ga_mutation_rate = 0)

  kids <- .ga_crossover_mutation(p1, p2, params, "UNIFORM")
  expect_equal(kids[[1]], logical(0))
  expect_equal(kids[[2]], logical(0))
})

test_that("`.ga_crossover_mutation()` -- UNIFORM crossover preserves each position's parent pair", {
  set.seed(1)
  p1 <- c(TRUE, TRUE, FALSE, FALSE, TRUE)
  p2 <- c(FALSE, TRUE, TRUE, FALSE, FALSE)
  params <- list(ga_crossover_rate = 1, ga_mutation_rate = 0)

  kids <- .ga_crossover_mutation(p1, p2, params, "UNIFORM")
  expect_position_preserving(p1, p2, kids[[1]], kids[[2]])
})

test_that("`.ga_crossover_mutation()` -- SINGLEPOINT crossover preserves each position's parent pair", {
  set.seed(2)
  p1 <- c(TRUE, TRUE, FALSE, FALSE, TRUE, TRUE)
  p2 <- c(FALSE, TRUE, TRUE, FALSE, FALSE, TRUE)
  params <- list(ga_crossover_rate = 1, ga_mutation_rate = 0)

  kids <- .ga_crossover_mutation(p1, p2, params, "SINGLEPOINT")
  expect_position_preserving(p1, p2, kids[[1]], kids[[2]])
})

test_that("`.ga_crossover_mutation()` -- MULTIPOINT crossover preserves each position's parent pair", {
  set.seed(3)
  p1 <- as.logical(rep(c(TRUE, FALSE), length.out = 25))
  p2 <- as.logical(rep(c(FALSE, TRUE), length.out = 25))
  params <- list(ga_crossover_rate = 1, ga_mutation_rate = 0)

  kids <- .ga_crossover_mutation(p1, p2, params, "MULTIPOINT")
  expect_position_preserving(p1, p2, kids[[1]], kids[[2]])
})

test_that("`.ga_crossover_mutation()` -- full mutation flips every child bit", {
  set.seed(4)
  p1 <- c(TRUE, FALSE, TRUE, FALSE, TRUE)
  p2 <- c(FALSE, FALSE, TRUE, TRUE, TRUE)
  # crossover passes (rate 1); mutation flips every bit (rate 1)
  params <- list(ga_crossover_rate = 1, ga_mutation_rate = 1)

  kids <- .ga_crossover_mutation(p1, p2, params, "UNIFORM")
  # after a guaranteed flip, the per-position pair is the complement of (p1, p2)
  for (i in seq_along(p1)) {
    expect_setequal(c(kids[[1]][i], kids[[2]][i]), c(!p1[i], !p2[i]))
  }
})


test_that("`.ga_make_individual()` -- stores the genome and the component-scorer output", {
  fake_scores <- c(3, 1)
  local_mocked_bindings(
    component_scores_sorted = function(rep_logical, csr_offsets, csr_nbrs,
                                       z_vec, means, stds) {
      fake_scores
    }
  )
  genome <- c(TRUE, FALSE, TRUE)
  ind <- .ga_make_individual(genome, 0L, 0L, c(1, 2, 3), 0, 1)
  expect_equal(ind$rep, genome)
  expect_equal(ind$scores, fake_scores)
})


# Real (small) network for the final .find_subnetworks() materialisation:
# A-B-C triangle plus an isolated, strongly negative node D.
make_network <- function() {
  g <- igraph::graph_from_literal(A - B, B - C, C - A, D)
  list(
    nodes       = c("A", "B", "C", "D"),
    g           = g,
    csr_offsets = integer(5), # only forwarded to the (stubbed) scorer
    csr_nbrs    = integer(0)
  )
}

make_sc <- function() {
  list(
    z     = c(A = 2, B = 2, C = 2, D = -5),
    means = rep(0, 4),
    stds  = rep(1, 4),
    nodes = c("A", "B", "C", "D")
  )
}

make_params <- function(...) {
  defaults <- list(
    ga_population_size       = 10,
    ga_iterations            = 5,
    gene_init_prob           = 0.5,
    ga_crossover_rate        = 1,
    ga_mutation_rate         = 0,
    start_with_all_positives = TRUE,
    seed                     = 1234
  )
  utils::modifyList(defaults, list(...))
}

# Deterministic stand-in for the C++ scorer: fitness = sum of on-node z-scores.
# Maximised uniquely by the all-positive set {A, B, C}.
sum_scorer <- function(rep_logical, csr_offsets, csr_nbrs, z_vec, means, stds) {
  on <- which(rep_logical)
  if (length(on) == 0L) {
    return(numeric(0))
  }
  sum(z_vec[on])
}

test_that("`.genetic_algorithm()` -- an empty network returns an empty list without scoring", {
  with_mocked_bindings(
    {
      expect_equal(
        .genetic_algorithm(list(nodes = character(0)), make_sc(), make_params()),
        list()
      )
    },
    component_scores_sorted = function(...) stop("should not be called")
  )
})

test_that("`.genetic_algorithm()` -- start_with_all_positives seeds and retains the optimal module", {
  with_mocked_bindings(
    {
      res <- .genetic_algorithm(make_network(), make_sc(), make_params())
      expect_true(length(res) >= 1L)
      expect_setequal(res[[1]]$nodes, c("A", "B", "C")) # D excluded (negative z)
      expect_true(res[[1]]$score > 0)
    },
    component_scores_sorted = sum_scorer
  )
})

test_that("`.genetic_algorithm()` -- is reproducible for a fixed seed", {
  with_mocked_bindings(
    {
      r1 <- .genetic_algorithm(make_network(), make_sc(), make_params())
      r2 <- .genetic_algorithm(make_network(), make_sc(), make_params())
      nodes_of <- function(r) lapply(r, function(s) sort(s$nodes))
      expect_equal(nodes_of(r1), nodes_of(r2))
      expect_equal(
        vapply(r1, function(s) s$score, numeric(1)),
        vapply(r2, function(s) s$score, numeric(1))
      )
    },
    component_scores_sorted = sum_scorer
  )
})

test_that("`.genetic_algorithm()` -- a stochastic (no-elite-seed) run still returns well-formed subnetworks", {
  with_mocked_bindings(
    {
      res <- .genetic_algorithm(
        make_network(), make_sc(),
        make_params(start_with_all_positives = FALSE)
      )
      expect_true(is.list(res))
      for (s in res) {
        expect_true(is.character(s$nodes))
        expect_true(is.numeric(s$score))
      }
    },
    component_scores_sorted = sum_scorer
  )
})

test_that("`.genetic_algorithm()` -- stops early after 50 unchanged generations", {
  with_mocked_bindings(
    {
      msgs <- testthat::capture_messages(
        .genetic_algorithm(
          make_network(), make_sc(),
          make_params(ga_iterations = 1000), # far more than convergence needs
          verbose = TRUE
        )
      )
      expect_true(any(grepl("did not improve in 50 generations", msgs)))

      gens <- as.integer(sub(
        ".*GA generation ([0-9]+).*", "\\1",
        grep("GA generation", msgs, value = TRUE)
      ))
      expect_true(max(gens) < 1000) # terminated on the no-improvement rule, not ga_iterations
    },
    component_scores_sorted = sum_scorer
  )
})

test_that("`.genetic_algorithm()` -- emits per-generation messages only when verbose", {
  with_mocked_bindings(
    {
      expect_message(
        .genetic_algorithm(make_network(), make_sc(), make_params(), verbose = TRUE),
        "GA generation"
      )
      expect_no_message(
        .genetic_algorithm(make_network(), make_sc(), make_params(), verbose = FALSE)
      )
    },
    component_scores_sorted = sum_scorer
  )
})

test_that("`.genetic_algorithm()` -- truncates the bred population when an odd size overshoots", {
  with_mocked_bindings(
    {
      # Children are bred two at a time, so an odd ga_population_size makes the
      # inner breeding loop overshoot pop_size by one. That triggers
      # `new_pop <- new_pop[seq_len(pop_size)]`. The run must still behave: with
      # the all-positive elite seeded, the optimal module {A, B, C} is retained.
      res <- .genetic_algorithm(
        make_network(), make_sc(),
        make_params(ga_population_size = 11) # odd -> overshoot by one
      )
      expect_true(length(res) >= 1L)
      expect_setequal(res[[1]]$nodes, c("A", "B", "C"))
      expect_true(res[[1]]$score > 0)
    },
    component_scores_sorted = sum_scorer
  )
})
