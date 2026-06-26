# A minimal subnetwork object.
sn <- function(nodes, score) list(nodes = nodes, score = score)

# Dummy inputs; contents are irrelevant because the strategies are mocked.
fake_network <- list(nodes = c("A", "B", "C"))
fake_sc <- list()
fake_params <- list()

test_that("`active_subnetwork_search()` -- an empty network returns an empty list without running any search", {
  with_mocked_bindings(
    {
      result <- active_subnetwork_search(
        list(nodes = character(0)), fake_sc,
        method = "GR", params = fake_params
      )
      expect_equal(result, list())
    },
    .greedy_search = function(...) stop("should not be called"),
    .simulated_annealing = function(...) stop("should not be called"),
    .genetic_algorithm = function(...) stop("should not be called")
  )
})

test_that("`active_subnetwork_search()` -- an invalid method is rejected by match.arg", {
  expect_error(
    active_subnetwork_search(fake_network, fake_sc, method = "XX", params = fake_params),
    "should be one of"
  )
})

test_that("`active_subnetwork_search()` -- dispatches to the strategy named by `method`", {
  with_mocked_bindings(
    {
      gr <- active_subnetwork_search(fake_network, fake_sc, method = "GR", params = fake_params)
      sa <- active_subnetwork_search(fake_network, fake_sc, method = "SA", params = fake_params)
      ga <- active_subnetwork_search(fake_network, fake_sc, method = "GA", params = fake_params)

      expect_equal(gr[[1]]$nodes, "from_GR")
      expect_equal(sa[[1]]$nodes, "from_SA")
      expect_equal(ga[[1]]$nodes, "from_GA")
    },
    .greedy_search = function(...) list(sn("from_GR", 1)),
    .simulated_annealing = function(...) list(sn("from_SA", 1)),
    .genetic_algorithm = function(...) list(sn("from_GA", 1))
  )
})

test_that("`active_subnetwork_search()` -- defaults to the greedy (GR) strategy", {
  with_mocked_bindings(
    {
      result <- active_subnetwork_search(fake_network, fake_sc, params = fake_params)
      expect_equal(result[[1]]$nodes, "from_GR")
    },
    .greedy_search = function(...) list(sn("from_GR", 1)),
    .simulated_annealing = function(...) stop("should not be called"),
    .genetic_algorithm = function(...) stop("should not be called")
  )
})

test_that("`active_subnetwork_search()` -- results are sorted by score descending", {
  with_mocked_bindings(
    {
      result <- active_subnetwork_search(fake_network, fake_sc, method = "GR", params = fake_params)
      scores <- vapply(result, function(s) s$score, numeric(1))
      expect_equal(scores, c(5, 3, 1))
      expect_equal(result[[1]]$nodes, "b")
    },
    .greedy_search = function(...) list(sn("a", 1), sn("b", 5), sn("c", 3))
  )
})

test_that("`active_subnetwork_search()` -- non-positive-scoring subnetworks are dropped", {
  with_mocked_bindings(
    {
      result <- active_subnetwork_search(fake_network, fake_sc, method = "GR", params = fake_params)
      scores <- vapply(result, function(s) s$score, numeric(1))
      expect_equal(scores, c(3, 1)) # 0 and -2 removed, remainder sorted desc
      expect_true(all(scores > 0))
    },
    .greedy_search = function(...) list(sn("a", -2), sn("b", 3), sn("c", 0), sn("d", 1))
  )
})

test_that("`active_subnetwork_search()` -- returns an empty list when no subnetwork scores above zero", {
  with_mocked_bindings(
    {
      result <- active_subnetwork_search(fake_network, fake_sc, method = "GR", params = fake_params)
      expect_equal(result, list())
    },
    .greedy_search = function(...) list(sn("a", 0), sn("b", -1))
  )
})

test_that("`active_subnetwork_search()` -- returns an empty list when the strategy finds nothing", {
  with_mocked_bindings(
    {
      result <- active_subnetwork_search(fake_network, fake_sc, method = "GR", params = fake_params)
      expect_equal(result, list())
    },
    .greedy_search = function(...) list()
  )
})

test_that("`active_subnetwork_search()` -- emits progress messages only when verbose", {
  with_mocked_bindings(
    {
      expect_message(
        active_subnetwork_search(fake_network, fake_sc,
          method = "GR",
          params = fake_params, verbose = TRUE
        ),
        "Searching subnetworks"
      )
      expect_no_message(
        active_subnetwork_search(fake_network, fake_sc,
          method = "GR",
          params = fake_params, verbose = FALSE
        )
      )
    },
    .greedy_search = function(...) list(sn("a", 1))
  )
})

test_that("`active_subnetwork_search()` -- forwards network, score_context, params and verbose to the strategy", {
  captured <- new.env(parent = emptyenv())
  with_mocked_bindings(
    {
      active_subnetwork_search(fake_network, fake_sc,
        method = "SA",
        params = fake_params, verbose = TRUE
      )
    },
    .simulated_annealing = function(network, score_context, params, verbose) {
      captured$network <- network
      captured$score_context <- score_context
      captured$params <- params
      captured$verbose <- verbose
      list(sn("x", 1))
    }
  )
  expect_identical(captured$network, fake_network)
  expect_identical(captured$score_context, fake_sc)
  expect_identical(captured$params, fake_params)
  expect_true(captured$verbose)
})
