# Triangle 0-1-2, neighbour ids are 1-based.
tri_nbr <- list(c(2L, 3L), c(1L, 3L), c(1L, 2L))
tri_names <- c("A", "B", "C")

run_tri <- function(z, max_depth = 1L, search_depth = 1L,
                    means = rep(0, 3), stds = rep(1, 3),
                    overlap = 0.5, subnetwork_num = 1000L) {
  run_greedy_search(
    nbr_idx = tri_nbr, z_vec = z, sc_means = means, sc_stds = stds,
    node_names = tri_names, max_depth = max_depth, search_depth = search_depth,
    n_nodes = 3L, overlap_threshold = overlap, subnetwork_num = subnetwork_num
  )
}

test_that("run_greedy_search() -- grows a positive triangle into one module (max_depth>0 uses the within-mask)", {
  res <- run_tri(z = c(2, 2, 2), max_depth = 1L, search_depth = 1L)
  expect_length(res, 1L)
  expect_setequal(res[[1]]$nodes, c("A", "B", "C"))
  expect_equal(res[[1]]$score, 6 / sqrt(3))
})

test_that("run_greedy_search() -- max_depth = 0 disables the within-mask but still finds the module", {
  res <- run_tri(z = c(2, 2, 2), max_depth = 0L, search_depth = 2L)
  expect_length(res, 1L)
  expect_setequal(res[[1]]$nodes, c("A", "B", "C"))
  expect_equal(res[[1]]$score, 6 / sqrt(3))
})

test_that("run_greedy_search() -- an all-negative network yields no positive-scoring subnetwork", {
  # every candidate component scores <= 0, so nothing is recorded
  res <- run_tri(z = c(-1, -1, -1))
  expect_length(res, 0L)
})

test_that("run_greedy_search() -- overlap filtering keeps distinct modules and drops duplicates", {
  # Two disjoint triangles {A,B,C} and {D,E,F}. Every seed rediscovers its own
  # triangle, producing duplicate candidates that the overlap filter (ratio 1.0)
  # collapses, leaving exactly the two non-overlapping modules.
  nbr <- list(
    c(2L, 3L), c(1L, 3L), c(1L, 2L), # triangle 1: A B C
    c(5L, 6L), c(4L, 6L), c(4L, 5L) # triangle 2: D E F
  )
  names6 <- c("A", "B", "C", "D", "E", "F")
  res <- run_greedy_search(
    nbr_idx = nbr, z_vec = rep(2, 6), sc_means = rep(0, 6), sc_stds = rep(1, 6),
    node_names = names6, max_depth = 1L, search_depth = 1L, n_nodes = 6L,
    overlap_threshold = 0.5, subnetwork_num = 1000L
  )
  expect_length(res, 2L)
  got <- lapply(res, function(s) sort(s$nodes))
  expect_true(any(vapply(got, function(g) setequal(g, c("A", "B", "C")), logical(1))))
  expect_true(any(vapply(got, function(g) setequal(g, c("D", "E", "F")), logical(1))))
  for (s in res) expect_equal(s$score, 6 / sqrt(3))
})

test_that("run_greedy_search() -- is deterministic", {
  expect_equal(run_tri(z = c(2, 2, 2)), run_tri(z = c(2, 2, 2)))
})
