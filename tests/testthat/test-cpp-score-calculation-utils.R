test_that("get_java_zscores() -- p <= 0 maps to +Inf", {
  z <- get_java_zscores(c(0, -0.5))
  expect_true(all(is.infinite(z) & z > 0))
})

test_that("get_java_zscores() -- p >= 1 maps to -Inf", {
  z <- get_java_zscores(c(1, 1.5))
  expect_true(all(is.infinite(z) & z < 0))
})

test_that("get_java_zscores() -- p = 0.5 maps to ~0", {
  expect_equal(get_java_zscores(0.5), 0, tolerance = 1e-4)
})

test_that("get_java_zscores() -- the two tails are antisymmetric", {
  # p in (0,0.5] uses oneMinus(p); p in (0.5,1) uses -oneMinus(1-p)
  z <- get_java_zscores(c(0.2, 0.8))
  expect_equal(z[2], -z[1])
})

test_that("get_java_zscores() -- approximates standard-normal upper quantiles and is monotone", {
  z <- get_java_zscores(c(0.05, 0.2, 0.5))
  expect_equal(z[1], qnorm(0.95), tolerance = 1e-2) # ~1.6449
  expect_true(z[1] > z[2] && z[2] > z[3]) # decreasing in p
})


test_that("get_java_mc_calibration() -- with identical z-scores the shuffle is irrelevant, giving exact means/stds", {
  # Every permutation of a constant vector has the same prefix sums, so for size
  # k>1 every trial scores k/sqrt(k) = sqrt(k): mean = sqrt(k), variance = 0.
  mc <- get_java_mc_calibration(rep(1, 4), trials = 2000, seed = 1234)

  expect_length(mc$means, 4L)
  expect_length(mc$stds, 4L)

  expect_equal(mc$means, c(0, sqrt(2), sqrt(3), sqrt(4)), tolerance = 1e-9)
  # size 1 is hard-coded to mean 0 / std 1; sizes >1 get sqrt(0 + 1e-7)
  expect_equal(mc$stds, c(1, sqrt(1e-7), sqrt(1e-7), sqrt(1e-7)), tolerance = 1e-6)
})

test_that("get_java_mc_calibration() -- size-1 entry is fixed at mean 0 / std 1", {
  mc <- get_java_mc_calibration(c(0.5, 1, 1.5, 2), trials = 100, seed = 7)
  expect_equal(mc$means[1], 0)
  expect_equal(mc$stds[1], 1)
})

test_that("get_java_mc_calibration() -- is reproducible for a fixed seed", {
  z <- c(0.5, 1.0, 1.5, 2.0, -0.5)
  a <- get_java_mc_calibration(z, trials = 500, seed = 99)
  b <- get_java_mc_calibration(z, trials = 500, seed = 99)
  expect_equal(a$means, b$means)
  expect_equal(a$stds, b$stds)
})
