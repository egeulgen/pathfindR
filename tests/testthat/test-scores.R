##################################################
## Project: pathfindR
## Script purpose: Testthat testing script for
## agglomerated term scoring functions
## Date: Oct 20, 2019
## Author: Ege Ulgen
##################################################

# calculate_scores --------------------------------------------------------
test_that("Score matrix is returned", {
  expect_is(calculate_scores(RA_output,
                             RA_exp_mat,
                             plot_hmap = FALSE),
            "matrix")
  expect_is(calculate_scores(RA_output,
                             RA_exp_mat,
                             plot_hmap = TRUE),
            "matrix")
  expect_is(calculate_scores(RA_output,
                             RA_exp_mat,
                             cases = colnames(RA_exp_mat)[1:3],
                             plot_hmap = TRUE),
            "matrix")
})

test_that("If cases do not match the ones in exp_mat,
          error is returned during score calculation", {
            expect_error(calculate_scores(RA_output,
                                          RA_exp_mat,
                                          cases = c("A", "B", "C"),
                                          plot_hmap = FALSE),
                         "Missing cases in the expression matrix!")
          })

test_that("Error thrown if class of `use_description` is not logical", {
  expect_error(calculate_scores(RA_output,
                                RA_exp_mat,
                                use_description = "WRONG"),
               "`use_description` must either be TRUE or FALSE!")
})

# plot_scores -------------------------------------------------------------
test_that("Term score heatmap ggplot object is created with correct labels", {
  score_mat <- calculate_scores(RA_output, RA_exp_mat, plot_hmap = FALSE)

  # default
  g <- plot_scores(score_mat)
  expect_is(g, "ggplot")
  expect_identical(g$labels$fill, "Score")
  expect_identical(g$labels$x, "Sample")
  expect_identical(g$labels$y, "Term")

  # cases provided
  g <- plot_scores(score_mat, cases = colnames(score_mat)[1:3])
  expect_is(g, "ggplot")
  expect_identical(g$labels$fill, "Score")
  expect_identical(g$labels$x, "Sample")
  expect_identical(g$labels$y, "Term")

  # default - label_samples = FALSE
  g <- plot_scores(score_mat, label_samples = FALSE)
  expect_is(g, "ggplot")
  expect_identical(g$labels$fill, "Score")
  expect_identical(g$labels$x, "Sample")
  expect_identical(g$labels$y, "Term")

  # cases provided - label_samples = FALSE
  g <- plot_scores(score_mat,
                   cases = colnames(score_mat)[1:3],
                   label_samples = FALSE)
  expect_is(g, "ggplot")
  expect_identical(g$labels$fill, "Score")
  expect_identical(g$labels$x, "Sample")
  expect_identical(g$labels$y, "Term")
})

test_that("Term score heatmap creation errors", {
  score_mat <- calculate_scores(RA_output, RA_exp_mat, plot_hmap = FALSE)

  expect_error(plot_scores(score_mat, case_control_titles = c("A", "B", "C")),
               "\"case_control_titles\" must contain exactly two elements!")
  expect_error(plot_scores(score_mat,
                           cases = c("A", "B", "C")),
               "Missing cases in the score matrix!")
})
