##################################################
## Project: pathfindR
## Script purpose: Testthat testing script for
## agglomerated term scoring functions
## Date: May 24, 2019
## Author: Ege Ulgen
##################################################

# score_terms --------------------------------------------------------
test_that("`score_terms()` returns score matrix", {
  tmp_res <- RA_output[1:3, ]
  expect_is(score_terms(enrichment_table = tmp_res,
                        exp_mat = RA_exp_mat,
                        plot_hmap = FALSE),
            "matrix")
  expect_is(score_terms(enrichment_table = tmp_res,
                        exp_mat = RA_exp_mat,
                        plot_hmap = TRUE),
            "matrix")
  expect_is(score_terms(enrichment_table = tmp_res,
                        exp_mat = RA_exp_mat,
                        cases = colnames(RA_exp_mat)[1:3],
                        plot_hmap = TRUE),
            "matrix")
})

test_that("`score_terms()` arg checks work", {
  expect_error(score_terms(enrichment_table = RA_output,
                           exp_mat = RA_exp_mat,
                           use_description = "INVALID"),
               "`use_description` should either be TRUE or FALSE")

  expect_error(score_terms(enrichment_table = RA_output,
                           exp_mat = RA_exp_mat,
                           plot_hmap = "INVALID"),
               "`plot_hmap` should either be TRUE or FALSE")

  expect_error(score_terms(enrichment_table = list(),
                           exp_mat = RA_exp_mat),
               "`enrichment_table` should be a data frame of enrichment results")

  tmp <- RA_output[, -c(1, 2)]
  nec_cols <- c("ID", "Up_regulated", "Down_regulated")
  expect_error(score_terms(enrichment_table = tmp,
                           exp_mat = RA_exp_mat),
               paste0("`enrichment_table` should contain all of ",
                      paste(dQuote(nec_cols), collapse = ", ")))
  nec_cols <- c("Term_Description", "Up_regulated", "Down_regulated")
  expect_error(score_terms(enrichment_table = tmp,
                           exp_mat = RA_exp_mat,
                           use_description = TRUE),
               paste0("`enrichment_table` should contain all of ",
                      paste(dQuote(nec_cols), collapse = ", ")))

  expect_error(score_terms(enrichment_table = RA_output,
                           exp_mat = list()),
               "`exp_mat` should be a matrix")

  expect_error(score_terms(enrichment_table = RA_output,
                           exp_mat = RA_exp_mat,
                           cases = list()),
               "`cases` should be a vector")
  expect_error(score_terms(enrichment_table = RA_output,
                           exp_mat = RA_exp_mat,
                           cases = LETTERS),
               "Missing `cases` in `exp_mat`")
})

# plot_scores -------------------------------------------------------------
test_that("`plot_scores()` creates term score heatmap ggplot object with correct labels", {
  score_mat <- score_terms(RA_output[1:3, ], RA_exp_mat, plot_hmap = FALSE)

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

test_that("`plot_scores()` arg checks work", {
  expect_error(plot_scores(score_matrix = list()),
               "`score_matrix` should be a matrix")
  mat <- matrix(1, nrow = 3, ncol = 2,
                dimnames = list(paste0("T", 1:3),
                                c("A", "B")))

  expect_error(plot_scores(score_matrix = mat,
                           cases = list()),
               "`cases` should be a vector")
  expect_error(plot_scores(score_matrix = mat,
                           cases = c("A", "B", "C")),
               "Missing `cases` in `score_matrix`")

  expect_error(plot_scores(score_matrix = mat,
                           label_samples = "INVALID"),
               "`label_samples` should be TRUE or FALSE")

  expect_error(plot_scores(score_matrix = mat,
                           case_title = 1),
               "`case_title` should be a single character value")
  expect_error(plot_scores(score_matrix = mat,
                           case_title = rep("z", 3)),
               "`case_title` should be a single character value")

  expect_error(plot_scores(score_matrix = mat,
                           control_title = 1),
               "`control_title` should be a single character value")
  expect_error(plot_scores(score_matrix = mat,
                           control_title = rep("z", 3)),
               "`control_title` should be a single character value")
})
