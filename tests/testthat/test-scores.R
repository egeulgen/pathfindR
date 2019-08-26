test_that("Score matrix is returned", {
  expect_is(calculate_pw_scores(RA_output,
                                RA_exp_mat,
                                plot_hmap = FALSE),
            "matrix")
  expect_is(calculate_pw_scores(RA_output,
                                RA_exp_mat,
                                plot_hmap = TRUE),
            "matrix")
  expect_is(calculate_pw_scores(RA_output,
                                RA_exp_mat,
                                cases = colnames(RA_exp_mat)[1:3],
                                plot_hmap = TRUE),
            "matrix")
})

test_that("If cases do not match the ones in exp_mat,
          error is returned during score calculation", {
            expect_error(calculate_pw_scores(RA_output,
                                             RA_exp_mat,
                                             cases = c("A", "B", "C"),
                                             plot_hmap = FALSE),
                         "Missing cases in the expression matrix!")
          })

test_that("Heatmap ggplot object is created with correct labels", {
  score_mat <- calculate_pw_scores(RA_output, RA_exp_mat, plot_hmap = FALSE)

  # default
  g <- plot_scores(score_mat)
  expect_is(g, "ggplot")
  expect_identical(g$labels$fill, "Pathway\nscore")
  expect_identical(g$labels$x, "Sample")
  expect_identical(g$labels$y, "Pathway")

  # cases provided
  g <- plot_scores(score_mat, cases = colnames(score_mat)[1:3])
  expect_is(g, "ggplot")
  expect_identical(g$labels$fill, "Pathway\nscore")
  expect_identical(g$labels$x, "Sample")
  expect_identical(g$labels$y, "Pathway")

})


test_that("Heatmap creation errors", {
  score_mat <- calculate_pw_scores(RA_output, RA_exp_mat, plot_hmap = FALSE)

  expect_error(plot_scores(score_mat, case_control_titles = c("A", "B", "C")),
               "\"case_control_titles\" must contain two elements!")
  expect_error(plot_scores(score_mat,
                           cases = c("A", "B", "C")),
               "Missing cases in the score matrix!")
})



