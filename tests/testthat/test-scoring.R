## Tests for agglomerated term scoring functions - Jan 2024

test_that("`score_terms()` -- returns score matrix", {
    mockery::stub(score_terms, "graphics::plot", NULL)

    small_result <- example_pathfindR_output[1:3, ]
    expect_is(score_terms(enrichment_table = small_result, exp_mat = example_experiment_matrix,
        plot_hmap = FALSE), "matrix")
    expect_is(score_terms(enrichment_table = small_result, exp_mat = example_experiment_matrix,
        plot_hmap = TRUE), "matrix")
    expect_is(score_terms(enrichment_table = small_result, exp_mat = example_experiment_matrix,
        cases = colnames(example_experiment_matrix)[1:3], plot_hmap = TRUE), "matrix")
})

test_that("`score_terms()` -- matches gene symbols correctly", {
  toy_result <- data.frame(
    ID = c("gset1", "gset2"),
    Term_Description = c("gset1", "gset2"),
    Up_regulated = "",
    Down_regulated = c(
      paste(paste0("Gene_", c(1, 3, 5)), collapse = ", "),
      paste(paste0("Gene_", c(6, 8)), collapse = ", ")
    )
  )
  toy_result2 <- data.frame(
    ID = c("gset1", "gset2"),
    Term_Description = c("gset1", "gset2"),
    Up_regulated = "",
    Down_regulated = c(
      paste(paste0("Dummy_", c(1, 3, 5)), collapse = ", "),
      paste(paste0("Gene_", c(6, 8)), collapse = ", ")
    )
  )
  toy_exp_mat <- matrix(
    rnorm(40), nrow = 10, ncol = 4, dimnames = list(paste0("gene_", 1:10), paste0("subject_", 1:4))
  )
  expect_is(res_mat <- score_terms(enrichment_table = toy_result, exp_mat = toy_exp_mat,
                        plot_hmap = FALSE), "matrix")
  expect_equal(nrow(res_mat), 2)


  expect_is(res_mat <- score_terms(enrichment_table = toy_result2, exp_mat = toy_exp_mat,
                                   plot_hmap = FALSE), "matrix")
  expect_equal(nrow(res_mat), 1)
})

test_that("`score_terms()` -- argument checks work", {
    expect_error(score_terms(enrichment_table = example_pathfindR_output, exp_mat = example_experiment_matrix,
        use_description = "INVALID"), "`use_description` should either be TRUE or FALSE")

    expect_error(score_terms(enrichment_table = example_pathfindR_output, exp_mat = example_experiment_matrix,
        plot_hmap = "INVALID"), "`plot_hmap` should either be TRUE or FALSE")

    expect_error(score_terms(enrichment_table = list(), exp_mat = example_experiment_matrix),
        "`enrichment_table` should be a data frame of enrichment results")

    tmp <- example_pathfindR_output[, -c(1, 2)]
    nec_cols <- c("ID", "Up_regulated", "Down_regulated")
    expect_error(score_terms(enrichment_table = tmp, exp_mat = example_experiment_matrix),
        paste0("`enrichment_table` should contain all of ", paste(dQuote(nec_cols),
            collapse = ", ")))
    nec_cols <- c("Term_Description", "Up_regulated", "Down_regulated")
    expect_error(score_terms(enrichment_table = tmp, exp_mat = example_experiment_matrix,
        use_description = TRUE), paste0("`enrichment_table` should contain all of ",
        paste(dQuote(nec_cols), collapse = ", ")))

    expect_error(score_terms(enrichment_table = example_pathfindR_output, exp_mat = list()),
        "`exp_mat` should be a matrix")

    expect_error(score_terms(enrichment_table = example_pathfindR_output, exp_mat = example_experiment_matrix,
        cases = list()), "`cases` should be a vector")
    expect_error(score_terms(enrichment_table = example_pathfindR_output, exp_mat = example_experiment_matrix,
        cases = LETTERS), "Missing `cases` in `exp_mat`")
})

test_that("duplicated term descriptions test", {
    small_result <- example_pathfindR_output[1:2, ]
    small_result$Term_Description <- small_result$Term_Description[1]
    expect_is(score_terms(enrichment_table = small_result, exp_mat = example_experiment_matrix,
        use_description = TRUE, plot_hmap = FALSE), "matrix")
})

test_that("`plot_scores()` -- creates term score heatmap ggplot object with correct labels",
    {
        score_mat <- score_terms(example_pathfindR_output[1:3, ], example_experiment_matrix,
            plot_hmap = FALSE)

        # default
        g <- plot_scores(score_mat)
        expect_is(g, "ggplot")
        
        labels <- ggplot2::get_labs(g)
        expect_identical(labels$fill, "Score")
        expect_identical(labels$x, "Sample")
        expect_identical(labels$y, "Term")

        # cases provided
        g <- plot_scores(score_mat, cases = colnames(score_mat)[1:3])
        expect_is(g, "ggplot")
        
        labels <- ggplot2::get_labs(g)
        expect_identical(labels$fill, "Score")
        expect_identical(labels$x, "Sample")
        expect_identical(labels$y, "Term")

        # default - label_samples = FALSE
        g <- plot_scores(score_mat, label_samples = FALSE)
        expect_is(g, "ggplot")
        
        labels <- ggplot2::get_labs(g)
        expect_identical(labels$fill, "Score")
        expect_identical(labels$x, "Sample")
        expect_identical(labels$y, "Term")

        # cases provided - label_samples = FALSE
        g <- plot_scores(score_mat, cases = colnames(score_mat)[1:3], label_samples = FALSE)
        expect_is(g, "ggplot")
        
        labels <- ggplot2::get_labs(g)
        expect_identical(labels$fill, "Score")
        expect_identical(labels$x, "Sample")
        expect_identical(labels$y, "Term")
    })

test_that("`plot_scores()` -- argument checks work", {
    expect_error(plot_scores(score_matrix = c()), "`score_matrix` should be a matrix")
    expect_error(plot_scores(score_matrix = data.frame()), "`score_matrix` should be a matrix")
    expect_error(plot_scores(score_matrix = list()), "`score_matrix` should be a matrix")

    mat <- matrix(1, nrow = 3, ncol = 2, dimnames = list(paste0("T", 1:3), c("A",
        "B")))

    expect_error(plot_scores(score_matrix = mat, cases = list()), "`cases` should be a vector")
    expect_error(plot_scores(score_matrix = mat, cases = c("A", "B", "C")), "Missing `cases` in `score_matrix`")

    expect_error(plot_scores(score_matrix = mat, label_samples = "INVALID"), "`label_samples` should be TRUE or FALSE")

    expect_error(plot_scores(score_matrix = mat, case_title = 1), "`case_title` should be a single character value")
    expect_error(plot_scores(score_matrix = mat, case_title = rep("z", 3)), "`case_title` should be a single character value")

    expect_error(plot_scores(score_matrix = mat, control_title = 1), "`control_title` should be a single character value")
    expect_error(plot_scores(score_matrix = mat, control_title = rep("z", 3)), "`control_title` should be a single character value")

    expect_error(plot_scores(score_matrix = mat, low = ""))
    expect_error(plot_scores(score_matrix = mat, mid = ""))
    expect_error(plot_scores(score_matrix = mat, high = ""))
})
