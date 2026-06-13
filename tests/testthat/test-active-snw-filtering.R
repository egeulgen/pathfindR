# set up input data
input_data_frame <- example_pathfindR_input[1:10, c(1, 3)]
colnames(input_data_frame) <- c("GENE", "P_VALUE")
example_snws_len <- 1000
example_snw_output <- system.file("extdata", "resultActiveSubnetworkSearch.txt",
  package = "pathfindR"
)

test_that("`filterActiveSnws()` -- returns expected list object", {
  snws_filtered <- filterActiveSnws(active_snw_path = example_snw_output, sig_genes_vec = input_data_frame$GENE)
  expect_is(snws_filtered, "list")
  expect_length(snws_filtered, 2)
  expect_is(snws_filtered$subnetworks, "list")
  expect_is(snws_filtered$scores, "numeric")

  expect_is(snws_filtered$subnetworks[[1]], "character")
  expect_true(length(snws_filtered$subnetworks) <= example_snws_len)

  # empty file case
  empty_path <- tempfile("empty", fileext = ".txt")
  file.create(empty_path)
  expect_null(suppressWarnings(filterActiveSnws(active_snw_path = empty_path, sig_genes_vec = input_data_frame$GENE)))
})

test_that("`filterActiveSnws()` -- `score_quan_thr` works", {
  snws_filtered <- filterActiveSnws(
    active_snw_path = example_snw_output, sig_genes_vec = example_pathfindR_input$Gene.symbol,
    score_quan_thr = -1, sig_gene_thr = 0
  )
  expect_length(snws_filtered$subnetworks, example_snws_len)

  for (q_thr in seq(0.1, 1, by = 0.1)) {
    snws_filtered <- filterActiveSnws(
      active_snw_path = example_snw_output, sig_genes_vec = example_pathfindR_input$Gene.symbol,
      score_quan_thr = q_thr, sig_gene_thr = 0
    )
    exp_len <- example_snws_len * (1 - q_thr)
    expect_length(snws_filtered$subnetworks, as.integer(exp_len + 0.5))
  }
})

test_that("`filterActiveSnws()` -- `sig_gene_thr` works", {
  snws_filtered1 <- filterActiveSnws(
    active_snw_path = example_snw_output, sig_genes_vec = example_pathfindR_input$Gene.symbol,
    sig_gene_thr = 0.02, score_quan_thr = -1
  )
  snws_filtered2 <- filterActiveSnws(
    active_snw_path = example_snw_output, sig_genes_vec = example_pathfindR_input$Gene.symbol,
    sig_gene_thr = 0.1, score_quan_thr = -1
  )

  expect_true(length(snws_filtered2$subnetworks) < example_snws_len)
  expect_true(length(snws_filtered1$subnetworks) > length(snws_filtered2$subnetworks))
})

test_that("`filterActiveSnws()` -- argument checks work", {
  expect_error(
    filterActiveSnws(active_snw_path = "this/is/not/a/valid/path"),
    "The active subnetwork file does not exist! Check the `active_snw_path` argument"
  )

  expect_error(
    filterActiveSnws(active_snw_path = example_snw_output, sig_genes_vec = list()),
    "`sig_genes_vec` should be a vector"
  )

  expect_error(filterActiveSnws(
    active_snw_path = example_snw_output, sig_genes_vec = example_pathfindR_input$Gene.symbol,
    score_quan_thr = "INVALID"
  ), "`score_quan_thr` should be numeric")
  expect_error(filterActiveSnws(
    active_snw_path = example_snw_output, sig_genes_vec = example_pathfindR_input$Gene.symbol,
    score_quan_thr = -2
  ), "`score_quan_thr` should be in \\[0, 1\\] or -1 \\(if not filtering\\)")
  expect_error(filterActiveSnws(
    active_snw_path = example_snw_output, sig_genes_vec = example_pathfindR_input$Gene.symbol,
    score_quan_thr = 2
  ), "`score_quan_thr` should be in \\[0, 1\\] or -1 \\(if not filtering\\)")

  expect_error(filterActiveSnws(
    active_snw_path = example_snw_output, sig_genes_vec = example_pathfindR_input$Gene.symbol,
    sig_gene_thr = "INVALID"
  ), "`sig_gene_thr` should be numeric")
  expect_error(filterActiveSnws(
    active_snw_path = example_snw_output, sig_genes_vec = example_pathfindR_input$Gene.symbol,
    sig_gene_thr = -1
  ), "`sig_gene_thr` should be in \\[0, 1\\]")
})
