##################################################
## Project: pathfindR
## Script purpose: Testthat testing script for
## active snw search functions
## Date: Aug 30, 2019
## Author: Ege Ulgen
##################################################

# filterActiveSnws --------------------------------------------------------
sample_path <- system.file("extdata/resultActiveSubnetworkSearch.txt",
                           package = "pathfindR")

test_that("Filter function returns list object", {

  tmp_filtered <- filterActiveSnws(sample_path, RA_input$Gene.symbol)
  expect_is(tmp_filtered, "list")
  expect_is(tmp_filtered[[1]], "character")

  # empty file case
  expect_equal(suppressWarnings(filterActiveSnws("", RA_input$Gene.symbol)),
               NULL)
})

test_that("Filter function arguments work OK", {

  # supplied num. of genes
  tmp_filtered1 <- filterActiveSnws(sample_path, RA_input$Gene.symbol)
  tmp_filtered2 <- filterActiveSnws(sample_path, RA_input$Gene.symbol[1:100])
  expect_equal(length(tmp_filtered1) > length(tmp_filtered2), TRUE)


  # score_quan_thr
  tmp_filtered <- filterActiveSnws(active_snw_path = sample_path,
                                    signif_genes = RA_input$Gene.symbol,
                                    score_quan_thr = 0.8) # default
  expect_equal(length(tmp_filtered1) <= 200 , TRUE)

  tmp_filtered <- filterActiveSnws(active_snw_path = sample_path,
                                    signif_genes = RA_input$Gene.symbol,
                                    score_quan_thr = 0.9)
  expect_equal(length(tmp_filtered) <= 100 , TRUE)

  # sig_gene_thr
  tmp_filtered <- filterActiveSnws(active_snw_path = sample_path,
                                   signif_genes = RA_input$Gene.symbol,
                                   sig_gene_thr = 0)
  expect_equal(length(tmp_filtered) == 200 , TRUE)

  tmp_filtered1 <- filterActiveSnws(active_snw_path = sample_path,
                                   signif_genes = RA_input$Gene.symbol,
                                   sig_gene_thr = 10) # default

  tmp_filtered2 <- filterActiveSnws(active_snw_path = sample_path,
                                    signif_genes = RA_input$Gene.symbol,
                                    sig_gene_thr = 20)

  expect_equal(length(tmp_filtered1) <= 200  , TRUE)
  expect_equal(length(tmp_filtered1) > length(tmp_filtered2)  , TRUE)
})




