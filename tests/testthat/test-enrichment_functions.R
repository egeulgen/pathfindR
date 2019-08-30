##################################################
## Project: pathfindR
## Script purpose: Testthat testing script for
## enrichment-related functions
## Date: Aug 30, 2019
## Author: Ege Ulgen
##################################################

# hyperg_test -------------------------------------------------------------
test_that("hyperg_test returns a p value"{
  expect_is(tmp_p <- hyperg_test(letters[1:5], letters[2:5], letters),
            "numeric")
  expect_true(tmp_p >= 0 & tmp_p <= 1)

  expect_is(tmp_p <- hyperg_test(letters[1:5], letters[2:10], letters),
            "numeric")
  expect_true(tmp_p >= 0 & tmp_p <= 1)
})

test_that("hyperg_test sanity checks work"{
  expect_error(hyperg_test(pw_genes = c(letters, "xx"),
                           chosen_genes = letters[1:3],
                           all_genes = letters),
               "`pw_genes` cannot be larger than `all_genes`!")

  expect_error(hyperg_test(pw_genes = letters[1:5],
                           chosen_genes = c(letters, "xx"),
                           all_genes = letters),
               "`chosen_genes` cannot be larger than `all_genes`!")
})

# enrichment --------------------------------------------------------------
tmp_gset <- kegg_genes[1:50]
tmp_gset_names <- kegg_pathways[names(tmp_gset)]
tmp_gset_genes <- unlist(tmp_gset)

test_that("enrichment function returns a data frame"{

  # default
  expect_is(tmp1 <- enrichment(genes_by_pathway = tmp_gset,
                               genes_of_interest = c("ACLY", "IDH1", "IDH3G",
                                                     "DLD", "GLUD1"),
                               pathways_list = tmp_gset_names,
                               adj_method = "bonferroni",
                               enrichment_threshold = 0.05,
                               DEG_vec = c("IDH1", "GLUD1"),
                               all_genes = tmp_gset_genes),
            "data.frame")

  # higher enrichment threshold
  expect_is(tmp2 <- enrichment(genes_by_pathway = tmp_gset,
                               genes_of_interest = c("ACLY", "IDH1", "IDH3G",
                                                     "DLD", "GLUD1"),
                               pathways_list = tmp_gset_names,
                               adj_method = "bonferroni",
                               enrichment_threshold = 0.8,
                               DEG_vec = c("IDH1", "GLUD1"),
                               all_genes = tmp_gset_genes),
            "data.frame")

  expect_true(nrow(tmp2) > nrow(tmp1))

  # no enrichment case
  expect_null(enrichment(genes_by_pathway = tmp_gset,
                         genes_of_interest = c("DLD", "GLUD1"),
                         pathways_list = tmp_gset_names,
                         adj_method = "bonferroni",
                         enrichment_threshold = 0.01,
                         DEG_vec = c("IDH1", "GLUD1"),
                         all_genes = tmp_gset_genes))
})

# enrichment_analyses -----------------------------------------------------
pin_path <- return_pin_path()

test_that("enrichment function returns a data frame"{

  # default
  expect_is(enrichment_analyses(snws = example_active_snws[1:2],
                                input_genes = RA_input$Gene.symbol,
                                gene_sets = "KEGG",
                                pin_path = pin_path,
                                adj_method = "bonferroni",
                                enrichment_threshold = 5e-2,
                                list_active_snw_genes = FALSE),
            "data.frame")

  # Custom gene sets
  expect_is(enrichment_analyses(snws = example_active_snws[1:2],
                                input_genes = RA_input$Gene.symbol,
                                gene_sets = "Custom",
                                custom_genes = kegg_genes,
                                custom_pathways = kegg_pathways,
                                pin_path = pin_path,
                                adj_method = "bonferroni",
                                enrichment_threshold = 5e-2,
                                list_active_snw_genes = FALSE),
            "data.frame")

  # list active snw genes
  expect_is(enr_res <- enrichment_analyses(snws = example_active_snws[3:4],
                                           input_genes = RA_input$Gene.symbol,
                                           gene_sets = "KEGG",
                                           pin_path = pin_path,
                                           adj_method = "bonferroni",
                                           enrichment_threshold = 5e-2,
                                           list_active_snw_genes = TRUE),
            "data.frame")

})

# summarize_enrichment_results --------------------------------------------
test_that("summarize_enrichment_results function returns a data frame"{

  # default
  expect_is(tmp <- summarize_enrichment_results(enr_res), "data.frame")
  expect_equal(ncol(tmp), 6)

  # list active snw genes
  expect_is(tmp <- summarize_enrichment_results(enr_res, TRUE), "data.frame")
  expect_equal(ncol(tmp), 7)
})

