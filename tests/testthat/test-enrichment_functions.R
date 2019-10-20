##################################################
## Project: pathfindR
## Script purpose: Testthat testing script for
## enrichment-related functions
## Date: Oct 20, 2019
## Author: Ege Ulgen
##################################################

# hyperg_test -------------------------------------------------------------
test_that("hyperg_test returns a p value", {
  expect_is(tmp_p <- hyperg_test(term_genes = LETTERS[1:10],
                                 chosen_genes = LETTERS[2:5],
                                 background_genes = LETTERS),
            "numeric")
  expect_true(tmp_p >= 0 & tmp_p <= 1)

  expect_is(tmp_p2 <- hyperg_test(term_genes = LETTERS[1:4],
                                  chosen_genes = LETTERS[3:10],
                                  background_genes = LETTERS),
            "numeric")
  expect_true(tmp_p2 >= 0 & tmp_p2 <= 1)
  expect_true(tmp_p2 > tmp_p)
})

test_that("hyperg_test arg checks work", {
  expect_error(hyperg_test(term_genes = c(LETTERS, LETTERS),
                           chosen_genes = LETTERS[1:3],
                           background_genes = LETTERS),
               "`term_genes` cannot be larger than `background_genes`!"
  )
  expect_error(hyperg_test(term_genes = LETTERS[1:10],
                           chosen_genes = c(LETTERS, LETTERS),
                           background_genes = LETTERS),
               "`chosen_genes` cannot be larger than `background_genes`!"
  )
})

# enrichment --------------------------------------------------------------
tmp_gset_obj <- fetch_gene_set(min_gset_size = 100, max_gset_size = 150)
tmp_gset <- tmp_gset_obj$genes_by_term
tmp_gset_descriptions <- tmp_gset_obj$term_descriptions
tmp_gset_genes <- unlist(tmp_gset)

test_that("enrichment function returns a data frame", {

  tmp_input_genes <- sample(unlist(tmp_gset[1:3]), 100)
  tmp_sig_vec <- c(sample(tmp_input_genes, 100), sample(unlist(tmp_gset), 400))

  # default
  expect_is(tmp1 <- enrichment(input_genes = tmp_input_genes,
                               genes_by_term = tmp_gset,
                               term_descriptions = tmp_gset_descriptions,
                               adj_method = "bonferroni",
                               enrichment_threshold = 0.05,
                               sig_genes_vec = tmp_sig_vec,
                               background_genes = tmp_gset_genes),
    "data.frame")

  # higher enrichment threshold
  expect_is(tmp2 <- enrichment(input_genes = tmp_input_genes,
                               genes_by_term = tmp_gset,
                               term_descriptions = tmp_gset_descriptions,
                               adj_method = "bonferroni",
                               enrichment_threshold = 0.5,
                               sig_genes_vec = tmp_sig_vec,
                               background_genes = tmp_gset_genes),
            "data.frame")

  expect_true(nrow(tmp2) > nrow(tmp1))

  # no enrichment case
  expect_null(enrichment(input_genes = paste0("X", tmp_input_genes),
                         genes_by_term = tmp_gset,
                         term_descriptions = tmp_gset_descriptions,
                         adj_method = "bonferroni",
                         enrichment_threshold = 0.5,
                         sig_genes_vec = tmp_sig_vec,
                         background_genes = tmp_gset_genes))
})

# enrichment_analyses -----------------------------------------------------
pin_path <- return_pin_path()
tmp <- fetch_gene_set()
tmp_gset_genes <- tmp$genes_by_term
tmp_gset_desc <- tmp$term_descriptions

test_that("enrichment function returns a data frame", {

  # default
  expect_is(tmp1 <- enrichment_analyses(snws = example_active_snws[1:3],
                                        sig_genes_vec = RA_input$Gene.symbol,
                                        pin_path = pin_path,
                                        genes_by_term = tmp_gset_genes,
                                        term_descriptions = tmp_gset_desc,
                                        adj_method = "bonferroni",
                                        enrichment_threshold = 0.05,
                                        list_active_snw_genes = FALSE),
            "data.frame")

  # list active snw genes
  expect_is(tmp2 <- enrichment_analyses(snws = example_active_snws[1:3],
                                        sig_genes_vec = RA_input$Gene.symbol,
                                        pin_path = pin_path,
                                        genes_by_term = tmp_gset_genes,
                                        term_descriptions = tmp_gset_desc,
                                        adj_method = "bonferroni",
                                        enrichment_threshold = 0.05,
                                        list_active_snw_genes = TRUE),
            "data.frame")
  expect_true(ncol(tmp2) == ncol(tmp1) + 1)
})

# summarize_enrichment_results --------------------------------------------
enr_res <- enrichment_analyses(snws = example_active_snws[1:10],
                               sig_genes_vec = RA_input$Gene.symbol,
                               pin_path = pin_path,
                               genes_by_term = tmp_gset_genes,
                               term_descriptions = tmp_gset_desc,
                               adj_method = "bonferroni",
                               enrichment_threshold = 0.05,
                               list_active_snw_genes = TRUE)

test_that("summarize_enrichment_results function returns a data frame", {
  # default
  expect_is(tmp <- summarize_enrichment_results(enrichment_res = enr_res),
            "data.frame")
  expect_equal(ncol(tmp), 6)
  expect_false("non_DEG_Active_Snw_Genes" %in% colnames(tmp))

  # list active snw genes
  expect_is(tmp <- summarize_enrichment_results(enr_res, TRUE), "data.frame")
  expect_equal(ncol(tmp), 7)
  expect_true("non_DEG_Active_Snw_Genes" %in% colnames(tmp))
  expect_true(nrow(tmp) <= nrow(enr_res))
})
