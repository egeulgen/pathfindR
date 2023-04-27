##################################################
## Package: pathfindR
## Script purpose: Unit testing script for
## functions related to enrichment
## Date: Apr 27, 2023
## Author: Ege Ulgen
##################################################

# hyperg_test -------------------------------------------------------------
test_that("`hyperg_test()` returns a p value", {
  expect_is(
    tmp_p <- hyperg_test(
      term_genes = LETTERS[1:10],
      chosen_genes = LETTERS[2:5],
      background_genes = LETTERS
    ),
    "numeric"
  )
  expect_true(tmp_p >= 0 & tmp_p <= 1)

  expect_is(
    tmp_p2 <- hyperg_test(
      term_genes = LETTERS[1:4],
      chosen_genes = LETTERS[3:10],
      background_genes = LETTERS
    ),
    "numeric"
  )
  expect_true(tmp_p2 >= 0 & tmp_p2 <= 1)
  expect_true(tmp_p2 > tmp_p)
})

test_that("`hyperg_test()` argument checks work", {
  expect_error(
    hyperg_test(term_genes = list()),
    "`term_genes` should be a vector"
  )

  expect_error(
    hyperg_test(term_genes = LETTERS, chosen_genes = list()),
    "`chosen_genes` should be a vector"
  )

  expect_error(
    hyperg_test(
      term_genes = LETTERS,
      chosen_genes = LETTERS[1:2],
      background_genes = list()
    ),
    "`background_genes` should be a vector"
  )

  expect_error(
    hyperg_test(
      term_genes = c(LETTERS, LETTERS),
      chosen_genes = LETTERS[1:3],
      background_genes = LETTERS
    ),
    "`term_genes` cannot be larger than `background_genes`!"
  )

  expect_error(
    hyperg_test(
      term_genes = LETTERS[1:10],
      chosen_genes = c(LETTERS, LETTERS),
      background_genes = LETTERS
    ),
    "`chosen_genes` cannot be larger than `background_genes`!"
  )
})

# enrichment --------------------------------------------------------------
tmp_gset_obj <- fetch_gene_set(min_gset_size = 100, max_gset_size = 150)
tmp_gset <- tmp_gset_obj$genes_by_term
tmp_gset_descriptions <- tmp_gset_obj$term_descriptions
tmp_gset_genes <- unlist(tmp_gset)

tmp_input_genes <- sample(unlist(tmp_gset[1:3]), 100)
tmp_sig_vec <- c(sample(tmp_input_genes, 100), sample(unlist(tmp_gset), 400))

test_that("`enrichment()` returns a data frame", {
  # default
  expect_is(
    tmp1 <- enrichment(
      input_genes = tmp_input_genes,
      genes_by_term = tmp_gset,
      term_descriptions = tmp_gset_descriptions,
      adj_method = "bonferroni",
      enrichment_threshold = 0.05,
      sig_genes_vec = tmp_sig_vec,
      background_genes = tmp_gset_genes
    ),
    "data.frame"
  )

  # higher enrichment threshold
  expect_is(
    tmp2 <- enrichment(
      input_genes = tmp_input_genes,
      genes_by_term = tmp_gset,
      term_descriptions = tmp_gset_descriptions,
      adj_method = "bonferroni",
      enrichment_threshold = 1,
      sig_genes_vec = tmp_sig_vec,
      background_genes = tmp_gset_genes
    ),
    "data.frame"
  )

  expect_true(nrow(tmp2) > nrow(tmp1))

  # no enrichment case
  expect_null(
    enrichment(
      input_genes = paste0("X", tmp_input_genes),
      genes_by_term = tmp_gset,
      term_descriptions = tmp_gset_descriptions,
      adj_method = "bonferroni",
      enrichment_threshold = 0.5,
      sig_genes_vec = tmp_sig_vec,
      background_genes = tmp_gset_genes
    )
  )
})

test_that("`enrichment()` argument checks work", {
  ## input genes
  expect_error(
    enrichment(
      input_genes = list(),
      sig_genes_vec = "PER1",
      background_genes = unlist(kegg_genes)
    ),
    "`input_genes` should be a vector of gene symbols"
  )

  ## gene sets data
  expect_error(
    enrichment(
      input_genes = tmp_input_genes,
      genes_by_term = "INVALID",
      sig_genes_vec = tmp_sig_vec,
      background_genes = unlist(kegg_genes)
    ),
    "`genes_by_term` should be a list of term gene sets"
  )
  expect_error(
    enrichment(
      input_genes = tmp_input_genes,
      genes_by_term = list(1:3),
      sig_genes_vec = tmp_sig_vec,
      background_genes = unlist(kegg_genes)
    ),
    "`genes_by_term` should be a named list \\(names are gene set IDs\\)"
  )

  expect_error(
    enrichment(
      input_genes = tmp_input_genes,
      term_descriptions = list(),
      sig_genes_vec = tmp_sig_vec,
      background_genes = unlist(kegg_genes)
    ),
    "`term_descriptions` should be a vector of term gene descriptions"
  )
  expect_error(
    enrichment(
      input_genes = tmp_input_genes,
      term_descriptions = 1:3,
      sig_genes_vec = tmp_sig_vec,
      background_genes = unlist(kegg_genes)
    ),
    "`term_descriptions` should be a named vector \\(names are gene set IDs\\)"
  )

  expect_error(
    enrichment(
      input_genes = tmp_input_genes,
      genes_by_term = list(A = 1:3),
      term_descriptions = c(A = "a", B = "b"),
      sig_genes_vec = tmp_sig_vec,
      background_genes = unlist(kegg_genes)
    ),
    "The lengths of `genes_by_term` and `term_descriptions` should be the same"
  )
  expect_error(
    enrichment(
      input_genes = tmp_input_genes,
      genes_by_term = list(A = 1:3, X = 1:3),
      term_descriptions = c(A = "a", B = "b"),
      sig_genes_vec = tmp_sig_vec,
      background_genes = unlist(kegg_genes)
    ),
    "The names of `genes_by_term` and `term_descriptions` should all be the same"
  )

  ## enrichment threshold
  expect_error(
    enrichment(
      input_genes = tmp_input_genes,
      sig_genes_vec = tmp_sig_vec,
      background_genes = unlist(kegg_genes),
      enrichment_threshold = "INVALID"
    ),
    "`enrichment_threshold` should be a numeric value between 0 and 1"
  )

  expect_error(
    enrichment(
      input_genes = tmp_input_genes,
      sig_genes_vec = tmp_sig_vec,
      background_genes = unlist(kegg_genes),
      enrichment_threshold = -1
    ),
    "`enrichment_threshold` should be between 0 and 1"
  )

  ## signif. genes and background (universal set) genes
  expect_error(
    enrichment(
      input_genes = tmp_input_genes,
      sig_genes_vec = list(),
      background_genes = unlist(kegg_genes)
    ),
    "`sig_genes_vec` should be a vector"
  )
  expect_error(
    enrichment(
      input_genes = tmp_input_genes,
      sig_genes_vec = tmp_sig_vec,
      background_genes = list()
    ),
    "`background_genes` should be a vector"
  )
})


# enrichment_analyses -----------------------------------------------------
tmp_gset_genes <- tmp_gset_obj$genes_by_term
tmp_gset_desc <- tmp_gset_obj$term_descriptions


test_that("`enrichment_analyses()` returns a data frame", {
  # default
  expect_is(
    tmp1 <- enrichment_analyses(
      snws = example_active_snws[1:3],
      sig_genes_vec = example_pathfindR_input$Gene.symbol,
      pin_name_path = "Biogrid",
      genes_by_term = tmp_gset_genes,
      term_descriptions = tmp_gset_desc,
      adj_method = "bonferroni",
      enrichment_threshold = 0.05,
      list_active_snw_genes = FALSE
    ),
    "data.frame"
  )

  # list active snw genes
  expect_is(
    tmp2 <- enrichment_analyses(
      snws = example_active_snws[1:3],
      sig_genes_vec = example_pathfindR_input$Gene.symbol,
      pin_name_path = "Biogrid",
      genes_by_term = tmp_gset_genes,
      term_descriptions = tmp_gset_desc,
      adj_method = "bonferroni",
      enrichment_threshold = 0.05,
      list_active_snw_genes = TRUE
    ),
    "data.frame"
  )
  expect_true(ncol(tmp2) == ncol(tmp1) + 1)
})

test_that("`enrichment_analyses()` argument check works", {
  expect_error(
    enrichment_analyses(
      snws = example_active_snws,
      list_active_snw_genes = "INVALID"
    ),
    "`list_active_snw_genes` should be either TRUE or FALSE"
  )
})

# summarize_enrichment_results --------------------------------------------
iter1_res <- enrichment_analyses(
  snws = example_active_snws[1:5],
  sig_genes_vec = example_pathfindR_input$Gene.symbol,
  genes_by_term = tmp_gset_genes,
  term_descriptions = tmp_gset_desc,
  list_active_snw_genes = TRUE
)

iter2_res <- enrichment_analyses(
  snws = example_active_snws[11:16],
  sig_genes_vec = example_pathfindR_input$Gene.symbol,
  genes_by_term = tmp_gset_genes,
  term_descriptions = tmp_gset_desc,
  list_active_snw_genes = TRUE
)
combined_res <- rbind(iter1_res, iter2_res)

test_that("`summarize_enrichment_results()` returns summarized enrichment results", {
  # default
  expect_is(
    tmp <- summarize_enrichment_results(enrichment_res = combined_res[, -6]),
    "data.frame"
  )
  expect_equal(ncol(tmp), 7)
  expect_false("non_Signif_Snw_Genes" %in% colnames(tmp))
  expect_true(nrow(tmp) <= nrow(combined_res))

  # list active snw genes
  expect_is(
    tmp <- summarize_enrichment_results(
      enrichment_res = combined_res,
      list_active_snw_genes = TRUE
    ),
    "data.frame"
  )
  expect_equal(ncol(tmp), 8)
  expect_true("non_Signif_Snw_Genes" %in% colnames(tmp))
  expect_true(nrow(tmp) <= nrow(combined_res))
})

test_that("summarize_enrichment_results() argument checks work", {
  expect_error(
    summarize_enrichment_results(
      enrichment_res = combined_res,
      list_active_snw_genes = "INVALID"
    ),
    "`list_active_snw_genes` should be either TRUE or FALSE"
  )

  expect_error(
    summarize_enrichment_results(enrichment_res = list()),
    "`enrichment_res` should be a data frame"
  )

  # list_active_snw_genes = FALSE
  nec_cols <- c("ID", "Term_Description", "Fold_Enrichment", "p_value", "adj_p", "support")

  expect_error(
    summarize_enrichment_results(enrichment_res = data.frame()),
    paste0("`enrichment_res` should have exactly ", length(nec_cols), " columns")
  )

  tmp <- as.data.frame(matrix(
    nrow = 1, ncol = length(nec_cols),
    dimnames = list(NULL, letters[seq_along(nec_cols)])
  ))
  expect_error(
    summarize_enrichment_results(enrichment_res = tmp),
    paste0(
      "`enrichment_res` should have column names ",
      paste(dQuote(nec_cols), collapse = ", ")
    )
  )

  # list_active_snw_genes = TRUE
  nec_cols <- c(
    "ID", "Term_Description", "Fold_Enrichment", "p_value", "adj_p", "support",
    "non_Signif_Snw_Genes"
  )

  expect_error(
    summarize_enrichment_results(
      enrichment_res = data.frame(),
      list_active_snw_genes = TRUE
    ),
    paste0("`enrichment_res` should have exactly ", length(nec_cols), " columns")
  )

  tmp <- as.data.frame(matrix(
    nrow = 1, ncol = length(nec_cols),
    dimnames = list(NULL, letters[seq_along(nec_cols)])
  ))
  expect_error(
    summarize_enrichment_results(
      enrichment_res = tmp,
      list_active_snw_genes = TRUE
    ),
    paste0(
      "`enrichment_res` should have column names ",
      paste(dQuote(nec_cols), collapse = ", ")
    )
  )
})
