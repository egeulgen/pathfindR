## Tests for functions related to enrichment analyses - Aug 2023
set.seed(123)

test_that("`hyperg_test()` -- returns an appropriate p value", {
    expect_is(tmp_p <- hyperg_test(term_genes = LETTERS[1:10], chosen_genes = LETTERS[2:5],
        background_genes = LETTERS), "numeric")
    expect_true(tmp_p >= 0 & tmp_p <= 1)

    expect_is(tmp_p2 <- hyperg_test(term_genes = LETTERS[1:4], chosen_genes = LETTERS[3:10],
        background_genes = LETTERS), "numeric")
    expect_true(tmp_p2 >= 0 & tmp_p2 <= 1)
    expect_true(tmp_p2 > tmp_p)
})

test_that("`hyperg_test()` -- argument checks work", {
    expect_error(hyperg_test(term_genes = list()), "`term_genes` should be a vector")
    expect_error(hyperg_test(term_genes = LETTERS, chosen_genes = list()), "`chosen_genes` should be a vector")
    expect_error(hyperg_test(term_genes = LETTERS, chosen_genes = LETTERS[1:2], background_genes = list()),
        "`background_genes` should be a vector")
    expect_error(hyperg_test(term_genes = c(LETTERS, LETTERS), chosen_genes = LETTERS[1:3],
        background_genes = LETTERS), "`term_genes` cannot be larger than `background_genes`!")
    expect_error(hyperg_test(term_genes = LETTERS[1:10], chosen_genes = c(LETTERS,
        LETTERS), background_genes = LETTERS), "`chosen_genes` cannot be larger than `background_genes`!")
})

test_that("`enrichment()` -- returns a data frame", {
    expected_num_significant <- 10
    gsets <- example_pathfindR_output$ID[1:50]
    p_val_vec <- c(runif(expected_num_significant, min = 1e-05, max = 0.001), runif(length(gsets) -
        expected_num_significant, min = 0.05, max = 1))
    names(p_val_vec) <- gsets
    mock_vapply <- mockery::mock(p_val_vec, 5, 2, cycle = TRUE)
    mockery::stub(enrichment, "vapply", mock_vapply)
    mockery::stub(enrichment, "base::setdiff", c("RPS6KA2", "HSPA2", "SCN4B", "PPP2R1B",
        "PTCH1", "CASP10", "TIRAP", "BEX3", "KIF5C", "TNFSF13B"))

    # default
    expect_is(enr_res <- enrichment(input_genes = example_pathfindR_input$Gene.symbol,
        sig_genes_vec = c("DummyGene"), background_genes = c("DummyGene")), "data.frame")
    expect_equal(nrow(enr_res), expected_num_significant)
    expect_true(any(enr_res$non_Signif_Snw_Genes != ""))
    expect_true(all(enr_res$Fold_Enrichment == 2.5))

    # higher threshold - no filter
    expect_is(enr_res2 <- enrichment(input_genes = example_pathfindR_input$Gene.symbol,
        sig_genes_vec = c("DummyGene"), background_genes = c("DummyGene"), enrichment_threshold = 1), "data.frame")
    expect_equal(nrow(enr_res2), 50)
    expect_true(any(enr_res2$non_Signif_Snw_Genes != ""))

    # no enrichment case
    mockery::stub(enrichment, "stats::p.adjust", rep(1, 50))
    expect_null(enr_res3 <- enrichment(input_genes = example_pathfindR_input$Gene.symbol,
        sig_genes_vec = c("DummyGene"), background_genes = c("DummyGene")))
})

test_that("`enrichment()` -- argument checks work", {
    tmp_input_genes <- example_pathfindR_input$Gene.symbol[1:6]
    tmp_sig_vec <- example_pathfindR_input$Gene.symbol[1:3]
    ## input genes
    expect_error(enrichment(input_genes = list(), sig_genes_vec = "PER1", background_genes = unlist(kegg_genes)),
        "`input_genes` should be a vector of gene symbols")

    ## gene sets data
    expect_error(enrichment(input_genes = tmp_input_genes, genes_by_term = "INVALID",
        sig_genes_vec = tmp_sig_vec, background_genes = unlist(kegg_genes)), "`genes_by_term` should be a list of term gene sets")
    expect_error(enrichment(input_genes = tmp_input_genes, genes_by_term = list(1:3),
        sig_genes_vec = tmp_sig_vec, background_genes = unlist(kegg_genes)), "`genes_by_term` should be a named list \\(names are gene set IDs\\)")

    expect_error(enrichment(input_genes = tmp_input_genes, term_descriptions = list(),
        sig_genes_vec = tmp_sig_vec, background_genes = unlist(kegg_genes)), "`term_descriptions` should be a vector of term gene descriptions")
    expect_error(enrichment(input_genes = tmp_input_genes, term_descriptions = 1:3,
        sig_genes_vec = tmp_sig_vec, background_genes = unlist(kegg_genes)), "`term_descriptions` should be a named vector \\(names are gene set IDs\\)")

    expect_error(enrichment(input_genes = tmp_input_genes, genes_by_term = list(A = 1:3),
        term_descriptions = c(A = "a", B = "b"), sig_genes_vec = tmp_sig_vec, background_genes = unlist(kegg_genes)),
        "The lengths of `genes_by_term` and `term_descriptions` should be the same")
    expect_error(enrichment(input_genes = tmp_input_genes, genes_by_term = list(A = 1:3,
        X = 1:3), term_descriptions = c(A = "a", B = "b"), sig_genes_vec = tmp_sig_vec,
        background_genes = unlist(kegg_genes)), "The names of `genes_by_term` and `term_descriptions` should all be the same")

    ## enrichment threshold
    expect_error(enrichment(input_genes = tmp_input_genes, sig_genes_vec = tmp_sig_vec,
        background_genes = unlist(kegg_genes), enrichment_threshold = "INVALID"),
        "`enrichment_threshold` should be a numeric value between 0 and 1")

    expect_error(enrichment(input_genes = tmp_input_genes, sig_genes_vec = tmp_sig_vec,
        background_genes = unlist(kegg_genes), enrichment_threshold = -1), "`enrichment_threshold` should be between 0 and 1")

    ## signif. genes and background (universal set) genes
    expect_error(enrichment(input_genes = tmp_input_genes, sig_genes_vec = list(),
        background_genes = unlist(kegg_genes)), "`sig_genes_vec` should be a vector")
    expect_error(enrichment(input_genes = tmp_input_genes, sig_genes_vec = tmp_sig_vec,
        background_genes = list()), "`background_genes` should be a vector")
})

tmp_gset_genes <- kegg_genes[example_pathfindR_output$ID[order(example_pathfindR_output$support,
    decreasing = TRUE)[1:10]]]
tmp_gset_desc <- kegg_descriptions[names(tmp_gset_genes)]

all_iter_enr_res <- list(NULL, NULL, NULL)
subnw_start_idx <- c(1, 4, 7)
for (idx in seq_along(subnw_start_idx)) {
    j <- subnw_start_idx[idx]
    res <- enrichment_analyses(snws = example_active_snws[j:j + 2], sig_genes_vec = example_pathfindR_input$Gene.symbol,
        genes_by_term = tmp_gset_genes, term_descriptions = tmp_gset_desc, list_active_snw_genes = TRUE)
    if (!is.null(res)) {
        all_iter_enr_res[[idx]] <- res
    }
}

combined_res <- do.call(rbind, all_iter_enr_res)

test_that("`enrichment_analyses()` -- returns a data frame", {
    toy_pin <- data.frame(V1 = paste("Gene", sample(1:50, 10)), V2 = "pp", V3 = paste("Gene",
        sample(1:50, 10)))
    mockery::stub(enrichment_analyses, "return_pin_path", NULL)
    mockery::stub(enrichment_analyses, "utils::read.delim", toy_pin)

    mock_lapply <- mockery::mock(c(), all_iter_enr_res, cycle = TRUE)
    mockery::stub(enrichment_analyses, "lapply", mock_lapply)

    # default
    expect_is(enr_res1 <- enrichment_analyses(snws = example_active_snws[1:3], sig_genes_vec = example_pathfindR_input$Gene.symbol,
        list_active_snw_genes = FALSE), "data.frame")
    total <- sum(vapply(all_iter_enr_res, function(x) ifelse(is.null(x), 0, nrow(x)),
        1))
    expect_true(nrow(enr_res1) <= total)

    # list active snw genes
    expect_is(enr_res2 <- enrichment_analyses(snws = example_active_snws[1:3], sig_genes_vec = example_pathfindR_input$Gene.symbol,
        list_active_snw_genes = TRUE), "data.frame")
    expect_true(ncol(enr_res2) == ncol(enr_res1) + 1)
})

test_that("`enrichment_analyses()` -- argument check works", {
    expect_error(enrichment_analyses(snws = example_active_snws, list_active_snw_genes = "INVALID"),
        "`list_active_snw_genes` should be either TRUE or FALSE")
})

test_that("`summarize_enrichment_results()` -- returns summarized enrichment results",
    {
        # default
        expect_is(summ_res <- summarize_enrichment_results(enrichment_res = combined_res[,
            -6]), "data.frame")
        expect_equal(ncol(summ_res), 7)
        expect_false("non_Signif_Snw_Genes" %in% colnames(summ_res))
        expect_true(nrow(summ_res) <= nrow(combined_res))

        # list active snw genes
        expect_is(summ_res2 <- summarize_enrichment_results(enrichment_res = combined_res,
            list_active_snw_genes = TRUE), "data.frame")
        expect_equal(ncol(summ_res2), 8)
        expect_true("non_Signif_Snw_Genes" %in% colnames(summ_res2))
        expect_true(nrow(summ_res2) <= nrow(combined_res))
    })

test_that("`summarize_enrichment_results()` -- argument checks work", {
    expect_error(summarize_enrichment_results(enrichment_res = combined_res, list_active_snw_genes = "INVALID"),
        "`list_active_snw_genes` should be either TRUE or FALSE")

    expect_error(summarize_enrichment_results(enrichment_res = list()), "`enrichment_res` should be a data frame")

    # list_active_snw_genes = FALSE
    nec_cols <- c("ID", "Term_Description", "Fold_Enrichment", "p_value", "adj_p",
        "support")

    expect_error(summarize_enrichment_results(enrichment_res = data.frame()), paste0("`enrichment_res` should have exactly ",
        length(nec_cols), " columns"))

    tmp <- as.data.frame(matrix(nrow = 1, ncol = length(nec_cols), dimnames = list(NULL,
        letters[seq_along(nec_cols)])))
    expect_error(summarize_enrichment_results(enrichment_res = tmp), paste0("`enrichment_res` should have column names ",
        paste(dQuote(nec_cols), collapse = ", ")))

    # list_active_snw_genes = TRUE
    nec_cols <- c("ID", "Term_Description", "Fold_Enrichment", "p_value", "adj_p",
        "support", "non_Signif_Snw_Genes")

    expect_error(summarize_enrichment_results(enrichment_res = data.frame(), list_active_snw_genes = TRUE),
        paste0("`enrichment_res` should have exactly ", length(nec_cols), " columns"))

    tmp <- as.data.frame(matrix(nrow = 1, ncol = length(nec_cols), dimnames = list(NULL,
        letters[seq_along(nec_cols)])))
    expect_error(summarize_enrichment_results(enrichment_res = tmp, list_active_snw_genes = TRUE),
        paste0("`enrichment_res` should have column names ", paste(dQuote(nec_cols),
            collapse = ", ")))
})
