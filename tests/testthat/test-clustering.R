## Tests for enriched term clustering functions - Aug 2023

enrichment_res <- example_pathfindR_output[1:5, ]
input_kappa_mat <- create_kappa_matrix(enrichment_res)

test_that("`create_kappa_matrix()` -- creates kappa matrix", {
    input_df <- enrichment_res
    kappa_mat <- create_kappa_matrix(input_df)
    expect_true(isSymmetric.matrix(kappa_mat))
    expect_true(all(kappa_mat >= 0 & kappa_mat <= 1 | kappa_mat >= -1 & kappa_mat <=
        0))
    expect_identical(colnames(kappa_mat), rownames(kappa_mat))
    expect_identical(colnames(kappa_mat), input_df$ID)

    # zero length excluded
    input_df2 <- input_df
    input_df2$Down_regulated[1] <- input_df2$Up_regulated[1] <- ""
    kappa_mat2 <- create_kappa_matrix(input_df2)
    expect_true(isSymmetric.matrix(kappa_mat2))
    expect_false(input_df2$ID[1] %in% colnames(kappa_mat2))

    input_df$non_Signif_Snw_Genes <- c("GeneA, GeneB", "GeneA", "GeneC", "GeneB, GeneC",
        "")
    kappa_mat3 <- create_kappa_matrix(input_df, use_active_snw_genes = TRUE)
    expect_true(isSymmetric.matrix(kappa_mat3))
    expect_true(!all(kappa_mat3 != kappa_mat))
})

test_that("`create_kappa_matrix()` -- argument checks works", {
    expect_error(create_kappa_matrix(example_pathfindR_output, use_description = "INVALID"),
        "`use_description` should be TRUE or FALSE")
    expect_error(create_kappa_matrix(example_pathfindR_output, use_active_snw_genes = "INVALID"),
        "`use_active_snw_genes` should be TRUE or FALSE")
    expect_error(create_kappa_matrix(list()), "`enrichment_res` should be a data frame of enrichment results")
    expect_error(create_kappa_matrix(example_pathfindR_output[1, ]), "`enrichment_res` should contain at least 2 rows")

    cr_cols <- function(use_description = FALSE, use_active_snw_genes = FALSE) {
        nec_cols <- c("Down_regulated", "Up_regulated")
        if (use_description) {
            nec_cols <- c("Term_Description", nec_cols)
        } else {
            nec_cols <- c("ID", nec_cols)
        }
        if (use_active_snw_genes) {
            nec_cols <- c(nec_cols, "non_Signif_Snw_Genes")
        }
        return(nec_cols)
    }

    # desc F
    nec_cols <- cr_cols()
    valid_res <- enrichment_res[, -2]
    expect_silent(create_kappa_matrix(valid_res))
    invalid_res <- enrichment_res[, -1]
    expect_error(create_kappa_matrix(invalid_res), paste0("`enrichment_res` should contain all of ",
        paste(dQuote(nec_cols), collapse = ", ")))
    # desc T
    nec_cols <- cr_cols(use_description = TRUE)
    valid_res <- enrichment_res[, -1]
    expect_silent(create_kappa_matrix(valid_res, use_description = TRUE))
    invalid_res <- enrichment_res[, -2]
    expect_error(create_kappa_matrix(invalid_res, use_description = TRUE), paste0("`enrichment_res` should contain all of ",
        paste(dQuote(nec_cols), collapse = ", ")))
    # snw_g T
    nec_cols <- cr_cols(use_active_snw_genes = TRUE)
    valid_res <- enrichment_res
    valid_res$non_Signif_Snw_Genes <- ""
    expect_silent(create_kappa_matrix(valid_res, use_active_snw_genes = TRUE))
    expect_error(create_kappa_matrix(enrichment_res, use_active_snw_genes = TRUE),
        paste0("`enrichment_res` should contain all of ", paste(dQuote(nec_cols),
            collapse = ", ")))
})

test_that("`hierarchical_term_clustering()` -- returns integer vector", {
    m <- mockery::mock(NULL, cycle = TRUE)
    mockery::stub(hierarchical_term_clustering, "graphics::plot", m)
    mockery::stub(hierarchical_term_clustering, "stats::heatmap", m)
    mockery::stub(hierarchical_term_clustering, "stats::rect.hclust", m)

    expected_message_regex <- "The maximum average silhouette width was -?(0\\.?\\d{0,2}|1) for k = \\d+ \n\n"
    expect_message(clu_res <- hierarchical_term_clustering(input_kappa_mat, enrichment_res,
        plot_hmap = TRUE, plot_dend = TRUE), expected_message_regex)
    expect_is(clu_res, "integer")
    expect_true(max(clu_res) <= nrow(input_kappa_mat))
    expect_identical(rownames(input_kappa_mat), names(clu_res))
})

test_that("`hierarchical_term_clustering()` -- `num_clusters` works", {
    for (selected_num_clusters in seq_len(nrow(enrichment_res))) {
        expect_is(res <- hierarchical_term_clustering(input_kappa_mat, enrichment_res,
            num_clusters = selected_num_clusters, plot_hmap = FALSE, plot_dend = FALSE),
            "integer")
        expect_equal(max(res), selected_num_clusters)
    }
})

test_that("`hierarchical_term_clustering()` -- `kseq` (sequence of number of clusters to try) is determined appropriately",
    {
        mockery::stub(hierarchical_term_clustering, "stats::hclust", NULL)
        mockery::stub(hierarchical_term_clustering, "isSymmetric.matrix", TRUE)

        mock_cutree <- function(tree, k, h = NULL) {
            return(k)
        }
        mockery::stub(hierarchical_term_clustering, "stats::cutree", mock_cutree)

        for (num_terms in c(3, 15, 153, 200, 204, 432)) {
            kmax <- max(num_terms%/%2, 2)
            num_expected_calls <- ifelse(kmax <= 20, kmax - 1, ifelse(kmax <= 100,
                18 + kmax%/%10 - 1, 26 + kmax%/%50 - 1))
            target_k <- ifelse(kmax <= 20, kmax, ifelse(kmax <= 100, round(kmax%/%10) *
                10, round(kmax%/%50) * 50))

            tmp_enr_res <- example_pathfindR_output[seq_len(num_terms), ]
            tmp_kappa_mat <- matrix(NA, nrow = num_terms, ncol = num_terms, dimnames = list(tmp_enr_res$ID,
                tmp_enr_res$ID))

            silwidth_out_vec <- vector("list", num_expected_calls)
            for (idx in seq_len(num_expected_calls)) {
                if (idx == length(silwidth_out_vec)) {
                  silwidth_out_vec[[idx]] <- list(avg.silwidth = 100)
                } else {
                  silwidth_out_vec[[idx]] <- list(avg.silwidth = -100)
                }
            }
            mock_cluster.stats <- do.call(mockery::mock, silwidth_out_vec)
            mockery::stub(hierarchical_term_clustering, "fpc::cluster.stats", mock_cluster.stats)

            expected_message <- paste0("The maximum average silhouette width was 100 for k = ",
                target_k, " \n\n")
            expect_message(res_k <- hierarchical_term_clustering(tmp_kappa_mat, tmp_enr_res,
                plot_hmap = FALSE, plot_dend = FALSE), expected_message)
            expect_equal(res_k, target_k)
            mockery::expect_called(mock_cluster.stats, num_expected_calls)
        }
    })

test_that("`hierarchical_term_clustering()` -- argument checks work", {
    expect_error(hierarchical_term_clustering(kappa_mat = list(), enrichment_res = data.frame()),
        "`kappa_mat` should be a symmetric matrix")
    expect_error(hierarchical_term_clustering(kappa_mat = matrix(nrow = 1, ncol = 2),
        enrichment_res = data.frame()), "`kappa_mat` should be a symmetric matrix")

    mat <- matrix(nrow = 3, ncol = 3, dimnames = list(1:3, 1:3))
    expect_error(hierarchical_term_clustering(kappa_mat = mat, enrichment_res = data.frame(ID = 4:5)),
        "All terms in `kappa_mat` should be present in `enrichment_res`")
    expect_error(hierarchical_term_clustering(kappa_mat = mat, enrichment_res = data.frame(ID = 1:3),
        plot_hmap = "INVALID"), "`plot_hmap` should be TRUE or FALSE")
    expect_error(hierarchical_term_clustering(kappa_mat = mat, enrichment_res = data.frame(ID = 1:3),
        plot_dend = "INVALID"), "`plot_dend` should be TRUE or FALSE")
})

test_that("`fuzzy_term_clustering()` -- returns matrix of cluster memberships", {
    expect_is(res_mat <- fuzzy_term_clustering(create_kappa_matrix(example_pathfindR_output[1:25, ]), example_pathfindR_output[1:25, ], kappa_threshold = 0.1),
        "matrix")
    expect_true(is.logical(res_mat))
})

test_that("`fuzzy_term_clustering()` -- argument checks work", {
    expect_error(fuzzy_term_clustering(kappa_mat = list(), enrichment_res = data.frame()),
        "`kappa_mat` should be a symmetric matrix")
    expect_error(fuzzy_term_clustering(kappa_mat = matrix(nrow = 1, ncol = 2), enrichment_res = data.frame()),
        "`kappa_mat` should be a symmetric matrix")

    mat <- matrix(nrow = 3, ncol = 3, dimnames = list(1:3, 1:3))
    expect_error(fuzzy_term_clustering(kappa_mat = mat, enrichment_res = data.frame(ID = 4:5)),
        "All terms in `kappa_mat` should be present in `enrichment_res`")
    expect_error(fuzzy_term_clustering(kappa_mat = mat, enrichment_res = data.frame(ID = 1:3),
        kappa_threshold = "INVALID"), "`kappa_threshold` should be numeric")
    expect_error(fuzzy_term_clustering(kappa_mat = mat, enrichment_res = data.frame(ID = 1:3),
        kappa_threshold = 1.5), "`kappa_threshold` should be at most 1 as kappa statistic is always <= 1")
})

test_that("`cluster_graph_vis()` -- graph visualization of clusters works OK", {
    mockery::stub(hierarchical_term_clustering, "graphics::plot", NULL)
    mockery::stub(hierarchical_term_clustering, "stats::rect.hclust", NULL)
    mock_plot.igraph <- mockery::mock(NULL, cycle = TRUE)
    mockery::stub(cluster_graph_vis, "igraph::plot.igraph", mock_plot.igraph)
    ## use_description = FALSE
    for (clustering_func in c(hierarchical_term_clustering, fuzzy_term_clustering)) {
        clu_obj <- clustering_func(input_kappa_mat, enrichment_res)
        expect_silent(cluster_graph_vis(clu_obj, input_kappa_mat, enrichment_res))
    }
    mockery::expect_called(mock_plot.igraph, 2)
})

test_that("`cluster_graph_vis()` -- coloring of 'extra' clusters work", {
    mockery::stub(cluster_graph_vis, "igraph::plot.igraph", NULL)
    ### more than 41 clusters (number of colors available)
    selected_num_terms <- 45
    clu_input_df <- example_pathfindR_output[seq_len(selected_num_terms), ]
    mock_kappa_mat <- matrix(NA, nrow = selected_num_terms, ncol = selected_num_terms,
        dimnames = list(clu_input_df$ID, clu_input_df$ID))

    # dummy hierarchical result
    hierarchical_clu_obj <- seq_len(selected_num_terms)
    names(hierarchical_clu_obj) <- clu_input_df$ID
    expect_silent(cluster_graph_vis(hierarchical_clu_obj, mock_kappa_mat, clu_input_df))

    # dummy fuzzy result
    fuzzy_clu_obj <- matrix(FALSE, nrow = selected_num_terms, ncol = selected_num_terms,
        dimnames = list(clu_input_df$ID, seq_len(selected_num_terms)))
    diag(fuzzy_clu_obj) <- TRUE
    expect_silent(cluster_graph_vis(fuzzy_clu_obj, mock_kappa_mat, clu_input_df))
})

test_that("`cluster_graph_vis()` -- check errors are raised appropriately", {
    expect_error(cluster_graph_vis(list(), matrix(), data.frame(ID = 1)), "Invalid class for `clu_obj`!")

    # hierarchical - missing terms in kappa matrix
    clu_obj <- hierarchical_term_clustering(input_kappa_mat, enrichment_res, plot_dend = FALSE)
    expect_error(cluster_graph_vis(c(clu_obj, EXTRA = 1L), input_kappa_mat, enrichment_res),
        "Not all terms in `clu_obj` present in `kappa_mat`!")

    # fuzzy - missing terms in kappa matrix
    clu_obj <- fuzzy_term_clustering(input_kappa_mat, enrichment_res)
    expect_error(cluster_graph_vis(rbind(clu_obj, EXTRA = rep(FALSE, ncol(clu_obj))),
        input_kappa_mat, enrichment_res), "Not all terms in `clu_obj` present in `kappa_mat`!")
})

test_that("`cluster_enriched_terms()` -- returns the input data frame
          with the additional columns `Cluster` and `Status`",
    {
        set.seed(123)
        num_clusters <- 3
        available_clus <- seq_len(num_clusters)
        toy_hierarchical_clu_obj <- sample(available_clus, size = nrow(enrichment_res) -
            1, replace = TRUE)
        missing_clu <- setdiff(available_clus, toy_hierarchical_clu_obj)
        toy_hierarchical_clu_obj <- c(toy_hierarchical_clu_obj, missing_clu)
        toy_hierarchical_clu_obj <- sample(toy_hierarchical_clu_obj)
        names(toy_hierarchical_clu_obj) <- enrichment_res$ID

        toy_fuzzy_clu_obj <- matrix(FALSE, nrow = nrow(enrichment_res), ncol = num_clusters,
            dimnames = list(enrichment_res$ID, available_clus))
        for (row_idx in seq_len(nrow(toy_fuzzy_clu_obj))) {
            num_memberships <- sample(available_clus, 1)
            new_cols <- sample(available_clus, size = num_memberships)
            toy_fuzzy_clu_obj[row_idx, new_cols] <- TRUE
        }

        mock_doCall <- function(...) {
            arguments <- list(...)
            if (arguments[1] == "hierarchical_term_clustering") {
                return(toy_hierarchical_clu_obj)
            }
            if (arguments[1] == "fuzzy_term_clustering") {
                return(toy_fuzzy_clu_obj)
            }
            return(NULL)
        }

        mockery::stub(cluster_enriched_terms, "create_kappa_matrix", input_kappa_mat)
        mockery::stub(cluster_enriched_terms, "R.utils::doCall", mock_doCall)

        # hierarchical
        expect_is(h_clu_res <- cluster_enriched_terms(enrichment_res), "data.frame")
        expect_true(all(c("Cluster", "Status") %in% colnames(h_clu_res)))
        expect_equal(max(h_clu_res$Cluster), num_clusters)
        # expect to have same number of rep. terms as the number of clusters
        expect_equal(max(h_clu_res$Cluster), sum(h_clu_res$Status == "Representative"))

        ## fuzzy
        expect_is(fuzzy_clu_res <- cluster_enriched_terms(enrichment_res, method = "fuzzy"),
            "data.frame")
        expect_true(all(c("Cluster", "Status") %in% colnames(fuzzy_clu_res)))
        expect_true(max(fuzzy_clu_res$Cluster) <= sum(fuzzy_clu_res$Status == "Representative"))
    })

test_that("`cluster_enriched_terms()` argument checks work", {
    expect_error(cluster_enriched_terms(enrichment_res, method = "INVALID"), "the clustering `method` must either be \"hierarchical\" or \"fuzzy\"")
    expect_error(cluster_enriched_terms(enrichment_res, plot_clusters_graph = "INVALID"),
        "`plot_clusters_graph` must be logical!")
})
