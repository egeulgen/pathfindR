##################################################
## Project: pathfindR
## Script purpose: Testthat testing script for
## enriched term clustering functions
## Date: Oct 29, 2019
## Author: Ege Ulgen
##################################################

# create_kappa_matrix -----------------------------------------------------
test_that("`create_kappa_matrix()` reates kappa matrix", {
  tmp <- RA_output[1:3, ]
  expect_is(create_kappa_matrix(tmp), "matrix")
  expect_is(create_kappa_matrix(tmp, use_description = TRUE), "matrix")

  tmp$Down_regulated[1] <- tmp$Up_regulated[1] <- ""
  expect_is(create_kappa_matrix(tmp), "matrix")

  tmp$non_Signif_Snw_Genes <- ""
  expect_is(create_kappa_matrix(tmp, use_active_snw_genes = TRUE), "matrix")
})

test_that("`create_kappa_matrix()` arg checks works", {
  expect_error(create_kappa_matrix(RA_output,
                                   use_description = "INVALID"),
               "`use_description` should be TRUE or FALSE")

  expect_error(create_kappa_matrix(RA_output,
                                   use_active_snw_genes = "INVALID"),
               "`use_active_snw_genes` should be TRUE or FALSE")

  expect_error(create_kappa_matrix("INVALID"),
               "`enrichment_res` should be a data frame of enrichment results")
  expect_error(create_kappa_matrix(RA_output[1, ]),
               "`enrichment_res` should contain at least 2 rows")

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
  valid_res <- RA_output[, -2]
  expect_silent(create_kappa_matrix(valid_res))
  invalid_res <- RA_output[, -1]
  expect_error(create_kappa_matrix(invalid_res),
               paste0("`enrichment_res` should contain all of ",
                      paste(dQuote(nec_cols), collapse = ", ")))
  # desc T
  nec_cols <- cr_cols(use_description = TRUE)
  valid_res <- RA_output[, -1]
  expect_silent(create_kappa_matrix(valid_res,
                                    use_description = TRUE))
  invalid_res <- RA_output[, -2]
  expect_error(create_kappa_matrix(invalid_res,
                                   use_description = TRUE),
               paste0("`enrichment_res` should contain all of ",
                      paste(dQuote(nec_cols), collapse = ", ")))
  # snw_g T
  nec_cols <- cr_cols(use_active_snw_genes = TRUE)
  valid_res <- RA_output
  valid_res$non_Signif_Snw_Genes <- ""
  expect_silent(create_kappa_matrix(valid_res,
                                    use_active_snw_genes = TRUE))
  expect_error(create_kappa_matrix(RA_output,
                                   use_active_snw_genes = TRUE),
               paste0("`enrichment_res` should contain all of ",
                      paste(dQuote(nec_cols), collapse = ", ")))
})

# hierarchical_term_clustering --------------------------------------------
test_that("H`hierarchical_term_clustering()` returns integer vector", {
  enrichment_res <- RA_output[1:5, ]
  kappa_mat <- create_kappa_matrix(enrichment_res)

  expect_is(hierarchical_term_clustering(kappa_mat, enrichment_res),
            "integer")
  expect_is(hierarchical_term_clustering(kappa_mat,
                                         enrichment_res,
                                         plot_hmap = TRUE),
            "integer")
})

test_that("Consistent `hierarchical_term_clustering()` results", {
  enrichment_res <- RA_output[1:5, ]
  kappa_mat1 <- create_kappa_matrix(enrichment_res, use_description = TRUE)
  kappa_mat2 <- create_kappa_matrix(enrichment_res, use_description = FALSE)

  expect_is(a <- hierarchical_term_clustering(kappa_mat1,
                                              enrichment_res,
                                              use_description = TRUE,
                                              plot_dend = FALSE),
            "integer")

  expect_is(b <- hierarchical_term_clustering(kappa_mat2,
                                              enrichment_res,
                                              use_description = FALSE,
                                              plot_dend = FALSE),
            "integer")

  expect_true(all(a == b))
  expect_true(all(names(a) != names(b)))
})

test_that("`hierarchical_term_clustering()` arg checks work", {
  expect_error(hierarchical_term_clustering(kappa_mat = list(),
                                            enrichment_res = data.frame()),
               "`kappa_mat` should be a symmetric matrix")
  expect_error(hierarchical_term_clustering(kappa_mat = matrix(nrow = 1,
                                                               ncol = 2),
                                            enrichment_res = data.frame()),
               "`kappa_mat` should be a symmetric matrix")

  mat <- matrix(nrow = 3, ncol = 3,
                dimnames = list(1:3, 1:3))
  expect_error(hierarchical_term_clustering(kappa_mat = mat,
                                            enrichment_res = data.frame(ID = 4:5)),
               "All terms in `kappa_mat` should be present in `enrichment_res`")

  expect_error(hierarchical_term_clustering(kappa_mat = mat,
                                            enrichment_res = data.frame(ID = 1:3),
                                            plot_hmap = "INVALID"),
               "`plot_hmap` should be logical")

  expect_error(hierarchical_term_clustering(kappa_mat = mat,
                                            enrichment_res = data.frame(ID = 1:3),
                                            plot_dend = "INVALID"),
               "`plot_dend` should be logical")
})

# fuzzy_term_clustering ---------------------------------------------------
test_that("`fuzzy_term_clustering()` returns matrix of cluster memberships", {
  kappa_mat <- create_kappa_matrix(RA_output)

  expect_is(fuzzy_term_clustering(kappa_mat, RA_output), "matrix")

  # lower kappa_threshold
  expect_is(fuzzy_term_clustering(kappa_mat,
                                  RA_output,
                                  kappa_threshold = -1), "matrix")
})

test_that("`fuzzy_term_clustering()` arg checks work", {
  expect_error(fuzzy_term_clustering(kappa_mat = list(),
                                     enrichment_res = data.frame()),
               "`kappa_mat` should be a symmetric matrix")
  expect_error(fuzzy_term_clustering(kappa_mat = matrix(nrow = 1,
                                                        ncol = 2),
                                     enrichment_res = data.frame()),
               "`kappa_mat` should be a symmetric matrix")

  mat <- matrix(nrow = 3, ncol = 3,
                dimnames = list(1:3, 1:3))
  expect_error(fuzzy_term_clustering(kappa_mat = mat,
                                     enrichment_res = data.frame(ID = 4:5)),
               "All terms in `kappa_mat` should be present in `enrichment_res`")

  expect_error(fuzzy_term_clustering(kappa_mat = mat,
                                     enrichment_res = data.frame(ID = 1:3),
                                     kappa_threshold = "INVALID"),
               "`kappa_threshold` should be numeric")
  expect_error(fuzzy_term_clustering(kappa_mat = mat,
                                     enrichment_res = data.frame(ID = 1:3),
                                     kappa_threshold = 1.5),
               "`kappa_threshold` should be at most 1 as kappa statistic is always <= 1")
})

test_that("Consistent fuzzy term clustering results", {
  kappa_mat1 <- create_kappa_matrix(RA_output, use_description = TRUE)
  kappa_mat2 <- create_kappa_matrix(RA_output, use_description = FALSE)

  expect_is(a <- fuzzy_term_clustering(kappa_mat1,
                                       RA_output,
                                       use_description = TRUE),
            "matrix")

  expect_is(b <- fuzzy_term_clustering(kappa_mat2,
                                       RA_output,
                                       use_description = FALSE),
            "matrix")

  expect_true(identical(unname(a), unname(b)))
  expect_true(all(rownames(a) != rownames(b)))
})

# cluster_graph_vis -------------------------------------------------------
test_that("Graph visualization of clusters via `cluster_graph_vis()` works OK", {

  ## use_description = FALSE
  enrichment_res <- RA_output[1:5, ]
  kappa_mat <- create_kappa_matrix(enrichment_res)

  # hierarchical
  clu_obj <- hierarchical_term_clustering(kappa_mat, enrichment_res)
  expect_silent(cluster_graph_vis(clu_obj, kappa_mat, enrichment_res))

  # fuzzy
  clu_obj <- fuzzy_term_clustering(kappa_mat, enrichment_res)
  expect_silent(cluster_graph_vis(clu_obj, kappa_mat, enrichment_res))

  ## use_description = TRUE
  kappa_mat <- create_kappa_matrix(enrichment_res, use_description = TRUE)

  # hierarchical
  clu_obj <- hierarchical_term_clustering(kappa_mat,
                                          enrichment_res,
                                          use_description = TRUE)
  expect_silent(cluster_graph_vis(clu_obj, kappa_mat,
                                  enrichment_res, use_description = TRUE))

  # fuzzy
  clu_obj <- fuzzy_term_clustering(kappa_mat, enrichment_res, use_description = TRUE)
  expect_silent(cluster_graph_vis(clu_obj, kappa_mat,
                                enrichment_res, use_description = TRUE))
})

test_that("`cluster_graph_vis()` coloring of 'extra' clusters work", {
  ### more than 41 clusters # better test case?
  N <- 45
  enrichment_res <- RA_output[1:N, ]
  kappa_mat <- create_kappa_matrix(enrichment_res)

  # simulated hierarchical result
  clu_obj <- 1:N
  names(clu_obj) <- RA_output$ID[1:N]
  expect_silent(cluster_graph_vis(clu_obj, kappa_mat, enrichment_res))

  # simulated fuzzy result
  clu_obj <- matrix(FALSE,
                    nrow = N, ncol = N,
                    dimnames = list(RA_output$ID[1:N], 1:N))
  diag(clu_obj) <- TRUE
  expect_silent(cluster_graph_vis(clu_obj, kappa_mat, enrichment_res))
})

test_that("Check errors of `cluster_graph_vis()`", {
  expect_error(cluster_graph_vis(list(), matrix(), data.frame(ID = 1)),
               "Invalid class for `clu_obj`!")

  enrichment_res <- RA_output[1:5, ]
  kappa_mat <- create_kappa_matrix(enrichment_res)

  # hierarchical - missing terms in kappa matrix
  clu_obj <- hierarchical_term_clustering(kappa_mat, enrichment_res)
  expect_error(cluster_graph_vis(c(clu_obj, EXTRA = 1L),
                                 kappa_mat, enrichment_res),
               "Not all terms in `clu_obj` present in `kappa_mat`!")

  # fuzzy - missing terms in kappa matrix
  clu_obj <- fuzzy_term_clustering(kappa_mat, enrichment_res)
  expect_error(cluster_graph_vis(rbind(clu_obj,
                                       EXTRA = rep(FALSE, ncol(clu_obj))),
                                 kappa_mat, enrichment_res),
               "Not all terms in `clu_obj` present in `kappa_mat`!")
})

# cluster_enriched_terms --------------------------------------------------
test_that("`cluster_enriched_terms()` returns the input data frame
          with the additional columns `Cluster` and `Status`", {

  #### use_description = FALSE
  ## hierarchical
  tmp <- cluster_enriched_terms(RA_output[1:5, ])
  expect_is(tmp, "data.frame")
  expect_true(all(c("Cluster", "Status") %in% colnames(tmp)))
  # expect to have same number of rep. terms as the number of clusters
  expect_true(max(tmp$Cluster) == sum(tmp$Status == "Representative"))

  ## fuzzy
  tmp <- cluster_enriched_terms(RA_output[1:5, ], method = "fuzzy")
  expect_is(tmp, "data.frame")
  expect_true(all(c("Cluster", "Status") %in% colnames(tmp)))
  # expect to have same number of rep. terms as the number of clusters
  expect_true(max(tmp$Cluster) == sum(tmp$Status == "Representative"))

  #### use_description = TRUE
  ## hierarchical
  tmp <- cluster_enriched_terms(RA_output[1:5, ], use_description = TRUE)
  expect_is(tmp, "data.frame")
  expect_true(all(c("Cluster", "Status") %in% colnames(tmp)))
  # expect to have same number of rep. terms as the number of clusters
  expect_true(max(tmp$Cluster) == sum(tmp$Status == "Representative"))

  ## fuzzy
  tmp <- cluster_enriched_terms(RA_output[1:5, ],
                                method = "fuzzy",
                                use_description = TRUE)
  expect_is(tmp, "data.frame")
  expect_true(all(c("Cluster", "Status") %in% colnames(tmp)))
  # expect to have same number of rep. terms as the number of clusters
  expect_true(max(tmp$Cluster) == sum(tmp$Status == "Representative"))
})

test_that("`cluster_enriched_terms()` arg checks work", {
  expect_error(cluster_enriched_terms(RA_output[1, 3, ], method = "WRONG"),
    "the clustering `method` must either be \"hierarchical\" or \"fuzzy\"")

  expect_error(cluster_enriched_terms(RA_output[1:3, ],
                                      plot_clusters_graph = "WRONG"),
               "`plot_clusters_graph` must be logical!")
})
