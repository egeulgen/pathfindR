##################################################
## Project: pathfindR
## Script purpose: Testthat testing script for
## enriched term clustering functions
## Date: Oct 20, 2019
## Author: Ege Ulgen
##################################################

# create_kappa_matrix -----------------------------------------------------
test_that("Creates kappa matrix", {
  tmp <- RA_output[1:3, ]
  expect_is(create_kappa_matrix(tmp), "matrix")
  expect_is(create_kappa_matrix(tmp, use_description = TRUE), "matrix")

  tmp$Down_regulated[1] <- tmp$Up_regulated[1] <- ""
  expect_is(create_kappa_matrix(tmp), "matrix")

  tmp$non_DEG_Active_Snw_Genes <- ""
  expect_is(create_kappa_matrix(tmp, use_active_snw_genes = TRUE), "matrix")
})

test_that("kappa matrix creator function error checks works", {
  expect_error(create_kappa_matrix(RA_output[1:3, ],
                                   use_active_snw_genes = TRUE),
  "No column named `non_DEG_Active_Snw_Genes`,
      please execute `run_pathfindR` with `list_active_snw_genes = TRUE`!")

  expect_error(create_kappa_matrix(RA_output[1:3, ],
                                   use_description = "WRONG"),
               "`use_description` must be logical!")
})

# hierarchical_term_clustering --------------------------------------------
test_that("Hierarchical term clustering returns integer vector", {
  enrichment_res <- RA_output[1:5, ]
  kappa_mat <- create_kappa_matrix(enrichment_res)

  expect_is(hierarchical_term_clustering(kappa_mat, enrichment_res),
            "integer")
  expect_is(hierarchical_term_clustering(kappa_mat,
                                         enrichment_res,
                                         plot_hmap = TRUE),
            "integer")
})

test_that("Consistent hierarchical term clustering results", {
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

# fuzzy_term_clustering ---------------------------------------------------
test_that("Fuzzy term clustering returns matrix of cluster memberships", {
  kappa_mat <- create_kappa_matrix(RA_output)

  expect_is(fuzzy_term_clustering(kappa_mat, RA_output), "matrix")

  # lower kappa_threshold
  expect_is(fuzzy_term_clustering(kappa_mat,
                                  RA_output,
                                  kappa_threshold = -1), "matrix")
})

test_that("Error raised if kappa threshold is not numeric", {
  expect_error(fuzzy_term_clustering(kappa_mat, enrichment_res,
                                     kappa_threshold = "test"),
               "`kappa_threshold` must be numeric!")
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
test_that("Graph visualization of clusters works OK", {

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

test_that("Test that coloring of 'extra' clusters work", {
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


test_that("Check errors of clustering visualizations", {
  expect_error(cluster_graph_vis(list(), matrix(), data.frame(ID = 1)),
               "Invalid class for `clu_obj`!")

  expect_error(cluster_graph_vis(c(1L), list(), data.frame()),
               "`kappa_mat` must be a matrix!")

  expect_error(cluster_graph_vis(c(1L), matrix(0, nrow = 2, ncol = 1),
                                 data.frame()),
    "`kappa_mat` must be a symmetric matrix!")

  expect_error(cluster_graph_vis(c(1L), matrix(), list()),
               "`enrichment_res` must be a data.frame!")

  expect_error(cluster_graph_vis(c(1L), matrix(), data.frame(a = 1:3)),
    "`kappa_mat` and `enrichment_res` must contain the same # of terms")

  expect_error(cluster_graph_vis(c(1L),
                                 matrix(nrow = 3, ncol = 3,
                                        dimnames = list(c("A", "B", "D"),
                                                        c("A", "B", "D"))),
                                 data.frame(ID = c("A", "B", "C"))),
               "Not all terms in `kappa_mat` and `enrichment_res` match!"
  )

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
test_that("Clustering wrapper returns the input data frame
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

test_that("Check errors of clustering wrapper function", {
  expect_error(cluster_enriched_terms(RA_output[1, 3, ], method = "WRONG"),
    "the clustering `method` must either be \"hierarchical\" or \"fuzzy\"")

  expect_error(cluster_enriched_terms(RA_output[1:3, ],
                                      use_description = "WRONG"),
               "`use_description` must be logical!")

  expect_error(cluster_enriched_terms(RA_output[1:3, ],
                                      plot_clusters_graph = "WRONG"),
               "`plot_clusters_graph` must be logical!")
})
