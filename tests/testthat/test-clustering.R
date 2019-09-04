##################################################
## Project: pathfindR
## Script purpose: Testthat testing script for
## pathway clustering functions
## Date: Aug 28, 2019
## Author: Ege Ulgen
##################################################

# create_kappa_matrix -----------------------------------------------------
test_that("Creates kappa matrix", {
  tmp <- RA_output[1:3, ]
  expect_is(create_kappa_matrix(tmp), "matrix")
  expect_is(create_kappa_matrix(tmp, use_names = TRUE), "matrix")

  tmp$Down_regulated[1] <- tmp$Up_regulated[1] <- ""
  expect_is(create_kappa_matrix(tmp), "matrix")

  tmp$non_DEG_Active_Snw_Genes <- ""
  expect_is(create_kappa_matrix(tmp, use_active_snw_genes = TRUE), "matrix")
})

test_that("kappa matrix function error checks", {
  expect_error(create_kappa_matrix(RA_output[1:3, ],
                                   use_active_snw_genes = TRUE),
               "No column named `non_DEG_Active_Snw_Genes`, please execute `run_pathfindR` with `list_active_snw_genes = TRUE`!")

  expect_error(create_kappa_matrix(RA_output[1:3, ],
                                   use_names = "WRONG"),
               "`use_names` must be logical!")
})

# hierarchical_pw_clustering ----------------------------------------------
test_that("Hierarchical pw clustering returns integer vector", {
  enrichment_res <- RA_output[1:5, ]
  kappa_mat <- create_kappa_matrix(enrichment_res)

  expect_is(hierarchical_pw_clustering(kappa_mat, enrichment_res), "integer")
  expect_is(hierarchical_pw_clustering(kappa_mat,
                                       enrichment_res,
                                       plot_hmap = TRUE),
            "integer")

  expect_is(hierarchical_pw_clustering(kappa_mat,
                                       enrichment_res,
                                       plot_dend = TRUE),
            "integer")
})


# fuzzy_pw_clustering -----------------------------------------------------
test_that("Fuzzy pw clustering returns matrix of cluster memberships", {
  kappa_mat <- create_kappa_matrix(RA_output)

  expect_is(fuzzy_pw_clustering(kappa_mat, RA_output), "matrix")

  # lower kappa_threshold
  expect_is(fuzzy_pw_clustering(kappa_mat,
                                RA_output,
                                kappa_threshold = -1), "matrix")
})

test_that("Error raised if kappa threshold is not numeric", {
  expect_error(fuzzy_pw_clustering(kappa_mat, enrichment_res,
                                   kappa_threshold = "test"),
               "`kappa_threshold` must be numeric!")
})


# cluster_graph_vis -------------------------------------------------------
test_that("Graph visualization of clusters works OK", {

  ## use_names = FALSE
  enrichment_res <- RA_output[1:5, ]
  kappa_mat <- create_kappa_matrix(enrichment_res)

  # hierarchical
  clu_obj <- hierarchical_pw_clustering(kappa_mat, enrichment_res)
  expect_identical(cluster_graph_vis(clu_obj, kappa_mat, enrichment_res), NULL)

  # fuzzy
  clu_obj <- fuzzy_pw_clustering(kappa_mat, enrichment_res)
  expect_identical(cluster_graph_vis(clu_obj, kappa_mat, enrichment_res), NULL)

  ## use_names = TRUE
  enrichment_res <- RA_output[1:5, ]
  kappa_mat <- create_kappa_matrix(enrichment_res, use_names = TRUE)

  # hierarchical
  clu_obj <- hierarchical_pw_clustering(kappa_mat,
                                        enrichment_res, use_names = TRUE)
  expect_identical(cluster_graph_vis(clu_obj, kappa_mat, enrichment_res,
                                     use_names = TRUE), NULL)

  # fuzzy
  clu_obj <- fuzzy_pw_clustering(kappa_mat, enrichment_res, use_names = TRUE)
  expect_identical(cluster_graph_vis(clu_obj, kappa_mat, enrichment_res,
                                     use_names = TRUE), NULL)

  ### more than 41 clusters # better test case?
  N <- 45
  enrichment_res <- RA_output[1:N, ]
  kappa_mat <- create_kappa_matrix(enrichment_res)

  # hierarchical
  clu_obj <- 1:N
  names(clu_obj) <- RA_output$ID[1:N]
  expect_silent(cluster_graph_vis(clu_obj, kappa_mat, enrichment_res))

  # fuzzy
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

  expect_error(cluster_graph_vis(c(1L), matrix(0, nrow = 2, ncol= 1),
                                 data.frame()),
               "`kappa_mat` must be symmetric!")

  expect_error(cluster_graph_vis(c(1L), matrix(), list()),
               "`enrichment_res` must be a data.frame!")

  expect_error(cluster_graph_vis(c(1L), matrix(), data.frame(a = 1:3)),
            "`kappa_mat` and `enrichment_res` must contain the same # of terms")

  expect_error(cluster_graph_vis(c(1L),
                                 matrix(nrow = 3, ncol = 3,
                                        dimnames = list(c("A", "B", "D"),
                                                        c("A", "B", "D"))),
                                 data.frame(ID = c("A", "B", "C"))),
               "Not all terms in `kappa_mat` and `enrichment_res` match!")


  enrichment_res <- RA_output[1:5, ]
  kappa_mat <- create_kappa_matrix(enrichment_res)

  # hierarchical - missing terms in kappa matrix
  clu_obj <- hierarchical_pw_clustering(kappa_mat, enrichment_res)
  expect_error(cluster_graph_vis(c(clu_obj, EXTRA = 1L),
                                 kappa_mat, enrichment_res),
               "Not all terms in `clu_obj` present in `kappa_mat`!")

  # fuzzy - missing terms in kappa matrix
  clu_obj <- fuzzy_pw_clustering(kappa_mat, enrichment_res)
  expect_error(cluster_graph_vis(rbind(clu_obj,
                                       EXTRA = rep(FALSE, ncol(clu_obj))),
                                 kappa_mat, enrichment_res),
               "Not all terms in `clu_obj` present in `kappa_mat`!")
})

# cluster_pathways --------------------------------------------------------
test_that("Clustering wrapper returns the input data frame
          with the additional columns `Cluster` and `Status`", {

            ## hierarchical
            tmp <- cluster_pathways(RA_output[1:3, ])
            expect_is(tmp, "data.frame")

            expect_equal(c("Cluster", "Status") %in% colnames(tmp),
                         c(TRUE, TRUE))

            ## fuzzy
            tmp <- cluster_pathways(RA_output[1:3, ], method = "fuzzy")
            expect_is(tmp, "data.frame")

            expect_equal(c("Cluster", "Status") %in% colnames(tmp),
                         c(TRUE, TRUE))

            ## use_names = TRUE
            ## hierarchical
            tmp <- cluster_pathways(RA_output[1:3, ], use_names = TRUE)
            expect_is(tmp, "data.frame")

            expect_equal(c("Cluster", "Status") %in% colnames(tmp),
                         c(TRUE, TRUE))

            ## fuzzy
            tmp <- cluster_pathways(RA_output[1:3, ],
                                    method = "fuzzy", use_names = TRUE)
            expect_is(tmp, "data.frame")

            expect_equal(c("Cluster", "Status") %in% colnames(tmp),
                         c(TRUE, TRUE))
          })

test_that("Check errors of clustering wrapper function", {
  expect_error(cluster_pathways(RA_output[1:3, ], method = "WRONG"),
        "the clustering `method` must either be \"hierarchical\" or \"fuzzy\"")

  expect_error(cluster_pathways(RA_output[1:3, ], use_names = "WRONG"),
               "`use_names` must be logical!")

  expect_error(cluster_pathways(RA_output[1:3, ],
                                plot_clusters_graph = "WRONG"),
               "`plot_clusters_graph` must be logical!")
})
