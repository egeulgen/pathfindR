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
  enrichment_res <- RA_output[1:5, ]
  kappa_mat <- create_kappa_matrix(enrichment_res)

  expect_is(fuzzy_pw_clustering(kappa_mat, enrichment_res), "matrix")

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

  # more than 41 clusters # better test case?
  enrichment_res <- RA_output[1:42, ]
  kappa_mat <- create_kappa_matrix(enrichment_res)
  clu_obj <- 1:42
  names(clu_obj) <- RA_output$ID[1:42]
  expect_identical(cluster_graph_vis(clu_obj, kappa_mat, enrichment_res), NULL)
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
test_that("Clustering returns appropriate object", {
  tmp <- RA_output[1:3, ]
  expect_is(create_kappa_matrix(tmp), "matrix")
  expect_is(create_kappa_matrix(tmp, use_names = TRUE), "matrix")

  tmp$non_DEG_Active_Snw_Genes <- ""
  expect_is(create_kappa_matrix(tmp, use_active_snw_genes = TRUE), "matrix")
})
