##################################################
## Project: pathfindR
## Script purpose: Testthat testing script for
## visualization-related functions
## Date: Nov 10, 2019
## Author: Ege Ulgen
##################################################

# visualize_terms ---------------------------------------------------------
tmp_res <- RA_output[1, ]

tmp_genes <- unlist(c(strsplit(tmp_res$Up_regulated, ", "),
                      strsplit(tmp_res$Down_regulated, ", ")))
input_processed <- suppressMessages(input_processing(RA_input, 0.05, "Biogrid"))

test_that("`visualize_terms()` creates expected png file(s)", {
  ## non-KEGG (visualize_term_interactions)
  expected_out_file <- file.path("term_visualizations",
                                 paste0(tmp_res$Term_Description, ".png"))
  suppressMessages(visualize_terms(result_df = tmp_res,
                                   input_processed = input_processed,
                                   hsa_KEGG = FALSE,
                                   pin_name_path = "Biogrid")) # default
  expect_true(file.exists(expected_out_file))
  unlink("term_visualizations", recursive = TRUE)


  suppressMessages(visualize_terms(result_df = tmp_res,
                                   input_processed = input_processed,
                                   hsa_KEGG = FALSE,
                                   pin_name_path = "GeneMania"))
  expect_true(file.exists(expected_out_file))
  unlink("term_visualizations", recursive = TRUE)

  skip_on_cran()
  ## hsa KEGG
  expected_out_file <- file.path("term_visualizations",
                                 paste0(tmp_res$ID, "_pathfindR.png"))

  suppressMessages(visualize_terms(result_df = tmp_res,
                                   input_processed = input_processed,
                                   hsa_KEGG = TRUE))
  expect_true(file.exists(expected_out_file))
  unlink("term_visualizations", recursive = TRUE)
})


test_that("`visualize_terms()` arg checks work", {
  expect_error(visualize_terms(result_df = "INVALID"),
               "`result_df` should be a data frame")

  # hsa_KEGG = TRUE
  nec_cols <- "ID"
  expect_error(visualize_terms(tmp_res[, -1],
                               hsa_KEGG = TRUE),
               paste0("`result_df` should contain the following columns: ",
                      paste(dQuote(nec_cols), collapse = ", ")))
  # hsa_KEGG = FALSE
  nec_cols <- c("Term_Description", "Up_regulated", "Down_regulated")
  expect_error(visualize_terms(tmp_res[, -2],
                               hsa_KEGG = FALSE),
               paste0("`result_df` should contain the following columns: ",
                      paste(dQuote(nec_cols), collapse = ", ")))

  expect_error(visualize_terms(result_df = tmp_res,
                               hsa_KEGG = TRUE),
               "`input_processed` should be specified when `hsa_KEGG = TRUE`")

  expect_error(visualize_terms(result_df = tmp_res,
                               hsa_KEGG = "INVALID"),
               "the argument `hsa_KEGG` should be either TRUE or FALSE")
})

# visualize_term_interactions ---------------------------------------------
test_that("`visualize_term_interactions()` creates expected png file(s)", {
  expected_out_file <- file.path("term_visualizations",
                                 paste0(tmp_res$Term_Description, ".png"))
  expect_null(visualize_term_interactions(tmp_res,
                                          pin_name_path = "Biogrid"))
  expect_true(file.exists(expected_out_file))
  unlink("term_visualizations", recursive = TRUE)

  tmp_res2 <- rbind(tmp_res, tmp_res)
  tmp_res2$Term_Description[2] <- "SKIP"
  tmp_res2$Up_regulated[2] <- ""
  tmp_res2$Down_regulated[2] <- ""
  expect_message(visualize_term_interactions(tmp_res2,
                                             pin_name_path = "KEGG"),
                 paste0("< 2 genes, skipping visualization of ",
                        tmp_res2$Term_Description[2]))
  unlink("term_visualizations", recursive = TRUE)

  # Non-empty non_Signif_Snw_Genes
  tmp_res3 <- tmp_res
  tmp_res3$non_Signif_Snw_Genes <- RA_output$Up_regulated[2]
  expect_null(visualize_term_interactions(tmp_res3,
                                          pin_name_path = "Biogrid"))
  expect_true(file.exists(expected_out_file))
  unlink("term_visualizations", recursive = TRUE)
})

# visualize_hsa_KEGG ------------------------------------------------------
test_that("`visualize_hsa_KEGG()` creates expected png file(s)", {
  skip_on_cran()

  expected_out_file <- file.path("term_visualizations",
                                 paste0(tmp_res$ID, "_pathfindR.png"))

  ###### Continuous change values
  ### Normalize = TRUE
  expect_null(visualize_hsa_KEGG(hsa_kegg_ids = tmp_res$ID,
                                 input_processed = input_processed))
  expect_true(file.exists(expected_out_file))
  unlink("term_visualizations", recursive = TRUE)

  expect_null(visualize_hsa_KEGG(hsa_kegg_ids = tmp_res$ID,
                                 input_processed = input_processed,
                                 node_cols = c("red", "green", "blue")))
  expect_true(file.exists(expected_out_file))
  unlink("term_visualizations", recursive = TRUE)

  ### Normalize = FALSE
  expect_null(visualize_hsa_KEGG(hsa_kegg_ids = tmp_res$ID,
                                 input_processed = input_processed,
                                 normalize_vals = FALSE))
  expect_true(file.exists(expected_out_file))
  unlink("term_visualizations", recursive = TRUE)

  ###### Constant change values (if no change values supplied in input)
  input_processed2 <- input_processed
  input_processed2$CHANGE <- 1e6
  expect_null(visualize_hsa_KEGG(hsa_kegg_ids = tmp_res$ID,
                                 input_processed = input_processed2))
  expect_true(file.exists(expected_out_file))
  unlink("term_visualizations", recursive = TRUE)

  expect_message(visualize_hsa_KEGG(hsa_kegg_ids = tmp_res$ID,
                                    input_processed = input_processed2,
                                    node_cols = "red"),
                 "all `change_vec` values are 1e6, using the first color in `node_cols`")
  expect_true(file.exists(expected_out_file))
  unlink("term_visualizations", recursive = TRUE)


  ###### max_to_plot works
  max_n <- 5
  expected_out_files <- file.path("term_visualizations",
                                  paste0(RA_output$ID[seq_len(max_n)],
                                         "_pathfindR.png"))
  expect_null(visualize_hsa_KEGG(hsa_kegg_ids = RA_output$ID,
                                 input_processed = input_processed,
                                 max_to_plot = max_n))
  expect_true(all(file.exists(expected_out_files)))
  unlink("term_visualizations", recursive = TRUE)

})

temp_ids <- RA_output$ID[1:3]
test_that("KEGML download error is handled properly", {
  skip_on_cran()
  temp_ids[2] <- "hsa00000"
  expected_out_files <- file.path("term_visualizations",
                                  paste0(temp_ids, "_pathfindR.png"))
  visualize_hsa_KEGG(hsa_kegg_ids = temp_ids,
                     input_processed = input_processed)

  expect_equal(file.exists(expected_out_files), c(TRUE, FALSE, TRUE))
  unlink("term_visualizations", recursive = TRUE)
})

test_that("arg checks for `visualize_hsa_KEGG()` work", {
  expect_error(visualize_hsa_KEGG(hsa_kegg_ids = list(),
                                  input_processed = input_processed),
               "`hsa_kegg_ids` should be a vector of hsa KEGG IDs")
  expect_error(visualize_hsa_KEGG(hsa_kegg_ids = c("X", "Y", "Z"),
                                  input_processed = input_processed),
               "`hsa_kegg_ids` should be a vector of valid hsa KEGG IDs")

  expect_error(visualize_hsa_KEGG(hsa_kegg_ids = tmp_res$ID,
                                  input_processed = list()),
               "`input_processed` should be a data frame")
  expect_error(visualize_hsa_KEGG(hsa_kegg_ids = tmp_res$ID,
                                  input_processed = input_processed[, -2]),
               paste0("`input_processed` should contain the following columns: ",
                      paste(dQuote(c("GENE", "CHANGE")), collapse = ", ")))

  expect_error(visualize_hsa_KEGG(hsa_kegg_ids = tmp_res$ID,
                                  input_processed = input_processed,
                                  max_to_plot = "INVALID"),
               "`max_to_plot` should be numeric or NULL")
  expect_error(visualize_hsa_KEGG(hsa_kegg_ids = tmp_res$ID,
                                  input_processed = input_processed,
                                  max_to_plot = 0),
               "`max_to_plot` should be >=1")




  expect_error(visualize_hsa_KEGG(hsa_kegg_ids = tmp_res$ID,
                                  input_processed = input_processed,
                                  normalize_vals = "INVALID"),
               "`normalize_vals` should be logical")

  expect_error(visualize_hsa_KEGG(hsa_kegg_ids = tmp_res$ID,
                                  input_processed = input_processed,
                                  node_cols = list()),
               "`node_cols` should be a vector of colors")
  expect_error(visualize_hsa_KEGG(hsa_kegg_ids = tmp_res$ID,
                                  input_processed = input_processed,
                                  node_cols = rep("red", 4)),
               "the length of `node_cols` should be 3")
  expect_error(visualize_hsa_KEGG(hsa_kegg_ids = tmp_res$ID,
                                  input_processed = input_processed,
                                  node_cols = c("red", "#FFFFFF", "INVALID")),
               "`node_cols` should be a vector of valid colors")
})

# enrichment_chart --------------------------------------------------------
test_that("enrichment_chart produces a ggplot object with correct labels", {

  # default
  expect_is(g <- enrichment_chart(RA_output), "ggplot")
  expect_equal(ggplot2::quo_name(g$mapping$x), "Fold_Enrichment")
  expect_equal(ggplot2::quo_name(g$mapping$y), "Term_Description")
  expect_equal(g$labels$size, "# genes")
  expect_equal(g$labels$colour, expression(-log[10](p)))
  expect_equal(g$labels$x, "Fold Enrichment")
  expect_equal(g$labels$y, "Term_Description")

  # plot_by_cluster
  expect_is(g <- enrichment_chart(RA_clustered, plot_by_cluster = TRUE),
            "ggplot")
  expect_equal(ggplot2::quo_name(g$mapping$x), "Fold_Enrichment")
  expect_equal(ggplot2::quo_name(g$mapping$y), "Term_Description")
  expect_equal(g$labels$size, "# genes")
  expect_equal(g$labels$colour, expression(-log[10](p)))
  expect_equal(g$labels$x, "Fold Enrichment")
  expect_equal(g$labels$y, "Term_Description")

  # chang top_terms
  expect_is(g <- enrichment_chart(RA_output, top_terms = NULL),
            "ggplot")
  expect_equal(ggplot2::quo_name(g$mapping$x), "Fold_Enrichment")
  expect_equal(ggplot2::quo_name(g$mapping$y), "Term_Description")
  expect_equal(g$labels$size, "# genes")
  expect_equal(g$labels$colour, expression(-log[10](p)))
  expect_equal(g$labels$x, "Fold Enrichment")
  expect_equal(g$labels$y, "Term_Description")

  expect_is(g <- enrichment_chart(RA_output, top_terms = 1e3),
            "ggplot")
  expect_equal(ggplot2::quo_name(g$mapping$x), "Fold_Enrichment")
  expect_equal(ggplot2::quo_name(g$mapping$y), "Term_Description")
  expect_equal(g$labels$size, "# genes")
  expect_equal(g$labels$colour, expression(-log[10](p)))
  expect_equal(g$labels$x, "Fold Enrichment")
  expect_equal(g$labels$y, "Term_Description")

  # change num_bubbles
  expect_is(g <- enrichment_chart(RA_clustered,
                                  num_bubbles = 30),
            "ggplot")
  expect_equal(ggplot2::quo_name(g$mapping$x), "Fold_Enrichment")
  expect_equal(ggplot2::quo_name(g$mapping$y), "Term_Description")
  expect_equal(g$labels$size, "# genes")
  expect_equal(g$labels$colour, expression(-log[10](p)))
  expect_equal(g$labels$x, "Fold Enrichment")
  expect_equal(g$labels$y, "Term_Description")

  # change even_breaks
  expect_is(g <- enrichment_chart(RA_clustered, even_breaks = FALSE),
            "ggplot")
  expect_equal(ggplot2::quo_name(g$mapping$x), "Fold_Enrichment")
  expect_equal(ggplot2::quo_name(g$mapping$y), "Term_Description")
  expect_equal(g$labels$size, "# genes")
  expect_equal(g$labels$colour, expression(-log[10](p)))
  expect_equal(g$labels$x, "Fold Enrichment")
  expect_equal(g$labels$y, "Term_Description")
})

test_that("enrichment_chart arg checks work", {
  necessary <- c("Term_Description", "Fold_Enrichment", "lowest_p",
                 "Up_regulated", "Down_regulated")
  expect_error(enrichment_chart(RA_output[, -2]),
               paste0("The input data frame must have the columns:\n",
                      paste(necessary, collapse = ", ")))

  expect_error(enrichment_chart(RA_output, plot_by_cluster = "INVALID"),
    "`plot_by_cluster` must be either TRUE or FALSE")

  expect_message(enrichment_chart(RA_output, plot_by_cluster = TRUE),
    "For plotting by cluster, there must a column named `Cluster` in the input data frame!")

  expect_error(enrichment_chart(RA_output, top_terms = "INVALID"),
               "`top_terms` must be either numeric or NULL")

  expect_error(enrichment_chart(RA_output, top_terms = 0),
               "`top_terms` must be > 1")
})

# term_gene_graph ---------------------------------------------------------
test_that("`term_gene_graph()` produces a ggplot object using the correct data", {

  # Top 10 (default)
  expect_is(p <- term_gene_graph(RA_output), "ggplot")
  expect_equal(sum(p$data$type == "term"), 10)

  # Top 3
  expect_is(p <- term_gene_graph(RA_output, num_terms = 3), "ggplot")
  expect_equal(sum(p$data$type == "term"), 3)

  # All terms
  expect_is(p <- term_gene_graph(RA_output, num_terms = NULL), "ggplot")
  expect_equal(sum(p$data$type == "term"), nrow(RA_output))

  # Top 1000, expect to plot top nrow(output)
  expect_is(p <- term_gene_graph(RA_output, num_terms = 1e3), "ggplot")
  expect_equal(sum(p$data$type == "term"), nrow(RA_output))

  # use_description = TRUE
  expect_is(p2 <- term_gene_graph(RA_output, use_description = TRUE), "ggplot")
  expect_equal(sum(p2$data$type == "term"), 10)

  # use_names = "p_val"
  expect_is(p <- term_gene_graph(RA_output, node_size = "p_val"), "ggplot")
  expect_equal(sum(p$data$type == "term"), 10)
})

test_that("`term_gene_graph()` arg checks work", {
  expect_error(term_gene_graph(RA_output, num_terms = "INVALID"),
    "`num_terms` must either be numeric or NULL!")

  expect_error(term_gene_graph(RA_output, use_description = "INVALID"),
    "`use_description` must either be TRUE or FALSE!")

  val_node_size <- c("num_genes", "p_val")
  expect_error(term_gene_graph(RA_output, node_size = "INVALID"),
               paste0("`node_size` should be one of ",
                      paste(dQuote(val_node_size), collapse = ", ")))

  wrong_df <- RA_output[, -c(1, 2)]
  ID_column <- "ID"
  necessary_cols <- c("Up_regulated", "Down_regulated", "lowest_p", ID_column)
  expect_error(term_gene_graph(wrong_df, use_description = FALSE),
               paste(c("All of", paste(necessary_cols, collapse = ", "),
                       "must be present in `results_df`!"), collapse = " "))

  ID_column <- "Term_Description"
  necessary_cols <- c("Up_regulated", "Down_regulated", "lowest_p", ID_column)
  expect_error(term_gene_graph(wrong_df, use_description = TRUE),
               paste(c("All of", paste(necessary_cols, collapse = ", "),
                       "must be present in `results_df`!"), collapse = " "))
})
