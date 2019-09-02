##################################################
## Project: pathfindR
## Script purpose: Testthat testing script for
## visualization-related functions
## Date: Aug 30, 2019
## Author: Ege Ulgen
##################################################

# visualize_pws -----------------------------------------------------------
processed_input <- suppressMessages(input_processing(RA_input[1:10, ],
                                                     0.05, return_pin_path()))
tmp_res <- RA_output[3, ]

test_that("visualize_pws creates expected png file(s)", {

  ## hsa KEGG (pathview)
  expected_out_file <- file.path(paste0("pathway_visualizations/",
                                        paste(tmp_res$ID,
                                              tmp_res$Pathway, "png",
                                              sep = ".")))
  suppressMessages(visualize_pws(result_df = tmp_res,
                                 input_processed = processed_input,
                                 gene_sets = "KEGG"))
  expect_true(file.exists(expected_out_file))

  ## non-KEGG (visualize_pw_interactions)
  expected_out_file <- file.path(paste0("pathway_visualizations/",
                                        paste(tmp_res$Pathway, "png",
                                              sep = ".")))
  suppressMessages(visualize_pws(result_df = tmp_res,
                                 input_processed = processed_input,
                                 gene_sets = "non-KEGG",
                                 pin_name_path = "Biogrid")) # default
  expect_true(file.exists(expected_out_file))

  suppressMessages(visualize_pws(result_df = tmp_res,
                                 input_processed = processed_input,
                                 gene_sets = "non-KEGG",
                                 pin_name_path = "GeneMania")) # default
  expect_true(file.exists(expected_out_file))
})

test_that("visualize_pws func. arg check works", {
  expect_error(visualize_pws(result_df = tmp_res),
               "`input_processed` must be specified when `gene_sets` is KEGG")
  expect_error(visualize_pws(result_df = tmp_res,
                             gene_sets = "WRONG"),
               "`gene_sets` must either be \"KEGG\" or \"non-KEGG\"")
})


# visualize_pw_interactions -----------------------------------------------
test_that("visualize_pw_interactions creates expected png file(s)", {
  expected_out_file <- file.path(paste0("pathway_visualizations/",
                                        paste(tmp_res$Pathway, "png",
                                              sep = ".")))
  expect_null(visualize_pw_interactions(tmp_res,
                                        pin_name_path = "Biogrid"))
  expect_true(file.exists(expected_out_file))

  tmp_res2 <- rbind(tmp_res, tmp_res)
  tmp_res2$Pathway[2] <- "SKIP"
  tmp_res2$Up_regulated[2] <- ""
  tmp_res2$Down_regulated[2] <- ""
  expect_message(visualize_pw_interactions(tmp_res2,
                                           pin_name_path = "KEGG"),
                 paste0("< 2 genes, skipping visualization of ",
                        tmp_res2$Pathway[2]))
})

# visualize_hsa_KEGG ------------------------------------------------------
test_that("visualize_hsa_KEGG creates expected png file(s)", {
  expected_out_file <- file.path(paste0("pathway_visualizations/",
                                        paste(tmp_res$ID,
                                              tmp_res$Pathway, "png",
                                              sep = ".")))

  ## Continuous change values
  genes_df <- processed_input[, c("GENE", "CHANGE")]
  rownames(genes_df) <- genes_df$GENE
  genes_df <- genes_df[, -1, drop = FALSE]

  expect_null(suppressMessages(visualize_hsa_KEGG(tmp_res, genes_df)))
  expect_true(file.exists(expected_out_file))

  ## Binary change values
  genes_df$CHANGE <- ifelse(genes_df$CHANGE < 0, -1, 1)

  expect_null(suppressMessages(visualize_hsa_KEGG(tmp_res, genes_df)))
  expect_true(file.exists(expected_out_file))

  ## Constant change values (if no change values supplied in input)
  genes_df$CHANGE <- 100

  expect_null(suppressMessages(visualize_hsa_KEGG(tmp_res, genes_df)))
  expect_true(file.exists(expected_out_file))
})

# enrichment_chart --------------------------------------------------------
test_that("enrichment_chart produces a ggplot object with correct labels", {

  # default
  expect_is(g <- enrichment_chart(RA_output), "ggplot")
  expect_equal(ggplot2::quo_name(g$mapping$x), "Fold_Enrichment")
  expect_equal(ggplot2::quo_name(g$mapping$y), "Pathway")
  expect_equal(g$labels$size, "# of DEGs")
  expect_equal(g$labels$colour, "-log10(lowest-p)")
  expect_equal(g$labels$x, "Fold Enrichment")
  expect_equal(g$labels$y, "")

  # plot_by_cluster
  expect_is(g <- enrichment_chart(RA_clustered,
                                  plot_by_cluster = TRUE),
            "ggplot")
  expect_equal(ggplot2::quo_name(g$mapping$x), "Fold_Enrichment")
  expect_equal(ggplot2::quo_name(g$mapping$y), "Pathway")
  expect_equal(g$labels$size, "# of DEGs")
  expect_equal(g$labels$colour, "-log10(lowest-p)")
  expect_equal(g$labels$x, "Fold Enrichment")
  expect_equal(g$labels$y, "")

  # change num_bubbles
  expect_is(g <- enrichment_chart(RA_clustered,
                                  num_bubbles = 2),
            "ggplot")
  expect_equal(ggplot2::quo_name(g$mapping$x), "Fold_Enrichment")
  expect_equal(ggplot2::quo_name(g$mapping$y), "Pathway")
  expect_equal(g$labels$size, "# of DEGs")
  expect_equal(g$labels$colour, "-log10(lowest-p)")
  expect_equal(g$labels$x, "Fold Enrichment")
  expect_equal(g$labels$y, "")

  # change even_breaks
  expect_is(g <- enrichment_chart(RA_clustered,
                                  even_breaks = FALSE),
            "ggplot")
  expect_equal(ggplot2::quo_name(g$mapping$x), "Fold_Enrichment")
  expect_equal(ggplot2::quo_name(g$mapping$y), "Pathway")
  expect_equal(g$labels$size, "# of DEGs")
  expect_equal(g$labels$colour, "-log10(lowest-p)")
  expect_equal(g$labels$x, "Fold Enrichment")
  expect_equal(g$labels$y, "")
})

test_that("enrichment_chart arg checks work", {

  necessary <- c("Pathway", "Fold_Enrichment", "lowest_p",
                 "Up_regulated", "Down_regulated")
  expect_error(enrichment_chart(RA_output[, -2]),
               paste0("The input data frame must have the columns:\n",
                      paste(necessary, collapse = ", ")))

  expect_error(enrichment_chart(RA_output, plot_by_cluster = "WRONG"),
               "`plot_by_cluster` must be either TRUE or FALSE")

  expect_message(enrichment_chart(RA_output, plot_by_cluster = TRUE),
                 "For plotting by cluster, there must a column named `Cluster` in the input data frame!")
})

# term_gene_graph ---------------------------------------------------------
test_that("term_gene_graph produces a ggplot object using the correct data", {

  # Top 10 (default)
  expect_is(p <- term_gene_graph(RA_output), "ggplot")
  expect_equal(sum(p$data$type == "pathway"), 10)

  # Top 3
  expect_is(p <- term_gene_graph(RA_output, num_terms = 3), "ggplot")
  expect_equal(sum(p$data$type == "pathway"), 3)

  # All terms
  expect_is(p <- term_gene_graph(RA_output, num_terms = NULL), "ggplot")
  expect_equal(sum(p$data$type == "pathway"), nrow(RA_output))

  # use_names = TRUE
  expect_is(p <- term_gene_graph(RA_output, use_names = TRUE), "ggplot")
  expect_equal(sum(p$data$type == "pathway"), 10)

  # use_names = "p_val"
  expect_is(p <- term_gene_graph(RA_output, node_size = "p_val"), "ggplot")
  expect_equal(sum(p$data$type == "pathway"), 10)
})

test_that("term_gene_graph arg checks work", {
  expect_error(term_gene_graph(RA_output, num_terms = "WRONG"),
               "`num_terms` must either be numeric or NULL!")

  expect_error(term_gene_graph(RA_output, use_names = "WRONG"),
               "`use_names` must either be TRUE or FALSE!")

  expect_error(term_gene_graph(RA_output, node_size = "WRONG"),
               "`node_size` must either be num_DEGs or p_val!")

  wrong_df <- RA_output[, -c(1,2)]

  ID_column <- "ID"
  necessary_cols <- c("Up_regulated", "Down_regulated", "lowest_p", ID_column)
  expect_error(term_gene_graph(wrong_df, use_names = FALSE),
               paste(c("All of", paste(necessary_cols, collapse = ", "),
                       "must be present in `results_df`!"), collapse = " "))

  ID_column <- "Pathway"
  necessary_cols <- c("Up_regulated", "Down_regulated", "lowest_p", ID_column)
  expect_error(term_gene_graph(wrong_df, use_names = TRUE),
               paste(c("All of", paste(necessary_cols, collapse = ", "),
                       "must be present in `results_df`!"), collapse = " "))

})
