##################################################
## Project: pathfindR
## Script purpose: Testthat testing script for
## core functions
## Date: Oct 28, 2019
## Author: Ege Ulgen
##################################################

# run_pathfindR -----------------------------------------------------------
test_that("`run_pathfindR()` works as expected", {
  skip_on_cran()
  ## GR
  expect_is(run_pathfindR(RA_input,
                          iterations = 1,
                          visualize_enriched_terms = FALSE),
            "data.frame")
  expect_is(run_pathfindR(RA_input,
                          iterations = 2,
                          n_processes = 2,
                          gene_sets = "BioCarta",
                          pin_name_path = "GeneMania",
                          plot_enrichment_chart = FALSE),
            "data.frame")

  skip("will test SA and GA if we can create a suitable (faster and non-empty) test case")
  ## SA
  expect_is(run_pathfindR(RA_input,
                          iterations = 1,
                          gene_sets = "GO-BP",
                          pin_name_path = "GeneMania",
                          search_method = "SA",
                          visualize_enriched_terms = FALSE,
                          plot_enrichment_chart = FALSE),
            "data.frame")

  ## GA
  expect_is(suppressWarnings(run_pathfindR(RA_input,
                                           iterations = 1,
                                           gene_sets = "BioCarta",
                                           pin_name_path = "GeneMania",
                                           search_method = "GA",
                                           visualize_enriched_terms = FALSE,
                                           plot_enrichment_chart = FALSE)),
            "data.frame")
})

test_that("Expect warning with empty result from `run_pathfindR()`", {
  expect_warning(res <- run_pathfindR(RA_input[1:3, ],
                                      iterations = 1),
                 "Did not find any enriched terms!")
  expect_identical(res, data.frame())
})


test_that("`run_pathfindR()` arg checks work", {
  expect_error(run_pathfindR(RA_input, search_method = "WRONG"),
               '`search_method` must be one of "GR", "SA", "GA"')

  expect_error(run_pathfindR(RA_input, use_all_positives = "WRONG"),
               "the argument `use_all_positives` must be either TRUE or FALSE")

  expect_error(run_pathfindR(RA_input, silent_option = "WRONG"),
    "the argument `silent_option` must be either TRUE or FALSE")

  all_gs_opts <- c("KEGG", "Reactome", "BioCarta",
                   "GO-All", "GO-BP", "GO-CC", "GO-MF",
                   "mmu_KEGG", "Custom")
  expect_error(run_pathfindR(RA_input, gene_sets = "WRONG"),
               paste0("`gene_sets` should be one of ", paste(dQuote(all_gs_opts), collapse = ", ")))

  expect_error(run_pathfindR(RA_input, gene_sets = "Custom"),
               "`custom_genes` and `custom_descriptions` must be provided if `gene_sets = \"Custom\"`")

  expect_error(run_pathfindR(RA_input, plot_enrichment_chart = "WRONG"),
               "the argument `plot_enrichment_chart` must be either TRUE or FALSE")
})

# fetch_gene_set ----------------------------------------------------------
test_that("`fetch_gene_set()` can fetch all gene set objects", {
  ###### KEGG
  expect_is(gset_obj <- fetch_gene_set(gene_sets = "KEGG",
                                       min_gset_size = 10,
                                       max_gset_size = 300),
            "list")
  expect_is(gset_obj$genes_by_term, "list")
  expect_is(gset_obj$term_descriptions, "character")
  expect_true(length(gset_obj$genes_by_term) == length(gset_obj$term_descriptions))
  tmp <- vapply(gset_obj$genes_by_term, length, 1L)
  expect_true(min(tmp) >= 10 & max(tmp) <= 300)

  ###### mmu KEGG
  expect_is(gset_obj <- fetch_gene_set(gene_sets = "mmu_KEGG",
                                       min_gset_size = 10,
                                       max_gset_size = 300),
            "list")
  expect_is(gset_obj$genes_by_term, "list")
  expect_is(gset_obj$term_descriptions, "character")
  expect_true(length(gset_obj$genes_by_term) == length(gset_obj$term_descriptions))
  tmp <- vapply(gset_obj$genes_by_term, length, 1L)
  expect_true(min(tmp) >= 10 & max(tmp) <= 300)

  ###### Reactome
  expect_is(gset_obj <- fetch_gene_set(gene_sets = "Reactome",
                                       min_gset_size = 10,
                                       max_gset_size = 300),
            "list")
  expect_is(gset_obj$genes_by_term, "list")
  expect_is(gset_obj$term_descriptions, "character")
  expect_true(length(gset_obj$genes_by_term) == length(gset_obj$term_descriptions))
  tmp <- vapply(gset_obj$genes_by_term, length, 1L)
  expect_true(min(tmp) >= 10 & max(tmp) <= 300)

  ###### BioCarta
  expect_is(gset_obj <- fetch_gene_set(gene_sets = "BioCarta",
                                       min_gset_size = 10,
                                       max_gset_size = 300),
            "list")
  expect_is(gset_obj$genes_by_term, "list")
  expect_is(gset_obj$term_descriptions, "character")
  expect_true(length(gset_obj$genes_by_term) == length(gset_obj$term_descriptions))
  tmp <- vapply(gset_obj$genes_by_term, length, 1L)
  expect_true(min(tmp) >= 10 & max(tmp) <= 300)

  ###### GO-All
  expect_is(gset_obj <- fetch_gene_set(gene_sets = "GO-All",
                                       min_gset_size = 10,
                                       max_gset_size = 300),
            "list")
  expect_is(gset_obj$genes_by_term, "list")
  expect_is(gset_obj$term_descriptions, "character")
  expect_true(length(gset_obj$genes_by_term) == length(gset_obj$term_descriptions))
  tmp <- vapply(gset_obj$genes_by_term, length, 1L)
  expect_true(min(tmp) >= 10 & max(tmp) <= 300)

  ###### GO-BP
  expect_is(gset_obj <- fetch_gene_set(gene_sets = "GO-BP",
                                       min_gset_size = 10,
                                       max_gset_size = 300),
            "list")
  expect_is(gset_obj$genes_by_term, "list")
  expect_is(gset_obj$term_descriptions, "character")
  expect_true(length(gset_obj$genes_by_term) == length(gset_obj$term_descriptions))
  tmp <- vapply(gset_obj$genes_by_term, length, 1L)
  expect_true(min(tmp) >= 10 & max(tmp) <= 300)

  ###### GO-CC
  expect_is(gset_obj <- fetch_gene_set(gene_sets = "GO-CC",
                                       min_gset_size = 10,
                                       max_gset_size = 300),
            "list")
  expect_is(gset_obj$genes_by_term, "list")
  expect_is(gset_obj$term_descriptions, "character")
  expect_true(length(gset_obj$genes_by_term) == length(gset_obj$term_descriptions))
  tmp <- vapply(gset_obj$genes_by_term, length, 1L)
  expect_true(min(tmp) >= 10 & max(tmp) <= 300)

  ###### GO-MF
  expect_is(gset_obj <- fetch_gene_set(gene_sets = "GO-MF",
                                       min_gset_size = 10,
                                       max_gset_size = 300),
            "list")
  expect_is(gset_obj$genes_by_term, "list")
  expect_is(gset_obj$term_descriptions, "character")
  expect_true(length(gset_obj$genes_by_term) == length(gset_obj$term_descriptions))
  tmp <- vapply(gset_obj$genes_by_term, length, 1L)
  expect_true(min(tmp) >= 10 & max(tmp) <= 300)

  ###### Custom
  fetch_gene_set(gene_sets = "Custom",
                 min_gset_size = 10,
                 max_gset_size = 300,
                 custom_genes = kegg_genes,
                 custom_descriptions = kegg_descriptions)
  expect_is(gset_obj$genes_by_term, "list")
  expect_is(gset_obj$term_descriptions, "character")
  expect_true(length(gset_obj$genes_by_term) == length(gset_obj$term_descriptions))
  tmp <- vapply(gset_obj$genes_by_term, length, 1L)
  expect_true(min(tmp) >= 10 & max(tmp) <= 300)
})

test_that("min/max_gset_size args in `fetch_gene_set()` correctly filter gene sets", {
  expect_is(gset_obj1 <- fetch_gene_set(gene_sets = "KEGG",
                                        min_gset_size = 10,
                                        max_gset_size = 300),
            "list")
  tmp <- vapply(gset_obj1$genes_by_term, length, 1L)
  expect_true(min(tmp) >= 10 & max(tmp) <= 300)

  expect_is(gset_obj2 <- fetch_gene_set(gene_sets = "KEGG",
                                        min_gset_size = 50,
                                        max_gset_size = 200),
            "list")
  tmp <- vapply(gset_obj2$genes_by_term, length, 1L)
  expect_true(min(tmp) >= 50 & max(tmp) <= 200)
  expect_true(length(gset_obj2$genes_by_term) < length(gset_obj1$genes_by_term))
})

test_that("In `fetch_gene_set()`, for 'Custom' gene set, check if the custom objects are provided", {
  expect_error(fetch_gene_set(gene_sets = "Custom"),
               "`custom_genes` and `custom_descriptions` must be provided if `gene_sets = \"Custom\"`")
  expect_error(fetch_gene_set(gene_sets = "Custom",
                              custom_genes = kegg_genes),
               "`custom_genes` and `custom_descriptions` must be provided if `gene_sets = \"Custom\"`")
  expect_error(fetch_gene_set(gene_sets = "Custom",
                              custom_descriptions = kegg_descriptions),
               "`custom_genes` and `custom_descriptions` must be provided if `gene_sets = \"Custom\"`")
})

# return_pin_path ---------------------------------------------------------
test_that("`return_pin_path()` returns the absolute path to PIN file", {

  # default PINs
  expect_true(file.exists(return_pin_path("Biogrid")))
  expect_true(file.exists(return_pin_path("GeneMania")))
  expect_true(file.exists(return_pin_path("IntAct")))
  expect_true(file.exists(return_pin_path("KEGG")))
  expect_true(file.exists(return_pin_path("mmu_STRING")))

  # custom PIN
  custom_sif_path <- file.path(tempdir(check = TRUE), "custom.sif")
  download.file("http://www.biofabric.org/sifFiles/full100.sif",
                custom_sif_path, quiet = TRUE)
  expect_true(file.exists(return_pin_path(custom_sif_path)))

  # invalid custom PIN - wrong format
  invalid_sif_path <- system.file(paste0("extdata/MYC.txt"),
                                  package = "pathfindR")
  expect_error(return_pin_path(invalid_sif_path),
               "The PIN file must have 3 columns and be tab-separated")

  # invalid custom PIN - invalid second column
  invalid_sif_path <- file.path(tempdir(check = TRUE), "custom.sif")
  invalid_custom_sif <- data.frame(P1 = "X", pp = "WRONG", P2 = "Y")
  write.table(invalid_custom_sif, invalid_sif_path, sep = "\t",
              col.names = FALSE, row.names = FALSE)
  expect_error(return_pin_path(invalid_sif_path),
               "The second column of the PIN file must all be \"pp\" ")

  # invalid option
  valid_opts <- c("Biogrid", "GeneMania", "IntAct", "KEGG", "mmu_STRING", "/path/to/custom/SIF")
  expect_error(return_pin_path("WRONG"),
               paste0("The chosen PIN must be one of:\n",
                      paste(dQuote(valid_opts), collapse = ", ")))
})

# input_testing -----------------------------------------------------------
test_that("`input_testing()` works", {
  expect_message(input_testing(input = RA_input,
                               p_val_threshold = 0.05),
                 "The input looks OK")

  expect_error(input_testing(input = matrix(),
                             p_val_threshold = 0.05),
               "the input is not a data frame")

  expect_error(input_testing(input = RA_input[1, ],
                             p_val_threshold = 0.05),
            "There must be at least 2 rows \\(genes\\) in the input data frame")

  expect_error(input_testing(input = RA_input[, 1, drop = FALSE],
                             p_val_threshold = 0.05),
               "There must be at least 2 columns in the input data frame")

  expect_error(input_testing(input = RA_input,
                             p_val_threshold = "WRONG"),
               "`p_val_threshold` must be a numeric value between 0 and 1")

  expect_error(input_testing(input = RA_input,
                             p_val_threshold = -1),
               "`p_val_threshold` must be between 0 and 1")

  tmp <- RA_input
  tmp$adj.P.Val <- NA
  expect_error(input_testing(input = tmp,
                             p_val_threshold = 0.05),
               "p values cannot contain NA values")

  tmp <- RA_input
  tmp$adj.P.Val <- "WRONG"
  expect_error(input_testing(input = tmp,
                             p_val_threshold = 0.05),
               "p values must all be numeric")

  tmp <- RA_input
  tmp$adj.P.Val[1] <- -1
  expect_error(input_testing(input = tmp,
                             p_val_threshold = 0.05),
               "p values must all be between 0 and 1")
})

# input_processing --------------------------------------------------------
test_that("`input_processing()` works", {
  # full df
  expect_is(tmp <- input_processing(input = RA_input,
                                    p_val_threshold = 0.05,
                                    pin_name_path = "Biogrid",
                                    convert2alias = TRUE),
            "data.frame")
  expect_true(ncol(tmp) == 4)
  expect_true(nrow(tmp) <= nrow(RA_input))

  expect_is(input_processing(RA_input[1:10, ],
                             p_val_threshold = 0.01,
                             pin_name_path = "Biogrid",
                             convert2alias = FALSE),
            "data.frame")

  # no change values provided
  input2 <- RA_input[, -2]
  expect_is(suppressWarnings(input_processing(input2,
                                              p_val_threshold = 0.05,
                                              pin_name_path = "Biogrid",
                                              convert2alias = TRUE)),
    "data.frame")

  # multiple mapping
  input_m <- RA_input
  input_m$Gene.symbol[1] <- "GIG24"
  input_m$Gene.symbol[2] <- "ACT"
  input_m$Gene.symbol[3] <- "AACT"
  input_m$Gene.symbol[4] <- "GIG25"
  expect_is(input_processing(input_m,
                             p_val_threshold = 0.05,
                             pin_name_path = "Biogrid",
                             convert2alias = TRUE),
            "data.frame")
})

test_that("`input_processing()` errors and warnings work", {
  input2 <- RA_input[1:10, ]
  input2$Gene.symbol <- as.factor(input2$Gene.symbol)
  expect_warning(input_processing(input2,
                                  p_val_threshold = 0.05,
                                  pin_name_path = "Biogrid",
                                  convert2alias = TRUE),
                 "The gene column was turned into character from factor.")

  expect_error(input_processing(RA_input,
                                p_val_threshold = 1e-100,
                                pin_name_path = "Biogrid"),
               "No input p value is lower than the provided threshold \\(1e-100\\)")

  input_dup <- RA_input[1:3, ]
  input_dup <- rbind(input_dup, input_dup[1, ])
  expect_warning(input_processing(input_dup,
                                  p_val_threshold = 5e-2,
                                  pin_name_path = "Biogrid"),
                 "Duplicated genes found! The lowest p value for each gene was selected")

  tmp_input <- RA_input[1:10, ]
  tmp_input$adj.P.Val <- 1e-15
  expect_message(tmp <- input_processing(tmp_input,
                                         p_val_threshold = 5e-2,
                                         pin_name_path = "Biogrid"),
                 "pathfindR cannot handle p values < 1e-13. These were changed to 1e-13")
  expect_true(all(tmp$P_VALUE == 1e-13))

  tmp_input$Gene.symbol <- paste0(LETTERS[seq_len(nrow(tmp_input))], "WRONG")
  expect_error(input_processing(tmp_input,
                                p_val_threshold = 5e-2,
                                pin_name_path = "Biogrid"),
               "None of the genes were in the PIN\nPlease check your gene symbols")

  tmp_input$Gene.symbol[1] <- "NRD1"
  tmp_input$Gene.symbol[2] <- "NRDC"
  expect_error(input_processing(tmp_input,
                                p_val_threshold = 5e-2,
                                pin_name_path = "Biogrid"),
               "After processing, 1 gene \\(or no genes\\) could be mapped to the PIN"
  )
})

# annotate_term_genes -----------------------------------------------------
example_gene_data <- RA_input[1:50, ]
colnames(example_gene_data) <- c("GENE", "CHANGE", "P_VALUE")
tmp_res <- RA_output[, -c(7, 8)]

test_that("`annotate_term_genes()` adds input genes for each term", {
  expect_is(annotated_result <- annotate_term_genes(result_df = tmp_res,
                                                    input_processed = example_gene_data,
                                                    genes_by_term = kegg_genes),
            "data.frame")
  expect_true("Up_regulated" %in% colnames(annotated_result) &
                "Down_regulated" %in% colnames(annotated_result))
  expect_true(nrow(annotated_result) == nrow(tmp_res))

})
