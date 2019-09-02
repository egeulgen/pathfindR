##################################################
## Project: pathfindR
## Script purpose: Testthat testing script for
## core functions
## Date: Sep 2, 2019
## Author: Ege Ulgen
##################################################

# run_pathfindR -----------------------------------------------------------
test_that("run_pathfindR works as expected", {

  ## GR
  expect_is(run_pathfindR(RA_input,
                          iterations = 1,
                          visualize_pathways = FALSE),
            "data.frame")
  expect_is(run_pathfindR(RA_input,
                          iterations = 2,
                          gene_sets = "BioCarta",
                          pin_name_path = "GeneMania"),
            "data.frame")
  ## SA
  expect_is(run_pathfindR(RA_input,
                          iterations = 1,
                          gene_sets = "BioCarta",
                          pin_name_path = "GeneMania",
                          search_method = "SA",
                          visualize_pathways = FALSE),
            "data.frame")

  ## GA
  expect_is(run_pathfindR(RA_input,
                          iterations = 1,
                          gene_sets = "BioCarta",
                          pin_name_path = "GeneMania",
                          search_method = "GA",
                          visualize_pathways = FALSE),
            "data.frame")

  expect_warning(run_pathfindR(RA_input[1:3,],
                               iterations = 1),
                 "Did not find any enriched pathways!")
})

test_that("run_pathfindR arg checks work", {
  expect_error(run_pathfindR(RA_input, search_method = "WRONG"),
               '`search_method` must be one of "GR", "SA", "GA"')

  expect_error(run_pathfindR(RA_input, use_all_positives = "WRONG"),
               "the argument `use_all_positives` must be either TRUE or FALSE")

  expect_error(run_pathfindR(RA_input, silent_option = "WRONG"),
               "the argument `silent_option` must be either TRUE or FALSE")

  expect_error(run_pathfindR(RA_input, gene_sets = "WRONG"),
               "`gene_sets` must be one of KEGG, Reactome, BioCarta, GO-All, GO-BP, GO-CC, GO-MF or Custom")

  expect_error(run_pathfindR(RA_input, gene_sets = "Custom"),
               "You must provide both `custom_genes` and `custom_pathways` if `gene_sets` is `Custom`!")

  expect_error(run_pathfindR(RA_input, bubble = "WRONG"),
               "the argument `bubble` must be either TRUE or FALSE")
})

# return_pin_path ---------------------------------------------------------
test_that("return_pin_path returns the absolute path to PIN file", {

  # default PINs
  expect_true(file.exists(return_pin_path()))
  expect_true(file.exists(return_pin_path("KEGG")))

  # custom PIN
  tmp_file_path <- system.file(paste0("extdata/IntAct.sif"),
                               package = "pathfindR")
  expect_true(file.exists(return_pin_path(tmp_file_path)))


  # invalid custom PIN - wrong format
  tmp_file_path <- system.file(paste0("extdata/MYC.txt"),
                               package = "pathfindR")

  expect_error(return_pin_path(tmp_file_path),
               "The PIN file must have 3 columns and be tab-separated")

  # invalid custom PIN - invalid second column
  custom_sif <- data.frame(P1 = "X", pp = "WRONG", P2 = "Y")
  write.table(custom_sif, "./tmp_sif.sif", sep = "\t",
              col.names = FALSE, row.names = FALSE)
  expect_error(return_pin_path("./tmp_sif.sif"),
               "The second column of the PIN file must all be \"pp\" ")
  unlink("./tmp_sif.sif")

  # invalid option
  expect_error(return_pin_path("WRONG"),
               paste0("The chosen PIN must be one of:\n",
                    "Biogrid, GeneMania, IntAct, KEGG or a valid /path/to/SIF"))
})

# input_testing -----------------------------------------------------------
test_that("input_testing works", {
  expect_message(input_testing(input = RA_input,
                               p_val_threshold = 0.05),
                 "The input looks OK\n\n")

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
test_that("input_processing works", {
  path2pin <- return_pin_path()
  # full df
  expect_is(input_processing(RA_input,
                             p_val_threshold = 0.05,
                             pin_path = path2pin,
                             human_genes = TRUE),
            "data.frame")

  expect_is(input_processing(RA_input[1:10,],
                             p_val_threshold = 0.05,
                             pin_path = path2pin,
                             human_genes = FALSE),
            "data.frame")

  # no change val.s df
  input2 <- RA_input[, -2]
  expect_is(input_processing(input2,
                             p_val_threshold = 0.05,
                             pin_path = path2pin,
                             human_genes = TRUE),
            "data.frame")
})

test_that("input_processing errors and warnings work", {
  input2 <- RA_input[, -2]
  expect_warning(input_processing(input2[1:10,],
                                  p_val_threshold = 0.05,
                                  pin_path = path2pin,
                                  human_genes = TRUE),
                 "The gene column was turned into character from factor.")

  expect_error(input_processing(RA_input,
                                p_val_threshold = 1e-100,
                                pin_path = path2pin),
               "No input p value is lower than the provided threshold \\(1e-100\\)")

  input_dup <- RA_input[1:3, ]
  input_dup <- rbind(input_dup, input_dup[1, ])
  expect_warning(input_processing(input_dup,
                                  p_val_threshold = 5e-2,
                                  pin_path = path2pin),
                 "Duplicated genes found! The lowest p value for each gene was selected")

  tmp_input <- RA_input[1:10,]
  tmp_input$adj.P.Val <- 1e-15
  expect_true(all(input_processing(tmp_input,
                                   p_val_threshold = 5e-2,
                                   pin_path = path2pin)$P_VALUE >= 1e-13))

  tmp_input$Gene.symbol <- paste0(LETTERS[seq_len(nrow(tmp_input))], "WRONG")
  expect_error(input_processing(tmp_input,
                                p_val_threshold = 5e-2,
                                pin_path = path2pin),
          "None of the genes were in the PIN\nPlease check your gene symbols")

  tmp_input$Gene.symbol[1] <- "B"
  tmp_input$Gene.symbol[2] <- "PPIB"
  expect_error(input_processing(tmp_input,
                                p_val_threshold = 5e-2,
                                pin_path = path2pin),
        "After processing, 1 gene \\(or no genes\\) could be mapped to the PIN")
})

# annotate_pathway_DEGs ---------------------------------------------------
example_gene_data <- RA_input[1:50, ]
colnames(example_gene_data) <- c("GENE", "CHANGE", "P_VALUE")
tmp_res <- RA_output[, -c(7, 8)]

test_that("annotate_pathway_DEGs adds input genes for each term", {

  ## Default
  expect_is(annotated_result <- annotate_pathway_DEGs(tmp_res,
                                                      example_gene_data),
            "data.frame")
  expect_true("Up_regulated" %in% colnames(annotated_result) &
                "Down_regulated" %in% colnames(annotated_result))
  expect_true(nrow(annotated_result) == nrow(tmp_res))

  ## Custom
  expect_is(annotated_result <- annotate_pathway_DEGs(tmp_res,
                                                      example_gene_data,
                                                      gene_sets = "Custom",
                                                      custom_genes = kegg_genes),
            "data.frame")
  expect_true("Up_regulated" %in% colnames(annotated_result) &
                "Down_regulated" %in% colnames(annotated_result))
  expect_true(nrow(annotated_result) == nrow(tmp_res))

  expect_error(annotate_pathway_DEGs(tmp_res,
                                     example_gene_data,
                                     gene_sets = "Custom"),
               '`custom_genes` must be provided if `gene_sets = "Custom"`')
})
