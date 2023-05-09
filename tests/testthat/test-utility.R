##################################################
## Package: pathfindR
## Script purpose: Unit testing script for
## utility functions
## Date: Apr 27, 2023
## Author: Ege Ulgen
##################################################

# active_snw_enrichment_wrapper -------------------------------------------
input_processed <- input_processing(example_pathfindR_input[1:2, ])
pin_path <- return_pin_path()
gset_list <- fetch_gene_set()

test_that("`active_snw_enrichment_wrapper()` works as expected", {
  org_dir <- getwd()
  test_directory <- file.path(tempdir(check = TRUE), "snw_wrapper_test")
  dir.create(test_directory)
  setwd(test_directory)
  on.exit(setwd(org_dir))
  on.exit(unlink(test_directory), add = TRUE)

  expect_is(
    pathfindR:::active_snw_enrichment_wrapper(
      input_processed = input_processed,
      pin_path = pin_path,
      gset_list = gset_list,
      enrichment_threshold = 0.05,
      list_active_snw_genes = FALSE,
      iterations = 1
    ),
    "NULL"
  )

  skip_on_cran()
  expect_is(
    pathfindR:::active_snw_enrichment_wrapper(
      input_processed = input_processed,
      pin_path = pin_path,
      gset_list = gset_list,
      enrichment_threshold = 0.05,
      list_active_snw_genes = FALSE,
      iterations = 2
    ),
    "NULL"
  )
})

test_that("`active_snw_enrichment_wrapper()` argument checks work", {
  valid_mets <- c("GR", "SA", "GA")
  expect_error(
    pathfindR:::active_snw_enrichment_wrapper(
      input_processed = input_processed,
      pin_path = pin_path,
      gset_list = gset_list,
      enrichment_threshold = 0.05,
      list_active_snw_genes = FALSE,
      search_method = "INVALID"
    ),
    paste0(
      "`search_method` should be one of ",
      paste(dQuote(valid_mets), collapse = ", ")
    )
  )

  expect_error(
    pathfindR:::active_snw_enrichment_wrapper(
      input_processed = input_processed,
      pin_path = pin_path,
      gset_list = gset_list,
      enrichment_threshold = 0.05,
      list_active_snw_genes = FALSE,
      use_all_positives = "INVALID"
    ),
    "`use_all_positives` should be either TRUE or FALSE"
  )

  expect_error(
    pathfindR:::active_snw_enrichment_wrapper(
      input_processed = input_processed,
      pin_path = pin_path,
      gset_list = gset_list,
      enrichment_threshold = 0.05,
      list_active_snw_genes = FALSE,
      silent_option = "INVALID"
    ),
    "`silent_option` should be either TRUE or FALSE"
  )

  expect_error(
    pathfindR:::active_snw_enrichment_wrapper(
      input_processed = input_processed,
      pin_path = pin_path,
      gset_list = gset_list,
      enrichment_threshold = 0.05,
      list_active_snw_genes = FALSE,
      iterations = "INVALID"
    ),
    "`iterations` should be a positive integer"
  )

  expect_error(
    pathfindR:::active_snw_enrichment_wrapper(
      input_processed = input_processed,
      pin_path = pin_path,
      gset_list = gset_list,
      enrichment_threshold = 0.05,
      list_active_snw_genes = FALSE,
      iterations = 0
    ),
    "`iterations` should be >= 1"
  )

  expect_error(
    pathfindR:::active_snw_enrichment_wrapper(
      input_processed = input_processed,
      pin_path = pin_path,
      gset_list = gset_list,
      enrichment_threshold = 0.05,
      list_active_snw_genes = FALSE,
      n_processes = "INVALID"
    ),
    "`n_processes` should be either NULL or a positive integer"
  )

  expect_error(
    pathfindR:::active_snw_enrichment_wrapper(
      input_processed = input_processed,
      pin_path = pin_path,
      gset_list = gset_list,
      enrichment_threshold = 0.05,
      list_active_snw_genes = FALSE,
      n_processes = 0
    ),
    "`n_processes` should be > 1"
  )
})


# configure_output_dir ----------------------------------------------------
test_that("`configure_output_dir()` works as expected", {
  expect_equal(
    pathfindR:::configure_output_dir(),
    file.path(tempdir(), "pathfindR_results")
  )

  test_out_dir <- file.path(tempdir(check = TRUE), "TEST")
  for (i in 1:3) {
    actual_dir <- pathfindR:::configure_output_dir(test_out_dir)
    dir_to_check <- test_out_dir
    if (i > 1) {
      dir_to_check <- paste0(dir_to_check, "(", i - 1, ")")
    }
    expect_equal(actual_dir, dir_to_check)
    dir.create(actual_dir)
  }
})

# create_HTML_report ------------------------------------------------------
test_that("`create_HTML_report()` works a expected", {
  test_directory <- tempdir(check = TRUE)

  input_processed <- input_processing(example_pathfindR_input)

  create_HTML_report(
    input = example_pathfindR_input,
    input_processed = input_processed,
    final_res = example_pathfindR_output,
    dir_for_report = test_directory
  )

  expect_true(file.exists(file.path(test_directory, "results.html")))
  expect_true(file.exists(file.path(test_directory, "enriched_terms.html")))
  expect_true(file.exists(file.path(test_directory, "conversion_table.html")))
})


# fetch_gene_set ----------------------------------------------------------
test_that("`fetch_gene_set()` can fetch all gene set objects", {
  skip_on_cran()
  ###### KEGG
  expect_is(
    gset_obj <- fetch_gene_set(
      gene_sets = "KEGG",
      min_gset_size = 10,
      max_gset_size = 300
    ),
    "list"
  )
  expect_is(gset_obj$genes_by_term, "list")
  expect_is(gset_obj$term_descriptions, "character")
  expect_true(length(gset_obj$genes_by_term) == length(gset_obj$term_descriptions))
  tmp <- vapply(gset_obj$genes_by_term, length, 1L)
  expect_true(min(tmp) >= 10 & max(tmp) <= 300)

  ###### mmu KEGG
  expect_is(
    gset_obj <- fetch_gene_set(
      gene_sets = "mmu_KEGG",
      min_gset_size = 10,
      max_gset_size = 300
    ),
    "list"
  )
  expect_is(gset_obj$genes_by_term, "list")
  expect_is(gset_obj$term_descriptions, "character")
  expect_true(length(gset_obj$genes_by_term) == length(gset_obj$term_descriptions))
  tmp <- vapply(gset_obj$genes_by_term, length, 1L)
  expect_true(min(tmp) >= 10 & max(tmp) <= 300)

  ###### Reactome
  expect_is(
    gset_obj <- fetch_gene_set(
      gene_sets = "Reactome",
      min_gset_size = 10,
      max_gset_size = 300
    ),
    "list"
  )
  expect_is(gset_obj$genes_by_term, "list")
  expect_is(gset_obj$term_descriptions, "character")
  expect_true(length(gset_obj$genes_by_term) == length(gset_obj$term_descriptions))
  tmp <- vapply(gset_obj$genes_by_term, length, 1L)
  expect_true(min(tmp) >= 10 & max(tmp) <= 300)

  ###### BioCarta
  expect_is(
    gset_obj <- fetch_gene_set(
      gene_sets = "BioCarta",
      min_gset_size = 10,
      max_gset_size = 300
    ),
    "list"
  )
  expect_is(gset_obj$genes_by_term, "list")
  expect_is(gset_obj$term_descriptions, "character")
  expect_true(length(gset_obj$genes_by_term) == length(gset_obj$term_descriptions))
  tmp <- vapply(gset_obj$genes_by_term, length, 1L)
  expect_true(min(tmp) >= 10 & max(tmp) <= 300)


  ###### cell_markers
  expect_is(
    gset_obj <- fetch_gene_set(
      gene_sets = "cell_markers",
      min_gset_size = 10,
      max_gset_size = 300
    ),
    "list"
  )
  expect_is(gset_obj$genes_by_term, "list")
  expect_is(gset_obj$term_descriptions, "character")
  expect_true(length(gset_obj$genes_by_term) == length(gset_obj$term_descriptions))
  tmp <- vapply(gset_obj$genes_by_term, length, 1L)
  expect_true(min(tmp) >= 10 & max(tmp) <= 300)

  ###### GO-All
  expect_is(
    gset_obj <- fetch_gene_set(
      gene_sets = "GO-All",
      min_gset_size = 10,
      max_gset_size = 300
    ),
    "list"
  )
  expect_is(gset_obj$genes_by_term, "list")
  expect_is(gset_obj$term_descriptions, "character")
  expect_true(length(gset_obj$genes_by_term) == length(gset_obj$term_descriptions))
  tmp <- vapply(gset_obj$genes_by_term, length, 1L)
  expect_true(min(tmp) >= 10 & max(tmp) <= 300)

  ###### GO-BP
  expect_is(
    gset_obj <- fetch_gene_set(
      gene_sets = "GO-BP",
      min_gset_size = 10,
      max_gset_size = 300
    ),
    "list"
  )
  expect_is(gset_obj$genes_by_term, "list")
  expect_is(gset_obj$term_descriptions, "character")
  expect_true(length(gset_obj$genes_by_term) == length(gset_obj$term_descriptions))
  tmp <- vapply(gset_obj$genes_by_term, length, 1L)
  expect_true(min(tmp) >= 10 & max(tmp) <= 300)

  ###### GO-CC
  expect_is(
    gset_obj <- fetch_gene_set(
      gene_sets = "GO-CC",
      min_gset_size = 10,
      max_gset_size = 300
    ),
    "list"
  )
  expect_is(gset_obj$genes_by_term, "list")
  expect_is(gset_obj$term_descriptions, "character")
  expect_true(length(gset_obj$genes_by_term) == length(gset_obj$term_descriptions))
  tmp <- vapply(gset_obj$genes_by_term, length, 1L)
  expect_true(min(tmp) >= 10 & max(tmp) <= 300)

  ###### GO-MF
  expect_is(
    gset_obj <- fetch_gene_set(
      gene_sets = "GO-MF",
      min_gset_size = 10,
      max_gset_size = 300
    ),
    "list"
  )
  expect_is(gset_obj$genes_by_term, "list")
  expect_is(gset_obj$term_descriptions, "character")
  expect_true(length(gset_obj$genes_by_term) == length(gset_obj$term_descriptions))
  tmp <- vapply(gset_obj$genes_by_term, length, 1L)
  expect_true(min(tmp) >= 10 & max(tmp) <= 300)

  ###### Custom
  gset_obj <- fetch_gene_set(
    gene_sets = "Custom",
    min_gset_size = 20,
    max_gset_size = 200,
    custom_genes = kegg_genes,
    custom_descriptions = kegg_descriptions
  )
  expect_is(gset_obj$genes_by_term, "list")
  expect_is(gset_obj$term_descriptions, "character")
  expect_true(length(gset_obj$genes_by_term) == length(gset_obj$term_descriptions))
  tmp <- vapply(gset_obj$genes_by_term, length, 1L)
  expect_true(min(tmp) >= 20 & max(tmp) <= 200)
})

test_that("min/max_gset_size args in `fetch_gene_set()` correctly filter gene sets", {
  skip_on_cran()
  expect_is(
    gset_obj1 <- fetch_gene_set(
      gene_sets = "KEGG",
      min_gset_size = 10,
      max_gset_size = 300
    ),
    "list"
  )
  tmp <- vapply(gset_obj1$genes_by_term, length, 1L)
  expect_true(min(tmp) >= 10 & max(tmp) <= 300)

  expect_is(
    gset_obj2 <- fetch_gene_set(
      gene_sets = "KEGG",
      min_gset_size = 50,
      max_gset_size = 200
    ),
    "list"
  )
  tmp <- vapply(gset_obj2$genes_by_term, length, 1L)
  expect_true(min(tmp) >= 50 & max(tmp) <= 200)
  expect_true(length(gset_obj2$genes_by_term) < length(gset_obj1$genes_by_term))
})

test_that("In `fetch_gene_set()`, for 'Custom' gene set, check if the custom objects are provided", {
  expect_error(
    fetch_gene_set(gene_sets = "Custom"),
    "`custom_genes` and `custom_descriptions` must be provided if `gene_sets = \"Custom\"`"
  )
  expect_error(
    fetch_gene_set(
      gene_sets = "Custom",
      custom_genes = kegg_genes
    ),
    "`custom_genes` and `custom_descriptions` must be provided if `gene_sets = \"Custom\"`"
  )
  expect_error(
    fetch_gene_set(
      gene_sets = "Custom",
      custom_descriptions = kegg_descriptions
    ),
    "`custom_genes` and `custom_descriptions` must be provided if `gene_sets = \"Custom\"`"
  )
})

test_that("`fetch_gene_set()` argument checks work", {
  all_gs_opts <- c(
    "KEGG", "Reactome", "BioCarta",
    "GO-All", "GO-BP", "GO-CC", "GO-MF",
    "cell_markers", "mmu_KEGG", "Custom"
  )
  expect_error(
    fetch_gene_set(gene_sets = "INVALID"),
    paste0(
      "`gene_sets` should be one of ",
      paste(dQuote(all_gs_opts), collapse = ", ")
    )
  )

  expect_error(
    fetch_gene_set(min_gset_size = "INVALID"),
    "`min_gset_size` should be numeric"
  )

  expect_error(
    fetch_gene_set(max_gset_size = "INVALID"),
    "`max_gset_size` should be numeric"
  )

  expect_error(
    fetch_gene_set(
      gene_sets = "Custom",
      custom_genes = "INVALID",
      custom_descriptions = ""
    ),
    "`custom_genes` should be a list of term gene sets"
  )
  expect_error(
    fetch_gene_set(
      gene_sets = "Custom",
      custom_genes = list(),
      custom_descriptions = ""
    ),
    "`custom_genes` should be a named list \\(names are gene set IDs\\)"
  )

  expect_error(
    fetch_gene_set(
      gene_sets = "Custom",
      custom_genes = kegg_genes,
      custom_descriptions = list()
    ),
    "`custom_descriptions` should be a vector of term gene descriptions"
  )
  expect_error(
    fetch_gene_set(
      gene_sets = "Custom",
      custom_genes = kegg_genes,
      custom_descriptions = 1:3
    ),
    "`custom_descriptions` should be a named vector \\(names are gene set IDs\\)"
  )
})

# return_pin_path ---------------------------------------------------------
test_that("`return_pin_path()` returns the absolute path to PIN file", {
  skip_on_cran()
  # default PINs
  expect_true(file.exists(return_pin_path("Biogrid")))
  expect_true(file.exists(return_pin_path("STRING")))
  expect_true(file.exists(return_pin_path("GeneMania")))
  expect_true(file.exists(return_pin_path("IntAct")))
  expect_true(file.exists(return_pin_path("KEGG")))
  expect_true(file.exists(return_pin_path("mmu_STRING")))

  # custom PIN
  custom_pin <- read.delim(return_pin_path("KEGG"),
    header = FALSE,
    stringsAsFactors = FALSE
  )
  custom_pin <- custom_pin[1:10, ]
  custom_pin$V1 <- tolower(custom_pin$V1)
  custom_sif_path <- file.path(tempdir(check = TRUE), "tmp_PIN.sif")
  utils::write.table(custom_pin,
    custom_sif_path,
    sep = "\t",
    row.names = FALSE, col.names = FALSE, quote = FALSE
  )
  expect_true(file.exists(return_pin_path(custom_sif_path)))
  # convert to uppercase works
  upper_case_custom <- read.delim(return_pin_path(custom_sif_path), header = FALSE)
  expect_true(all(toupper(upper_case_custom[, 1]) == upper_case_custom[, 1]))
  expect_true(all(toupper(upper_case_custom[, 3]) == upper_case_custom[, 3]))


  # invalid custom PIN - wrong format
  invalid_sif_path <- system.file(paste0("extdata/MYC.txt"),
    package = "pathfindR"
  )
  expect_error(
    return_pin_path(invalid_sif_path),
    "The PIN file must have 3 columns and be tab-separated"
  )

  # invalid custom PIN - invalid second column
  invalid_sif_path <- file.path(tempdir(check = TRUE), "custom.sif")
  invalid_custom_sif <- data.frame(P1 = "X", pp = "INVALID", P2 = "Y")
  write.table(invalid_custom_sif, invalid_sif_path,
    sep = "\t",
    col.names = FALSE, row.names = FALSE
  )
  expect_error(
    return_pin_path(invalid_sif_path),
    "The second column of the PIN file must all be \"pp\" "
  )

  # invalid option
  valid_opts <- c(
    "Biogrid", "STRING", "GeneMania", "IntAct", "KEGG",
    "mmu_STRING", "/path/to/custom/SIF"
  )
  expect_error(
    return_pin_path("INVALID"),
    paste0(
      "The chosen PIN must be one of:\n",
      paste(dQuote(valid_opts), collapse = ", ")
    )
  )
})

# input_testing -----------------------------------------------------------
test_that("`input_testing()` works", {
  expect_message(
    input_testing(
      input = example_pathfindR_input,
      p_val_threshold = 0.05
    ),
    "The input looks OK"
  )

  expect_error(
    input_testing(
      input = matrix(),
      p_val_threshold = 0.05
    ),
    "the input is not a data frame"
  )

  expect_error(
    input_testing(
      input = example_pathfindR_input[, 1, drop = FALSE],
      p_val_threshold = 0.05
    ),
    "the input should have 2 or 3 columns"
  )

  expect_error(
    input_testing(
      input = example_pathfindR_input[1, ],
      p_val_threshold = 0.05
    ),
    "There must be at least 2 rows \\(genes\\) in the input data frame"
  )

  expect_error(
    input_testing(
      input = example_pathfindR_input,
      p_val_threshold = "INVALID"
    ),
    "`p_val_threshold` must be a numeric value between 0 and 1"
  )

  expect_error(
    input_testing(
      input = example_pathfindR_input,
      p_val_threshold = -1
    ),
    "`p_val_threshold` must be between 0 and 1"
  )

  tmp <- example_pathfindR_input
  tmp$adj.P.Val <- NA
  expect_error(
    input_testing(
      input = tmp,
      p_val_threshold = 0.05
    ),
    "p values cannot contain NA values"
  )

  tmp <- example_pathfindR_input
  tmp$adj.P.Val <- "INVALID"
  expect_error(
    input_testing(
      input = tmp,
      p_val_threshold = 0.05
    ),
    "p values must all be numeric"
  )

  tmp <- example_pathfindR_input
  tmp$adj.P.Val[1] <- -1
  expect_error(
    input_testing(
      input = tmp,
      p_val_threshold = 0.05
    ),
    "p values must all be between 0 and 1"
  )
})

# input_processing --------------------------------------------------------
test_that("`input_processing()` works", {
  skip_on_cran()
  expect_is(
    tmp <- input_processing(
      input = example_pathfindR_input[1:5, ],
      p_val_threshold = 0.05,
      pin_name_path = "KEGG",
      convert2alias = TRUE
    ),
    "data.frame"
  )
  expect_true(ncol(tmp) == 4)
  expect_true(nrow(tmp) <= nrow(example_pathfindR_input))

  expect_is(
    input_processing(example_pathfindR_input[5:10, ],
      p_val_threshold = 0.01,
      pin_name_path = "KEGG",
      convert2alias = FALSE
    ),
    "data.frame"
  )

  # no change values provided
  input2 <- example_pathfindR_input[1:5, -2]
  expect_is(
    tmp <- suppressWarnings(input_processing(input2,
      p_val_threshold = 0.05,
      pin_name_path = "Biogrid",
      convert2alias = TRUE
    )),
    "data.frame"
  )
  expect_true(all(tmp$CHANGE == 1e6))

  # multiple mapping
  input_m <- example_pathfindR_input[1:5, ]
  input_m$Gene.symbol[1] <- "GIG24"
  input_m$Gene.symbol[2] <- "ACT"
  input_m$Gene.symbol[3] <- "AACT"
  input_m$Gene.symbol[4] <- "GIG25"
  expect_is(
    tmp <- input_processing(input_m,
      p_val_threshold = 0.05,
      pin_name_path = "Biogrid",
      convert2alias = TRUE
    ),
    "data.frame"
  )
})

test_that("`input_processing()` errors and warnings work", {
  input2 <- example_pathfindR_input[1:5, ]
  input2$Gene.symbol <- as.factor(input2$Gene.symbol)
  expect_warning(
    input_processing(input2,
      p_val_threshold = 0.05,
      pin_name_path = "Biogrid",
      convert2alias = TRUE
    ),
    "The gene column was turned into character from factor."
  )

  expect_error(
    input_processing(example_pathfindR_input,
      p_val_threshold = 1e-100,
      pin_name_path = "Biogrid"
    ),
    "No input p value is lower than the provided threshold \\(1e-100\\)"
  )

  input_dup <- example_pathfindR_input[1:3, ]
  input_dup <- rbind(input_dup, input_dup[1, ])
  expect_warning(
    input_processing(input_dup,
      p_val_threshold = 5e-2,
      pin_name_path = "Biogrid"
    ),
    "Duplicated genes found! The lowest p value for each gene was selected"
  )

  tmp_input <- example_pathfindR_input[1:3, ]
  tmp_input$adj.P.Val <- 1e-15
  expect_message(
    tmp <- input_processing(tmp_input,
      p_val_threshold = 5e-2,
      pin_name_path = "Biogrid"
    ),
    "pathfindR cannot handle p values < 1e-13. These were changed to 1e-13"
  )
  expect_true(all(tmp$P_VALUE == 1e-13))

  tmp_input$Gene.symbol <- paste0(LETTERS[seq_len(nrow(tmp_input))], "INVALID")
  expect_error(
    input_processing(tmp_input,
      p_val_threshold = 5e-2,
      pin_name_path = "Biogrid"
    ),
    "None of the genes were in the PIN\nPlease check your gene symbols"
  )

  tmp_input$Gene.symbol[1] <- "NRD1"
  tmp_input$Gene.symbol[2] <- "NRDC"
  expect_error(
    input_processing(tmp_input,
      p_val_threshold = 5e-2,
      pin_name_path = "Biogrid"
    ),
    "After processing, 1 gene \\(or no genes\\) could be mapped to the PIN"
  )

  expect_error(
    input_processing(tmp_input,
      p_val_threshold = 5e-2,
      pin_name_path = "Biogrid",
      convert2alias = "INVALID"
    ),
    "`convert2alias` should be either TRUE or FALSE"
  )
})

# annotate_term_genes -----------------------------------------------------
example_gene_data <- example_pathfindR_input[1:10, ]
colnames(example_gene_data) <- c("GENE", "CHANGE", "P_VALUE")
tmp_res <- example_pathfindR_output[1:5, -c(7, 8)]

test_that("`annotate_term_genes()` adds input genes for each term", {
  expect_is(
    annotated_result <- annotate_term_genes(
      result_df = tmp_res,
      input_processed = example_gene_data
    ),
    "data.frame"
  )
  expect_true("Up_regulated" %in% colnames(annotated_result) &
    "Down_regulated" %in% colnames(annotated_result))
  expect_true(nrow(annotated_result) == nrow(tmp_res))
})

test_that("annotate_term_genes() argument checks work", {
  expect_error(
    annotate_term_genes(
      result_df = list(),
      input_processed = example_gene_data
    ),
    "`result_df` should be a data frame"
  )
  expect_error(
    annotate_term_genes(
      result_df = tmp_res[, -1],
      input_processed = example_gene_data
    ),
    "`result_df` should contain an \"ID\" column"
  )

  expect_error(
    annotate_term_genes(
      result_df = tmp_res,
      input_processed = list()
    ),
    "`input_processed` should be a data frame"
  )
  expect_error(
    annotate_term_genes(
      result_df = tmp_res,
      input_processed = example_gene_data[, -1]
    ),
    "`input_processed` should contain the columns \"GENE\" and \"CHANGE\""
  )


  expect_error(
    annotate_term_genes(
      result_df = tmp_res,
      input_processed = example_gene_data,
      genes_by_term = "INVALID"
    ),
    "`genes_by_term` should be a list of term gene sets"
  )
  expect_error(
    annotate_term_genes(
      result_df = tmp_res,
      input_processed = example_gene_data,
      genes_by_term = list(1)
    ),
    "`genes_by_term` should be a named list \\(names are gene set IDs\\)"
  )
})
