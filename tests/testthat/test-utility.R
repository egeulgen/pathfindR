## Tests for various utility functions - Aug 2023

set.seed(123)

test_that("`active_snw_enrichment_wrapper()` -- works as expected", {
    input_df <- example_pathfindR_input[, c(1, 3)]
    colnames(input_df) <- c("GENE", "P_VALUE")

    org_dir <- getwd()
    test_directory <- file.path(tempdir(check = TRUE), "snw_wrapper_test")
    dir.create(test_directory)
    setwd(test_directory)
    on.exit(setwd(org_dir))
    on.exit(unlink(test_directory), add = TRUE)

    with_mocked_bindings({
        expect_is(active_snw_enrichment_wrapper(input_processed = input_df, pin_path = "Biogrid",
            gset_list = list(), enrichment_threshold = 0.05, list_active_snw_genes = FALSE,
            iterations = 1), "data.frame")

        expect_is(active_snw_enrichment_wrapper(input_processed = input_df, pin_path = "Biogrid",
            gset_list = list(), enrichment_threshold = 0.05, list_active_snw_genes = FALSE,
            iterations = 2, disable_parallel = TRUE), "data.frame")

        expect_warning(active_snw_enrichment_wrapper(input_processed = input_df,
            pin_path = "Biogrid", gset_list = list(), enrichment_threshold = 0.05,
            list_active_snw_genes = FALSE, search_method = "GA", iterations = 2))
    }, single_iter_wrapper = function(...) example_pathfindR_output, .package = "pathfindR")

    skip_on_cran()
    expect_is(active_snw_enrichment_wrapper(input_processed = input_df[1:10, ], pin_path = "Biogrid",
        gset_list = list(genes_by_term = kegg_genes[1:2], term_descriptions = kegg_descriptions[names(kegg_genes[1:2])]),
        enrichment_threshold = 0.05, list_active_snw_genes = FALSE, iterations = 2),
        "NULL")
})

test_that("`active_snw_enrichment_wrapper()` -- argument checks work", {
    valid_mets <- c("GR", "SA", "GA")
    expect_error(active_snw_enrichment_wrapper(input_processed = input_processed,
        pin_path = pin_path, gset_list = list(), enrichment_threshold = 0.05, list_active_snw_genes = FALSE,
        search_method = "INVALID"), paste0("`search_method` should be one of ", paste(dQuote(valid_mets),
        collapse = ", ")))

    expect_error(active_snw_enrichment_wrapper(input_processed = input_processed,
        pin_path = pin_path, gset_list = list(), enrichment_threshold = 0.05, list_active_snw_genes = FALSE,
        use_all_positives = "INVALID"), "`use_all_positives` should be either TRUE or FALSE")

    expect_error(active_snw_enrichment_wrapper(input_processed = input_processed,
        pin_path = pin_path, gset_list = list(), enrichment_threshold = 0.05, list_active_snw_genes = FALSE,
        silent_option = "INVALID"), "`silent_option` should be either TRUE or FALSE")

    expect_error(active_snw_enrichment_wrapper(input_processed = input_processed,
        pin_path = pin_path, gset_list = list(), enrichment_threshold = 0.05, list_active_snw_genes = FALSE,
        disable_parallel = "INVALID"), "`disable_parallel` should be either TRUE or FALSE")

    expect_error(active_snw_enrichment_wrapper(input_processed = input_processed,
        pin_path = pin_path, gset_list = list(), enrichment_threshold = 0.05, list_active_snw_genes = FALSE,
        iterations = "INVALID"), "`iterations` should be a positive integer")

    expect_error(active_snw_enrichment_wrapper(input_processed = input_processed,
        pin_path = pin_path, gset_list = list(), enrichment_threshold = 0.05, list_active_snw_genes = FALSE,
        iterations = 0), "`iterations` should be >= 1")

    expect_error(active_snw_enrichment_wrapper(input_processed = input_processed,
        pin_path = pin_path, gset_list = list(), enrichment_threshold = 0.05, list_active_snw_genes = FALSE,
        n_processes = "INVALID"), "`n_processes` should be either NULL or a positive integer")

    expect_error(active_snw_enrichment_wrapper(input_processed = input_processed,
        pin_path = pin_path, gset_list = list(), enrichment_threshold = 0.05, list_active_snw_genes = FALSE,
        n_processes = 0), "`n_processes` should be > 1")
})

test_that("`configure_output_dir()` -- works as expected", {
    expected_dir <- file.path(tempdir(), "test_pathfindR_results")
    mockery::stub(configure_output_dir, "file.path", expected_dir)
    expect_equal(configure_output_dir(), expected_dir)

    test_out_dir <- file.path(tempdir(), "TEST")
    for (i in 1:3) {
        actual_dir <- configure_output_dir(test_out_dir)
        dir_to_check <- test_out_dir
        if (i > 1) {
            dir_to_check <- paste0(dir_to_check, "(", i - 1, ")")
        }
        expect_equal(actual_dir, dir_to_check)
        dir.create(actual_dir)
    }
})

test_that("`fetch_gene_set()` -- can fetch all gene set objects", {
    skip_on_cran()
    for (gset_name in c("KEGG", "mmu_KEGG", "Reactome", "BioCarta", "cell_markers",
        "GO-All", "GO-BP", "GO-CC", "GO-MF")) {
        expect_is(gset_obj <- fetch_gene_set(gene_sets = gset_name, min_gset_size = 10,
            max_gset_size = 300), "list")
        expect_is(gset_obj$genes_by_term, "list")
        expect_is(gset_obj$term_descriptions, "character")
        expect_true(length(gset_obj$genes_by_term) == length(gset_obj$term_descriptions))
        tmp <- vapply(gset_obj$genes_by_term, length, 1L)
        expect_true(min(tmp) >= 10 & max(tmp) <= 300)
    }
    # Custom
    gset_obj <- fetch_gene_set(gene_sets = "Custom", min_gset_size = 20, max_gset_size = 200,
        custom_genes = kegg_genes, custom_descriptions = kegg_descriptions)
    expect_is(gset_obj$genes_by_term, "list")
    expect_is(gset_obj$term_descriptions, "character")
    expect_true(length(gset_obj$genes_by_term) == length(gset_obj$term_descriptions))
    tmp <- vapply(gset_obj$genes_by_term, length, 1L)
    expect_true(min(tmp) >= 20 & max(tmp) <= 200)
})

test_that("`create_HTML_report()` -- works a expected", {
    mock_render <- mockery::mock(NULL, cycle = TRUE)
    mockery::stub(create_HTML_report, "rmarkdown::render", mock_render)

    create_HTML_report(input = data.frame(), input_processed = data.frame(), final_res = data.frame(),
        dir_for_report = "/path/to/report/dir")
    mockery::expect_called(mock_render, 3)
})

test_that("`fetch_gene_set()` -- min/max_gset_size args correctly filter gene sets",
    {
        skip_on_cran()
        min_max_pairs <- list(c(min = 10, max = 300), c(min = 50, max = 200))
        num_of_terms_after_size_filtering <- c()
        for (idx in seq_along(min_max_pairs)) {
            cur_vals <- min_max_pairs[[idx]]
            expect_is(gset_obj <- fetch_gene_set(gene_sets = "KEGG", min_gset_size = cur_vals["min"],
                max_gset_size = cur_vals["max"]), "list")
            sizes_of_terms <- vapply(gset_obj$genes_by_term, length, 1L)
            expect_true(min(sizes_of_terms) >= cur_vals["min"] & max(sizes_of_terms) <=
                cur_vals["max"])
            num_of_terms_after_size_filtering <- c(num_of_terms_after_size_filtering,
                length(gset_obj$genes_by_term))
        }

        expect_true(num_of_terms_after_size_filtering[2] < num_of_terms_after_size_filtering[1])
    })

test_that("`fetch_gene_set()` -- for 'Custom' gene set, check if the custom objects are provided",
    {
        expect_error(fetch_gene_set(gene_sets = "Custom"), "`custom_genes` and `custom_descriptions` must be provided if `gene_sets = \"Custom\"`")
        expect_error(fetch_gene_set(gene_sets = "Custom", custom_genes = kegg_genes),
            "`custom_genes` and `custom_descriptions` must be provided if `gene_sets = \"Custom\"`")
        expect_error(fetch_gene_set(gene_sets = "Custom", custom_descriptions = kegg_descriptions),
            "`custom_genes` and `custom_descriptions` must be provided if `gene_sets = \"Custom\"`")
    })

test_that("`fetch_gene_set()` -- argument checks work", {
    all_gs_opts <- c("KEGG", "Reactome", "BioCarta", "GO-All", "GO-BP", "GO-CC",
        "GO-MF", "cell_markers", "mmu_KEGG", "Custom")
    expect_error(fetch_gene_set(gene_sets = "INVALID"), paste0("`gene_sets` should be one of ",
        paste(dQuote(all_gs_opts), collapse = ", ")))

    expect_error(fetch_gene_set(min_gset_size = "INVALID"), "`min_gset_size` should be numeric")

    expect_error(fetch_gene_set(max_gset_size = "INVALID"), "`max_gset_size` should be numeric")

    expect_error(fetch_gene_set(gene_sets = "Custom", custom_genes = "INVALID", custom_descriptions = ""),
        "`custom_genes` should be a list of term gene sets")
    expect_error(fetch_gene_set(gene_sets = "Custom", custom_genes = list(), custom_descriptions = ""),
        "`custom_genes` should be a named list \\(names are gene set IDs\\)")

    expect_error(fetch_gene_set(gene_sets = "Custom", custom_genes = kegg_genes,
        custom_descriptions = list()), "`custom_descriptions` should be a vector of term gene descriptions")
    expect_error(fetch_gene_set(gene_sets = "Custom", custom_genes = kegg_genes,
        custom_descriptions = 1:3), "`custom_descriptions` should be a named vector \\(names are gene set IDs\\)")
})

test_that("`return_pin_path()` -- returns the absolute path to PIN file", {
    mockery::stub(return_pin_path, "utils::getFromNamespace", list())
    mockery::stub(return_pin_path, "lapply", list(data.frame(V1 = paste0("G", 1:10),
        V2 = "pp", V3 = paste0("G", 2:11)), data.frame(V1 = paste0("G", 3:5), V2 = "pp",
        V3 = paste0("G", 5:7))))
    expect_silent(path2file <- return_pin_path("Biogrid"))
    expect_true(file.exists(path2file))

    custom_pin <- read.delim(path2file, header = FALSE)
    custom_pin$V1 <- tolower(custom_pin$V1)
    custom_sif_path <- file.path(tempdir(check = TRUE), "tmp_PIN.sif")
    utils::write.table(custom_pin, custom_sif_path, sep = "\t", row.names = FALSE,
        col.names = FALSE, quote = FALSE)
    expect_silent(final_custom_path <- return_pin_path(custom_sif_path))
    expect_true(file.exists(final_custom_path))

    # convert to uppercase works
    upper_case_custom <- read.delim(final_custom_path, header = FALSE)
    expect_true(all(toupper(upper_case_custom[, 1]) == upper_case_custom[, 1]))
    expect_true(all(toupper(upper_case_custom[, 3]) == upper_case_custom[, 3]))


    # invalid custom PIN - wrong format
    invalid_sif_path <- system.file(paste0("extdata/MYC.txt"), package = "pathfindR")
    expect_error(return_pin_path(invalid_sif_path), "The PIN file must have 3 columns and be tab-separated")

    # invalid custom PIN - invalid second column
    invalid_sif_path <- file.path(tempdir(check = TRUE), "custom.sif")
    invalid_custom_sif <- data.frame(P1 = "X", pp = "INVALID", P2 = "Y")
    write.table(invalid_custom_sif, invalid_sif_path, sep = "\t", col.names = FALSE,
        row.names = FALSE)
    expect_error(return_pin_path(invalid_sif_path), "The second column of the PIN file must all be \"pp\" ")

    # invalid option
    valid_opts <- c("Biogrid", "STRING", "GeneMania", "IntAct", "KEGG", "mmu_STRING",
        "/path/to/custom/SIF")
    expect_error(return_pin_path("INVALID"), paste0("The chosen PIN must be one of:\n",
        paste(dQuote(valid_opts), collapse = ", ")))
})

test_that("`input_testing()` -- works as expected", {
    expect_message(input_testing(input = example_pathfindR_input, p_val_threshold = 0.05),
        "The input looks OK")

    expect_error(input_testing(input = matrix(), p_val_threshold = 0.05), "the input is not a data frame")

    expect_error(input_testing(input = example_pathfindR_input[, 1, drop = FALSE],
        p_val_threshold = 0.05), "the input should have 2 or 3 columns")

    expect_error(input_testing(input = example_pathfindR_input[1, ], p_val_threshold = 0.05),
        "There must be at least 2 rows \\(genes\\) in the input data frame")

    expect_error(input_testing(input = example_pathfindR_input, p_val_threshold = "INVALID"),
        "`p_val_threshold` must be a numeric value between 0 and 1")

    expect_error(input_testing(input = example_pathfindR_input, p_val_threshold = -1),
        "`p_val_threshold` must be between 0 and 1")

    tmp <- example_pathfindR_input
    tmp$adj.P.Val <- NA
    expect_error(input_testing(input = tmp, p_val_threshold = 0.05), "p values cannot contain NA values")

    tmp <- example_pathfindR_input
    tmp$adj.P.Val <- "INVALID"
    expect_error(input_testing(input = tmp, p_val_threshold = 0.05), "p values must all be numeric")

    tmp <- example_pathfindR_input
    tmp$adj.P.Val[1] <- -1
    expect_error(input_testing(input = tmp, p_val_threshold = 0.05), "p values must all be between 0 and 1")
})

test_that("`input_processing()` -- works as expected", {
    input_df <- example_pathfindR_input[1:10, ]
    toy_PIN <- data.frame(V1 = sample(example_pathfindR_input$Gene.symbol, 100),
        V2 = "pp", V3 = sample(example_pathfindR_input$Gene.symbol, 100))
    mockery::stub(input_processing, "return_pin_path", NULL)
    mockery::stub(input_processing, "utils::read.delim", toy_PIN)

    expect_is(processed_df <- input_processing(input_df), "data.frame")
    expect_true(ncol(processed_df) == 4)
    expect_true(nrow(processed_df) <= nrow(example_pathfindR_input))

    # no change values provided
    input_df2 <- input_df[, -2]
    expect_is(processed_df2 <- suppressWarnings(input_processing(input_df2)), "data.frame")
    expect_true(ncol(processed_df2) == 4)
    expect_true(all(processed_df2$CHANGE == 1e+06))

    toy_PIN2 <- rbind(toy_PIN, data.frame(V1 = c("SERPINA3", "ARHGAP17"), V2 = "pp",
        V3 = c("ACT", "GIG25")))
    mockery::stub(input_processing, "utils::read.delim", toy_PIN2)

    # multiple mapping
    input_multimap <- input_df
    input_multimap$Gene.symbol[1] <- "GIG24"
    input_multimap$Gene.symbol[2] <- "ACT"
    input_multimap$Gene.symbol[3] <- "AACT"
    input_multimap$Gene.symbol[4] <- "GIG25"
    expect_is(processed_df3 <- input_processing(input_multimap), "data.frame")
})

test_that("`input_processing()` -- errors and warnings work", {
    input_df <- example_pathfindR_input[1:10, ]

    toy_PIN <- data.frame(V1 = sample(input_df$Gene.symbol, 7), V2 = "pp", V3 = sample(input_df$Gene.symbol,
        7))
    mockery::stub(input_processing, "return_pin_path", NULL)
    mockery::stub(input_processing, "utils::read.delim", toy_PIN)

    input_df$Gene.symbol <- as.factor(input_df$Gene.symbol)
    expect_warning(input_processing(input_df, p_val_threshold = 0.05, pin_name_path = "Biogrid",
        convert2alias = TRUE), "The gene column was turned into character from factor.")

    expect_error(input_processing(example_pathfindR_input, p_val_threshold = 1e-100,
        pin_name_path = "Biogrid"), "No input p value is lower than the provided threshold \\(1e-100\\)")

    input_dup <- example_pathfindR_input[1:3, ]
    input_dup <- rbind(input_dup, input_dup[1, ])
    expect_warning(input_processing(input_dup, p_val_threshold = 0.05, pin_name_path = "Biogrid"),
        "Duplicated genes found! The lowest p value for each gene was selected")

    low_sig_input <- example_pathfindR_input[1:3, ]
    low_sig_input$adj.P.Val <- 1e-15
    expect_message(res <- input_processing(low_sig_input, p_val_threshold = 0.05,
        pin_name_path = "Biogrid"), "pathfindR cannot handle p values < 1e-13. These were changed to 1e-13")
    expect_true(all(res$P_VALUE == 1e-13))

    invalid_genes_input <- low_sig_input
    invalid_genes_input$Gene.symbol <- paste0(LETTERS[seq_len(nrow(invalid_genes_input))],
        "INVALID")
    expect_error(input_processing(invalid_genes_input, p_val_threshold = 0.05, pin_name_path = "Biogrid"),
        "None of the genes were in the PIN\nPlease check your gene symbols")

    low_sig_input$Gene.symbol[1] <- "INVALID_A"
    low_sig_input$Gene.symbol[2] <- "INVALID_B"
    low_sig_input$Gene.symbol[3] <- toy_PIN$V1[1]
    expect_error(input_processing(low_sig_input, p_val_threshold = 0.05, pin_name_path = "Biogrid"),
        "After processing, 1 gene \\(or no genes\\) could be mapped to the PIN")

    expect_error(input_processing(low_sig_input, p_val_threshold = 0.05, pin_name_path = "Biogrid",
        convert2alias = "INVALID"), "`convert2alias` should be either TRUE or FALSE")
})

example_gene_data <- example_pathfindR_input[1:10, ]
colnames(example_gene_data) <- c("GENE", "CHANGE", "P_VALUE")
tmp_res <- example_pathfindR_output[1:5, -c(7, 8)]

test_that("`annotate_term_genes()` -- adds input genes for each term", {
    expect_is(annotated_result <- annotate_term_genes(result_df = tmp_res, input_processed = example_gene_data),
        "data.frame")
    expect_true("Up_regulated" %in% colnames(annotated_result) & "Down_regulated" %in%
        colnames(annotated_result))
    expect_true(nrow(annotated_result) == nrow(tmp_res))
})

test_that("annotate_term_genes() -- argument checks work", {
    expect_error(annotate_term_genes(result_df = list(), input_processed = example_gene_data),
        "`result_df` should be a data frame")
    expect_error(annotate_term_genes(result_df = tmp_res[, -1], input_processed = example_gene_data),
        "`result_df` should contain an \"ID\" column")

    expect_error(annotate_term_genes(result_df = tmp_res, input_processed = list()),
        "`input_processed` should be a data frame")
    expect_error(annotate_term_genes(result_df = tmp_res, input_processed = example_gene_data[,
        -1]), "`input_processed` should contain the columns \"GENE\" and \"CHANGE\"")


    expect_error(annotate_term_genes(result_df = tmp_res, input_processed = example_gene_data,
        genes_by_term = "INVALID"), "`genes_by_term` should be a list of term gene sets")
    expect_error(annotate_term_genes(result_df = tmp_res, input_processed = example_gene_data,
        genes_by_term = list(1)), "`genes_by_term` should be a named list \\(names are gene set IDs\\)")
})
