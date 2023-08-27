## Tests for core function - Aug 2023

# set up input data
input_data_frame <- example_pathfindR_input[1:10, c(1, 3)]
colnames(input_data_frame) <- c("GENE", "P_VALUE")

test_that("`run_pathfindR()` -- works as expected", {
    mock_fetch_gene_set <- mockery::mock(list(), cycle = TRUE)
    mock_return_pin_path <- mockery::mock("/path/to/some/PIN/SIF", cycle = TRUE)
    mock_input_processing <- mockery::mock(input_data_frame, cycle = TRUE)
    mock_active_snw_enrichment_wrapper <- mockery::mock(data.frame(), c())
    mock_summarize_enrichment_results <- mockery::mock(data.frame())
    mock_annotate_term_genes <- mockery::mock(example_pathfindR_output)
    mock_plot <- mockery::mock(NULL)

    mockery::stub(run_pathfindR, "fetch_gene_set", mock_fetch_gene_set)
    mockery::stub(run_pathfindR, "return_pin_path", mock_return_pin_path)
    mockery::stub(run_pathfindR, "input_processing", mock_input_processing)
    mockery::stub(run_pathfindR, "active_snw_enrichment_wrapper", mock_active_snw_enrichment_wrapper)
    mockery::stub(run_pathfindR, "summarize_enrichment_results", mock_summarize_enrichment_results)
    mockery::stub(run_pathfindR, "annotate_term_genes", mock_annotate_term_genes)
    mockery::stub(run_pathfindR, "graphics::plot", mock_plot)
    mockery::stub(run_pathfindR, "create_HTML_report", NULL)

    expected_messages <- paste(c("The input looks OK", "Plotting the enrichment bubble chart",
        paste(c(paste0("Found ", nrow(example_pathfindR_output), " enriched terms\n"),
            "You may run:", "- cluster_enriched_terms() for clustering enriched terms",
            "- visualize_terms() for visualizing enriched term diagrams\n"), collapse = "\n")),
        collapse = "|")
    # wrapper functions correctly - with output_dir provided
    out_dir <- file.path(tempdir(check = TRUE), "core_test")
    expect_message(res <- run_pathfindR(input_data_frame, output_dir = out_dir),
        expected_messages)
    expect_is(res, "data.frame")
    expect_identical(res, example_pathfindR_output)
    expect_true(dir.exists(out_dir))
    mockery::expect_called(mock_fetch_gene_set, 1)
    mockery::expect_called(mock_return_pin_path, 1)
    mockery::expect_called(mock_input_processing, 1)
    mockery::expect_called(mock_active_snw_enrichment_wrapper, 1)
    mockery::expect_called(mock_summarize_enrichment_results, 1)
    mockery::expect_called(mock_annotate_term_genes, 1)
    mockery::expect_called(mock_plot, 1)

    # warning raised as expected when no results found
    expect_warning(res <- run_pathfindR(input_data_frame), "Did not find any enriched terms!")
    expect_identical(res, data.frame())
})

test_that("`run_pathfindR()` argument checks work", {
    expect_error(run_pathfindR(input_data_frame, plot_enrichment_chart = "INVALID"),
        "`plot_enrichment_chart` should be either TRUE or FALSE")
    expect_error(run_pathfindR(input_data_frame, list_active_snw_genes = "INVALID"),
        "`list_active_snw_genes` should be either TRUE or FALSE")
})
