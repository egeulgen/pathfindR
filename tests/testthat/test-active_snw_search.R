## Tests for functions related to active subnetwork search - Aug 2023

# set up input data
input_data_frame <- example_pathfindR_input[1:10, c(1, 3)]
colnames(input_data_frame) <- c("GENE", "P_VALUE")

example_snws_len <- 20
example_snw_output <- system.file("extdata", "resultActiveSubnetworkSearch.txt",
    package = "pathfindR")
mock_file_path <- function(...) {
    args <- list(...)
    if (args[[1]] == "active_snw_search") {
        return(example_snw_output)
    }
    return(file.path(...))
}

test_that("`active_snw_search()` -- returns a list object", {
    mockery::stub(active_snw_search, "dir.exists", TRUE)
    mockery::stub(active_snw_search, "file.exists", TRUE)
    mockery::stub(active_snw_search, "normalizePath", NULL)
    mockery::stub(active_snw_search, "system", NULL)
    mockery::stub(active_snw_search, "file.path", mock_file_path)
    mockery::stub(active_snw_search, "file.rename", NULL)

    # Expect > 0 active snws
    expect_message(snw_list <- active_snw_search(input_for_search = input_data_frame),
        "Found [1-9]\\d* active subnetworks")
    expect_is(snw_list, "list")
    expect_is(snw_list[[1]], "character")
    expect_true(length(snw_list) > 0)

    # Expect no active snws
    mockery::stub(active_snw_search, "filterActiveSnws", NULL)
    expect_message(snw_list <- active_snw_search(input_for_search = input_data_frame),
        "Found 0 active subnetworks")
    expect_identical(snw_list, list())
})

test_that("`active_snw_search()` -- `dir_for_parallel_run` arg is used when provided",
    {
        mockery::stub(active_snw_search, "dir.exists", TRUE)
        mockery::stub(active_snw_search, "file.exists", TRUE)
        mockery::stub(active_snw_search, "normalizePath", NULL)
        mockery::stub(active_snw_search, "system", NULL)
        mockery::stub(active_snw_search, "file.path", mock_file_path)
        mockery::stub(active_snw_search, "file.rename", NULL)

        m <- mockery::mock(NULL, cycle = TRUE)
        mockery::stub(active_snw_search, "setwd", m)
        res <- active_snw_search(input_for_search = input_data_frame, dir_for_parallel_run = tempdir())
        mockery::expect_called(m, 2)
    })

test_that("`active_snw_search()` -- argument checks work", {
    # input_for_search
    expect_error(snw_list <- active_snw_search(input_for_search = list()), "`input_for_search` should be data frame")

    invalid_input <- input_data_frame
    colnames(invalid_input) <- c("A", "B")
    expect_error(snw_list <- active_snw_search(input_for_search = invalid_input),
        paste0("`input_for_search` should contain the columns ", paste(dQuote(c("GENE",
            "P_VALUE")), collapse = ",")))

    # snws_file
    expect_error(snw_list <- active_snw_search(input_for_search = input_data_frame,
        snws_file = "[/]"), "`snws_file` may be containing forbidden characters. Please change and try again")

    # search_method
    valid_mets <- c("GR", "SA", "GA")
    expect_error(active_snw_search(input_for_search = input_data_frame, search_method = "INVALID"),
        paste0("`search_method` should be one of ", paste(dQuote(valid_mets), collapse = ", ")))

    # silent_option
    expect_error(active_snw_search(input_for_search = input_data_frame, silent_option = "WRONG"),
        "`silent_option` should be either TRUE or FALSE")

    expect_error(active_snw_search(input_for_search = input_data_frame, use_all_positives = "INVALID"),
        "`use_all_positives` should be either TRUE or FALSE")
})

test_that("`active_snw_search()` -- all search methods work", {
    skip_on_cran()
    ## GR
    expect_message(snw_list <- active_snw_search(input_for_search = input_data_frame,
        pin_name_path = "Biogrid", search_method = "GR", dir_for_parallel_run = tempdir(check = TRUE)),
        "Found [1-9]\\d* active subnetworks")
    expect_is(snw_list, "list")
    expect_is(snw_list[[1]], "character")

    skip("will test SA and GA if we can create a suitable (faster and non-empty) test case")
    ## SA
    expect_message(snw_list <- active_snw_search(input_for_search = input_data_frame,
        pin_name_path = "Biogrid", search_method = "SA", dir_for_parallel_run = tempdir(check = TRUE)),
        "Found [1-9]\\d* active subnetworks")
    expect_is(snw_list, "list")
    expect_is(snw_list[[1]], "character")

    ## GA
    expect_message(snw_list <- active_snw_search(input_for_search = input_data_frame,
        pin_name_path = "Biogrid", search_method = "GA", dir_for_parallel_run = tempdir(check = TRUE)),
        "Found [1-9]\\d* active subnetworks")
    expect_is(snw_list, "list")
    expect_is(snw_list[[1]], "character")
})

test_that("`active_snw_search()` -- results are reproducible", {
    skip_on_cran()
    snw_lists <- list()
    seed_vals <- c(123, 123, 456)
    for (idx in 1:3) {
        seed <- seed_vals[idx]
        snw_lists[[idx]] <- active_snw_search(input_for_search = input_data_frame,
            seedForRandom = seed, dir_for_parallel_run = tempdir(check = TRUE))
    }
    expect_identical(snw_lists[[1]], snw_lists[[2]])
    expect_false(identical(snw_lists[[1]], snw_lists[[3]]))
})


test_that("`filterActiveSnws()` -- returns expected list object", {
    snws_filtered <- filterActiveSnws(active_snw_path = example_snw_output, sig_genes_vec = input_data_frame$GENE)
    expect_is(snws_filtered, "list")
    expect_length(snws_filtered, 2)
    expect_is(snws_filtered$subnetworks, "list")
    expect_is(snws_filtered$scores, "numeric")

    expect_is(snws_filtered$subnetworks[[1]], "character")
    expect_true(length(snws_filtered$subnetworks) <= example_snws_len)

    # empty file case
    empty_path <- tempfile("empty", fileext = ".txt")
    file.create(empty_path)
    expect_null(suppressWarnings(filterActiveSnws(active_snw_path = empty_path, sig_genes_vec = input_data_frame$GENE)))
})

test_that("`filterActiveSnws()` -- `score_quan_thr` works", {
    snws_filtered <- filterActiveSnws(active_snw_path = example_snw_output, sig_genes_vec = example_pathfindR_input$Gene.symbol,
        score_quan_thr = -1, sig_gene_thr = 0)
    expect_length(snws_filtered$subnetworks, example_snws_len)

    for (q_thr in seq(0.1, 1, by = 0.1)) {
        snws_filtered <- filterActiveSnws(active_snw_path = example_snw_output, sig_genes_vec = example_pathfindR_input$Gene.symbol,
            score_quan_thr = q_thr, sig_gene_thr = 0)
        exp_len <- example_snws_len * (1 - q_thr)
        expect_length(snws_filtered$subnetworks, as.integer(exp_len + 0.5))
    }
})

test_that("`filterActiveSnws()` -- `sig_gene_thr` works", {
    snws_filtered1 <- filterActiveSnws(active_snw_path = example_snw_output, sig_genes_vec = example_pathfindR_input$Gene.symbol,
        sig_gene_thr = 0.02, score_quan_thr = -1)
    snws_filtered2 <- filterActiveSnws(active_snw_path = example_snw_output, sig_genes_vec = example_pathfindR_input$Gene.symbol,
        sig_gene_thr = 0.1, score_quan_thr = -1)

    expect_true(length(snws_filtered2$subnetworks) < example_snws_len)
    expect_true(length(snws_filtered1$subnetworks) > length(snws_filtered2$subnetworks))
})

test_that("`filterActiveSnws()` -- argument checks work", {
    expect_error(filterActiveSnws(active_snw_path = "this/is/not/a/valid/path"),
        "The active subnetwork file does not exist! Check the `active_snw_path` argument")

    expect_error(filterActiveSnws(active_snw_path = example_snw_output, sig_genes_vec = list()),
        "`sig_genes_vec` should be a vector")

    expect_error(filterActiveSnws(active_snw_path = example_snw_output, sig_genes_vec = example_pathfindR_input$Gene.symbol,
        score_quan_thr = "INVALID"), "`score_quan_thr` should be numeric")
    expect_error(filterActiveSnws(active_snw_path = example_snw_output, sig_genes_vec = example_pathfindR_input$Gene.symbol,
        score_quan_thr = -2), "`score_quan_thr` should be in \\[0, 1\\] or -1 \\(if not filtering\\)")
    expect_error(filterActiveSnws(active_snw_path = example_snw_output, sig_genes_vec = example_pathfindR_input$Gene.symbol,
        score_quan_thr = 2), "`score_quan_thr` should be in \\[0, 1\\] or -1 \\(if not filtering\\)")

    expect_error(filterActiveSnws(active_snw_path = example_snw_output, sig_genes_vec = example_pathfindR_input$Gene.symbol,
        sig_gene_thr = "INVALID"), "`sig_gene_thr` should be numeric")
    expect_error(filterActiveSnws(active_snw_path = example_snw_output, sig_genes_vec = example_pathfindR_input$Gene.symbol,
        sig_gene_thr = -1), "`sig_gene_thr` should be in \\[0, 1\\]")
})

test_that("`visualize_active_subnetworks()` -- returns list of ggraph objects", {
    # empty file case
    empty_path <- tempfile("empty", fileext = ".txt")
    file.create(empty_path)
    expect_null(visualize_active_subnetworks(active_snw_path = empty_path, genes_df = input_data_frame))

    skip_on_cran()
    # default
    g_list <- visualize_active_subnetworks(example_snw_output, input_data_frame)
    expect_is(g_list, "list")
    expect_is(g_list[[1]], "ggraph")
    expect_true(length(g_list) <= example_snws_len)

    # set `num_snws` to larger than actual number
    g_list <- visualize_active_subnetworks(example_snw_output, input_data_frame,
        num_snws = 21)
    expect_is(g_list, "list")
    expect_is(g_list[[1]], "ggraph")
    expect_length(g_list, 3)
})
