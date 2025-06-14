## Tests for functions related to data generation - June 2025

set.seed(123)
gene_pool <- paste0("Gene", 1:100)
toy_biogrid_pin <- data.frame(A = sample(gene_pool, 25), B = sample(gene_pool, 25))
colnames(toy_biogrid_pin) <- c("Official Symbol Interactor A", "Official Symbol Interactor B")

test_that("`process_pin()` -- removes self-interactions and duplicated interactions",
    {
        input_pin_df <- toy_biogrid_pin
        colnames(input_pin_df) <- c("Interactor_A", "Interactor_B")
        input_pin_df <- rbind(input_pin_df, data.frame(Interactor_A = input_pin_df$Interactor_B[1:5],
            Interactor_B = input_pin_df$Interactor_A[1:5]))

        processed_df <- process_pin(input_pin_df)

        expect_true(nrow(processed_df) < nrow(input_pin_df))
    })

test_that("`get_biogrid_pin()` -- returns a path to a valid PIN file", {
    mockery::stub(get_biogrid_pin, "utils::download.file", NULL)
    mockery::stub(get_biogrid_pin, "utils::unzip", list(Name = "BIOGRID-ORGANISM-Homo_sapiens-4.4.211.tab3.txt"))
    mockery::stub(get_biogrid_pin, "utils::read.delim", toy_biogrid_pin)

    expected_biogrid_pin_df <- toy_biogrid_pin
    colnames(expected_biogrid_pin_df) <- c("Interactor_A", "Interactor_B")
    expected_biogrid_pin_df <- process_pin(expected_biogrid_pin_df)
    expected_biogrid_pin_df <- data.frame(V1 = expected_biogrid_pin_df$Interactor_A,
        V2 = "pp", V3 = expected_biogrid_pin_df$Interactor_B)

    pin_path <- get_biogrid_pin(release = "4.4.211")
    pin_df <- read.delim(pin_path, header = FALSE)
    expect_true(ncol(pin_df) == 3)
    expect_true(all(pin_df[, 2] == "pp"))
    expect_identical(pin_df, expected_biogrid_pin_df)
})

test_that("`get_biogrid_pin()` -- determines and downloads the latest version", {
  mockery::stub(get_biogrid_pin, "utils::download.file", NULL)
  mockery::stub(get_biogrid_pin, "utils::unzip", list(Name = "BIOGRID-ORGANISM-Homo_sapiens-X.X.X.tab3.txt"))
  mockery::stub(get_biogrid_pin, "utils::read.delim", toy_biogrid_pin)

  expected_biogrid_pin_df <- toy_biogrid_pin
  colnames(expected_biogrid_pin_df) <- c("Interactor_A", "Interactor_B")
  expected_biogrid_pin_df <- process_pin(expected_biogrid_pin_df)
  expected_biogrid_pin_df <- data.frame(V1 = expected_biogrid_pin_df$Interactor_A,
                                        V2 = "pp", V3 = expected_biogrid_pin_df$Interactor_B)

  pin_path <- get_biogrid_pin()
  pin_df <- read.delim(pin_path, header = FALSE)
  expect_true(ncol(pin_df) == 3)
  expect_true(all(pin_df[, 2] == "pp"))
  expect_identical(pin_df, expected_biogrid_pin_df)
})

test_that("`get_biogrid_pin()` -- error check works", {
    # invalid organism error
    expect_error(get_biogrid_pin(org = "Hsapiens"), paste("Hsapiens is not a valid Biogrid organism.",
        "Available organisms are listed on: https://wiki.thebiogrid.org/doku.php/statistics"))
})

test_that("`get_pin_file()` -- works as expected", {
    with_mocked_bindings({
        pin_path <- get_pin_file()
        expect_identical(pin_path, "/path/to/some/PIN/file")
    }, get_biogrid_pin = function(...) "/path/to/some/PIN/file", .package = "pathfindR")

    expect_error(get_pin_file(source = "STRING"), "As of this version, this function is implemented to get data from BioGRID only")
})

test_that("`gset_list_from_gmt()` -- works as expected", {
    gmt_list <- list(GSA = sample(gene_pool, 80), GSB = sample(gene_pool, 100), GSC = sample(gene_pool,
        33))
    description_vec <- c(GSA = "gene set A", GSB = "gene set B", GSC = "gene set C")

    gmt_df <- c()
    for (gset in names(gmt_list)) {
        tmp <- c(gset, description_vec[gset])
        tmp <- c(tmp, gmt_list[[gset]], rep("", 100 - length(gmt_list[[gset]])))
        gmt_df <- rbind(gmt_df, tmp)
    }

    path2gmt <- tempfile()
    write.table(gmt_df, path2gmt, sep = "\t", col.names = FALSE, row.names = FALSE,
        quote = FALSE)

    expect_is(res <- gset_list_from_gmt(path2gmt), "list")
    expect_identical(res$gene_sets, gmt_list)
    expect_identical(res$descriptions, description_vec)
})


test_that("`get_kegg_gsets()` -- works as expected", {
  skip_on_cran()
  mock_response <- "eco00010\tdescription\neco00071\tdescription2"

  # mocked binding to manage sequential responses
  with_mocked_bindings(
    {
      expect_is(toy_eco_kegg <- pathfindR:::get_kegg_gsets(), "list")
    },
    content = function(...) mock_response, .package = "httr"
  )

  expect_length(toy_eco_kegg, 2)
  expect_true(all(names(toy_eco_kegg) == c("gene_sets", "descriptions")))
  expect_true(all(names(toy_eco_kegg[["gene_sets"]]) %in% names(toy_eco_kegg[["descriptions"]])))
  expect_length(toy_eco_kegg[["gene_sets"]], 2)
  expect_length(toy_eco_kegg[["descriptions"]], 2)

  expect_true(toy_eco_kegg[["descriptions"]]["eco00010"] == "description")
  expect_true(toy_eco_kegg[["descriptions"]]["eco00071"] == "description2")

  expect_length(toy_eco_kegg[["gene_sets"]][["eco00010"]], 47)
  expect_length(toy_eco_kegg[["gene_sets"]][["eco00071"]], 15)
})

test_that("`get_reactome_gsets()` -- works as expected", {
  skip_on_cran()

  pw1 <- "Pathway1"
  pw2 <- "Pathway2"
  desc1 <- "Description1"
  desc2 <- "Description2"
  genes1 <- c("GeneA", "GeneB")
  genes2 <- c("GeneC", "GeneD", "GeneE")

  gmt_content <- paste(
    c(
      paste(c(desc1, pw1, genes1), collapse = "\t"),
      paste(c(desc2, pw2, genes2), collapse = "\t")
    ),
    collapse = "\n"
  )

  mockery::stub(get_reactome_gsets, "utils::download.file", NULL)

  unz_mock <- function(zipfile, filename, ...) {
    textConnection(gmt_content)
  }
  mockery::stub(get_reactome_gsets, "unz", unz_mock)

  expected_gsets <- list(genes1, genes2)
  names(expected_gsets) <- c(pw1, pw2)
  expect_descriptions <- c(desc1, desc2)
  names(expect_descriptions) <- c(pw1, pw2)

  expect_is(reactome <- get_reactome_gsets(), "list")
  expect_length(reactome, 2)
  expect_length(reactome$gene_sets, 2)
  expect_length(reactome$descriptions, 2)
  expect_equal(names(reactome$gene_sets), names(reactome$descriptions))
  expect_equal(reactome$gene_sets, expected_gsets)
  expect_equal(reactome$descriptions, expect_descriptions)
})

test_that("`get_mgsigdb_gsets()` -- works as expected", {
    toy_msigdb_df <- c()
    for (gs_idx in 1:5) {
        toy_msigdb_df <- rbind(toy_msigdb_df, data.frame(gene_symbol = sample(gene_pool,
            sample(25:75, 1)), gs_id = paste0("GS", gs_idx), gs_name = paste("Gene Set",
            gs_idx)))
    }
    mockery::stub(get_mgsigdb_gsets, "msigdbr::msigdbr", toy_msigdb_df)

    expect_is(res_msig_db <- get_mgsigdb_gsets(collection = "C1"), "list")
    expect_length(res_msig_db, 2)
    expect_true(all(names(res_msig_db) == c("gene_sets", "descriptions")))
    expect_true(all(names(res_msig_db[["gene_sets"]] %in% names(res_msig_db[["descriptions"]]))))
})

test_that("`get_gene_sets_list()` works", {
    expect_error(gsets <- get_gene_sets_list("Wiki"), "As of this version, this function is implemented to get data from KEGG, Reactome and MSigDB only")

    mockery::stub(get_gene_sets_list, "get_kegg_gsets", NULL)
    mockery::stub(get_gene_sets_list, "get_reactome_gsets", NULL)
    mockery::stub(get_gene_sets_list, "get_mgsigdb_gsets", NULL)
    expect_silent(kegg <- get_gene_sets_list(org_code = "vcn"))
    expect_message(rctm <- get_gene_sets_list("Reactome"))
    expect_silent(msig <- get_gene_sets_list("MSigDB", species = "Mus musculus", db_species = "MS",
        collection = "C3", subcollection = "MIR:MIR_Legacy"))
})
