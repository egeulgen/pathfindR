##################################################
## Project: pathfindR
## Script purpose: Testthat testing script for
## data generation functions
## Date: Dec 2, 2021
## Author: Ege Ulgen
##################################################

org_met <- getOption("download.file.method")
org_extra <- getOption("download.file.extra")
options(download.file.method="curl", download.file.extra="-k -L")

# get_biogrid_pin ---------------------------------------------------------
test_that("`get_biogrid_pin()` returns a path to a valid PIN file", {
  skip_on_cran()
  # v4.4.198
  pin_path <- pathfindR:::get_biogrid_pin(org="Pan_troglodytes")
  pin_df <- read.delim(pin_path, header = FALSE)
  expect_true(ncol(pin_df) == 3)
  expect_true(all(pin_df[, 2] == "pp"))
})

test_that("`get_biogrid_pin()` error check works", {
  # invalid organism error
  expect_error(pathfindR:::get_biogrid_pin(org = "Hsapiens"),
               paste("Hsapiens is not a valid Biogrid organism.",
                     "Available organisms are listed on: https://wiki.thebiogrid.org/doku.php/statistics"))
})

# get_pin_file ------------------------------------------------------------
test_that("`get_pin_file()` works", {
  # unimplemented error
  expect_error(get_pin_file(source = "STRING"),
"As of this version, this function is implemented to get data from BioGRID only")

  skip_on_cran()
  pin_path <- get_pin_file(source = "BioGRID",
                           org = "Pan_troglodytes",
                           release = "3.5.179")
  pin_df <- read.delim(pin_path, header = FALSE)
  expect_true(ncol(pin_df) == 3)
  expect_true(all(pin_df[, 2] == "pp"))
})

options(download.file.method = org_met, download.file.extra=org_extra)
# get_kegg_gsets ----------------------------------------------------------
test_that("`get_kegg_gsets() works`", {
  skip_on_cran()
  # hsa - default
  expect_silent(hsa_kegg <- pathfindR:::get_kegg_gsets())
  expect_length(hsa_kegg, 2)
  expect_true(all(names(hsa_kegg) == c("gene_sets", "descriptions")))
  expect_true(all(names(hsa_kegg[["gene_sets"]] %in% names(hsa_kegg[["descriptions"]]))))
})

options(download.file.method="curl", download.file.extra="-k -L")
# get_reactome_gsets ------------------------------------------------------
test_that("`get_reactome_gsets()` works", {
  skip_on_cran()
  expect_silent(reactome <- pathfindR:::get_reactome_gsets())
  expect_length(reactome, 2)
  expect_true(all(names(reactome) == c("gene_sets", "descriptions")))
  expect_true(all(names(reactome[["gene_sets"]] %in% names(reactome[["descriptions"]]))))
})

# get_mgsigdb_gsets -------------------------------------------------------
test_that("`get_mgsigdb_gsets()` works", {
  skip_on_cran()
  expect_silent(hsa_C2_cgp <- pathfindR:::get_mgsigdb_gsets(collection = "C3",
                                                            subcollection = "MIR:MIR_Legacy"))
  expect_length(hsa_C2_cgp, 2)
  expect_true(all(names(hsa_C2_cgp) == c("gene_sets", "descriptions")))
  expect_true(all(names(hsa_C2_cgp[["gene_sets"]] %in% names(hsa_C2_cgp[["descriptions"]]))))
})

test_that("`get_mgsigdb_gsets()` errors work", {
  all_collections <- c("H", "C1", "C2", "C3", "C4", "C5", "C6", "C7")
  expect_error(pathfindR:::get_mgsigdb_gsets(collection = "INVALID"),
               paste0("`collection` should be one of ",
                      paste(dQuote(all_collections), collapse = ", ")))
  skip_on_cran()
  species <- "Homo sapiens"
  collection <- "C2"
  subcollection <- "INVALID"
  expect_error(pathfindR:::get_mgsigdb_gsets(species = species,
                                             collection = collection,
                                             subcollection = subcollection),
               "unknown subcategory")
})

# get_gene_sets_list ------------------------------------------------------
test_that("`get_gene_sets_list()` works", {
  expect_error(gsets <- get_gene_sets_list("Wiki"),
               "As of this version, this function is implemented to get data from KEGG, Reactome and MSigDB only")

  skip_on_cran()
  expect_silent(kegg <- get_gene_sets_list(org_code = "vcn"))
  expect_message(rctm <- get_gene_sets_list("Reactome"))
  expect_silent(msig <- get_gene_sets_list("MSigDB",
                                           species = "Mus musculus",
                                           collection = "C3",
                                           subcollection = "MIR:MIR_Legacy"))
})
options(download.file.method = org_met, download.file.extra=org_extra)

