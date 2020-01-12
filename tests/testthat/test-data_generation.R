##################################################
## Project: pathfindR
## Script purpose: Testthat testing script for
## data generation functions
## Date: Jan 12, 2020
## Author: Ege Ulgen
##################################################

# process_pin -------------------------------------------------------------
test_that("`process_pin()` works correctly", {
  # test removal of self-interactions
  toy_pin <- data.frame(Interactor_A = LETTERS[1:6],
                        Interactor_B = LETTERS[c(7:10, 5, 6)],
                        stringsAsFactors = FALSE)
  expect_true(nrow(pathfindR:::process_pin(toy_pin)) == nrow(toy_pin) - 2)

  # test removal of duplicated interactions
  toy_pin <- data.frame(Interactor_A = c(rep(LETTERS[1:3], 2), "X", "Y"),
                        Interactor_B = c(rep(LETTERS[4:6], 2), "Y", "X"),
                        stringsAsFactors = FALSE)
  expect_true(nrow(pathfindR:::process_pin(toy_pin)) == nrow(toy_pin) / 2)
})

# get_biogrid_pin ---------------------------------------------------------
test_that("`get_biogrid_pin()` returns a path to a valid PIN file", {
  skip_on_cran()
  # Latest release
  pin_path <- get_biogrid_pin()
  pin_df <- read.delim(pin_path, header = FALSE, stringsAsFactors = FALSE)
  expect_true(ncol(pin_df) == 3)
  expect_true(all(pin_df[, 2] == "pp"))

  # Specific release
  pin_path <- get_biogrid_pin(release = "3.5.180")
  pin_df <- read.delim(pin_path, header = FALSE, stringsAsFactors = FALSE)
  expect_true(ncol(pin_df) == 3)
  expect_true(all(pin_df[, 2] == "pp"))
})

test_that("`get_biogrid_pin()` error check works", {
  # invalid organism error
  expect_error(get_biogrid_pin(org = "Hsapiens"),
               paste("Hsapiens is not a valid Biogrid organism.",
                     "Available organisms are listed on: https://wiki.thebiogrid.org/doku.php/statistics"))
})

# get_pin_file ------------------------------------------------------------
test_that("`get_pin_file()` works", {
  # unimplemented error
  expect_error(get_pin_file(source = "STRING"),
"As of this version, this function is implemented to get data from BioGRID only")

  skip_on_cran()
  pin_path <- get_pin_file()
  pin_df <- read.delim(pin_path, header = FALSE, stringsAsFactors = FALSE)
  expect_true(ncol(pin_df) == 3)
  expect_true(all(pin_df[, 2] == "pp"))
})

# get_kegg_gsets ----------------------------------------------------------
test_that("`get_kegg_gsets() works`", {
  skip_on_cran()
  # hsa - default
  expect_silent(hsa_kegg <- get_kegg_gsets())
  expect_length(hsa_kegg, 2)
  expect_true(all(names(hsa_kegg) == c("gene_sets", "descriptions")))
  expect_true(length(hsa_kegg[["gene_sets"]]) == length(hsa_kegg[["descriptions"]]))

  # ptr
  expect_silent(ptr_kegg <- get_kegg_gsets("ptr"))
  expect_length(ptr_kegg, 2)
  expect_true(all(names(ptr_kegg) == c("gene_sets", "descriptions")))
  expect_true(length(ptr_kegg[["gene_sets"]]) == length(ptr_kegg[["descriptions"]]))
})

# get_reactome_gsets ------------------------------------------------------
test_that("`get_reactome_gsets()` works", {
  expect_silent(reactome <- get_reactome_gsets())
  expect_length(reactome, 2)
  expect_true(all(names(reactome) == c("gene_sets", "descriptions")))
  expect_true(length(reactome[["gene_sets"]]) == length(reactome[["descriptions"]]))
})

# get_gene_sets_list ------------------------------------------------------
test_that("`get_gene_sets_list()` works", {
  expect_error(rctm <- get_gene_sets_list("Wiki"),
               "As of this version, this function is implemented to get data from KEGG and Reactome only")

  skip_on_cran()
  expect_silent(rctm <- get_gene_sets_list())
  expect_silent(rctm <- get_gene_sets_list("Reactome"))
})

