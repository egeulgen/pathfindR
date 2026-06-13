test_that("`return_pin_path()` -- returns the absolute path to PIN file", {
  mockery::stub(return_pin_path, "utils::getFromNamespace", list())
  mockery::stub(return_pin_path, "lapply", list(data.frame(
    V1 = paste0("G", 1:10),
    V2 = "pp", V3 = paste0("G", 2:11)
  ), data.frame(
    V1 = paste0("G", 3:5), V2 = "pp",
    V3 = paste0("G", 5:7)
  )))
  expect_silent(path2file <- return_pin_path("Biogrid"))
  expect_true(file.exists(path2file))

  custom_pin <- read.delim(path2file, header = FALSE)
  custom_pin$V1 <- tolower(custom_pin$V1)
  custom_sif_path <- file.path(tempdir(check = TRUE), "tmp_PIN.sif")
  utils::write.table(custom_pin, custom_sif_path,
    sep = "\t", row.names = FALSE,
    col.names = FALSE, quote = FALSE
  )
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
  write.table(invalid_custom_sif, invalid_sif_path,
    sep = "\t", col.names = FALSE,
    row.names = FALSE
  )
  expect_error(return_pin_path(invalid_sif_path), "The second column of the PIN file must all be \"pp\" ")

  # invalid option
  valid_opts <- c(
    "Biogrid", "STRING", "GeneMania", "IntAct", "KEGG", "mmu_STRING",
    "/path/to/custom/SIF"
  )
  expect_error(return_pin_path("INVALID"), paste0(
    "The chosen PIN must be one of:\n",
    paste(dQuote(valid_opts), collapse = ", ")
  ))
})
