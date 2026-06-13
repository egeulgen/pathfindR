#' Configure Output Directory Name
#'
#' @inheritParams run_pathfindR
#'
#' @return /path/to/output/dir
configure_output_dir <- function(output_dir = NULL) {
  output_dir_init <- output_dir
  output_dir <- ifelse(is.null(output_dir), file.path(tempdir(check = TRUE), "pathfindR_results"),
    output_dir
  )
  dir_changed <- FALSE
  while (dir.exists(output_dir)) {
    output_dir <- sub("/$", "", output_dir)
    if (grepl("\\(\\d+\\)$", output_dir)) {
      output_dir <- unlist(strsplit(output_dir, "\\("))
      suffix <- as.numeric(sub("\\)", "", output_dir[2])) + 1
      output_dir <- paste0(output_dir[1], "(", suffix, ")")
    } else {
      output_dir <- paste0(output_dir, "(1)")
    }
    dir_changed <- TRUE
  }

  if (dir_changed & !is.null(output_dir_init)) {
    message(paste0(
      "There is already a directory named \"", output_dir_init,
      "\".\nWriting the result to \"", output_dir, "\" not to overwrite any previous results."
    ))
  }
  return(output_dir)
}
