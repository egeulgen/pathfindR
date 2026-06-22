#' Return The Path to Given Protein-Protein Interaction Network (PIN)
#'
#' This function returns the absolute path/to/PIN.sif. While the default PINs are
#' 'Biogrid', 'STRING', 'GeneMania', 'IntAct', 'KEGG' and 'mmu_STRING'. The user can also
#' use any other PIN by specifying the 'path/to/PIN.sif'. All PINs to be used
#' in this package must formatted as SIF files: i.e. have 3 columns with no
#' header, no row names and be tab-separated. Columns 1 and 3 must be
#' interactors' gene symbols, column 2 must be a column with all
#' rows consisting of 'pp'.
#'
#' @param pin_name_path Name of the chosen PIN or absolute/path/to/PIN.sif. If PIN name,
#'   must be one of c('Biogrid', 'STRING', 'GeneMania', 'IntAct', 'KEGG', 'mmu_STRING'). If
#'   path/to/PIN.sif, the file must comply with the PIN specifications. (Default = 'Biogrid')
#'
#' @return The absolute path to chosen PIN.
#'
#' @export
#' @seealso See \code{\link{run_pathfindR}} for the wrapper function of the
#'   pathfindR workflow
#' @examples
#' \dontrun{
#' pin_path <- return_pin_path("GeneMania")
#' }
return_pin_path <- function(pin_name_path = "Biogrid") {
  ## Default PINs
  valid_opts <- c(
    "Biogrid", "STRING", "GeneMania", "IntAct", "KEGG", "mmu_STRING",
    "/path/to/custom/SIF"
  )
  if (pin_name_path %in% valid_opts[-length(valid_opts)]) {
    path <- file.path(tempdir(check = TRUE), paste0(pin_name_path, ".sif"))
    if (!file.exists(path)) {
      adj_list <- utils::getFromNamespace(paste0(tolower(pin_name_path), "_adj_list"),
        ns = "pathfindR.data"
      )

      pin_df <- lapply(seq_along(adj_list), function(i, nm, val) {
        data.frame(base::toupper(nm[[i]]), "pp", base::toupper(val[[i]]))
      }, val = adj_list, nm = names(adj_list))
      pin_df <- base::do.call("rbind", pin_df)
      utils::write.table(pin_df, path,
        sep = "\t", row.names = FALSE, col.names = FALSE,
        quote = FALSE
      )
    }
    path <- normalizePath(path)

    ## Custom PIN
  } else if (file.exists(suppressWarnings(normalizePath(pin_name_path)))) {
    path <- normalizePath(pin_name_path)
    pin <- utils::read.delim(file = path, quote = "", header = FALSE)
    if (ncol(pin) != 3) {
      stop("The PIN file must have 3 columns and be tab-separated")
    }

    if (any(pin[, 2] != "pp")) {
      stop("The second column of the PIN file must all be \"pp\" ")
    }

    if (any(grepl("[a-z]", pin[, 1])) | any(grepl("[a-z]", pin[, 3]))) {
      pin[, 1] <- base::toupper(pin[, 1])
      pin[, 3] <- base::toupper(pin[, 3])

      path <- file.path(tempdir(check = TRUE), "custom_PIN.sif")
      utils::write.table(pin, path,
        sep = "\t", row.names = FALSE, col.names = FALSE,
        quote = FALSE
      )
      path <- normalizePath(path)
    }
  } else {
    stop("The chosen PIN must be one of:\n", paste(dQuote(valid_opts), collapse = ", "))
  }
  return(path)
}
