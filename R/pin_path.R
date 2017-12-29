#' Return Path to Given Protein Interaction Network (PIN)
#'
#' @param pin_name Name of the chosen PIN. Must be one of c("Biogrid", "STRING",
#'   "GeneMania", "BioPlex"). Defaults to "GeneMania".
#'
#' @return A character value that contains the path to chosen PIN.
#'
#' @export
#'
#' @examples
#' pin_path <- return_pin_path("Biogrid")
return_pin_path <- function(pin_name = "BioPlex") {
  if (!pin_name %in% c("Biogrid", "STRING", "GeneMania", "BioPlex"))
    stop(paste0("The chosen PIN must be one of:\n",
                "Biogrid, STRING, GeneMania, BioPlex"))

  path <- normalizePath(system.file(paste0("data/", pin_name, ".sif"),
                            package = "pathfindr"))
  return(path)
}
