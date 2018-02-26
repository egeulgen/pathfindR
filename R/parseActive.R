#' Parse Active Subnetwork Search Output File
#'
#' @param output_path path to the output of an Active Subnetwork Search.
#' @param signif_genes the vector of significant genes.
#'
#' @return A list of genes in every active subnetwork that has a score > 3 and
#'   that has at least 2 significant genes.
#' @export
#'
#' @seealso See \code{\link{run_pathfindR}} for the wrapper function of the
#'   pathfindR workflow
#'
#' @examples
#' \dontshow{
#' parseActiveSnwSearch(normalizePath(system.file("extdata/resultActiveSubnetworkSearch.txt",
#' package = "pathfindR")), RA_input$Gene.symbol)
#' }
#' \dontrun{
#' parseActiveSnwSearch("path/to/output", significant_genes)
#' }
parseActiveSnwSearch <- function(output_path, signif_genes) {

  output <- readLines(output_path)

  score <- c()
  subnetworks <- list()
  for (i in 1:length(output)) {
    snw <- output[[i]]

    snw <- unlist(strsplit(snw, " "))

    score <- c(score, snw[1])
    subnetworks[[i]] <- snw[-1]
  }

  # keep subnetworks with score > 3
  subnetworks <- subnetworks[score > 3]
  # select subnetworks with at least 2 significant genes
  cond <- sapply(subnetworks, function(x) sum(x %in% signif_genes)) >= 2
  subnetworks <- subnetworks[cond]

  return(subnetworks)
}
