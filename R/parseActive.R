#' Parse Active Subnetwork Search Output File
#'
#' @param output_path path to the output of an Active Subnetwork Search.
#' @param signif_genes the vector of significant genes.
#' @param score_thr active subnetwork score threshold (Default = 3)
#' @param sig_gene_thr threshold for minimum number of significant genes (Default = 2)
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
parseActiveSnwSearch <- function(output_path, signif_genes,
                                 score_thr = 3, sig_gene_thr = 2) {

  output <- readLines(output_path)

  score <- c()
  subnetworks <- list()
  for (i in 1:length(output)) {
    snw <- output[[i]]

    snw <- unlist(strsplit(snw, " "))

    score <- c(score, snw[1])
    subnetworks[[i]] <- snw[-1]
  }

  # keep subnetworks with score > score_thr
  subnetworks <- subnetworks[score > score_thr]
  # select subnetworks with at least sig_gene_thr significant genes
  cond <- sapply(subnetworks, function(x) sum(x %in% signif_genes)) >= sig_gene_thr
  subnetworks <- subnetworks[cond]

  return(subnetworks)
}
