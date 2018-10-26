#' Parse Active Subnetwork Search Output File
#'
#' @param active_snw_path path to the output of an Active Subnetwork Search.
#' @param signif_genes the vector of significant genes.
#' @param score_quan_thr active subnetwork score quantile threshold (Default = 0.80)
#' @param sig_gene_thr threshold for minimum number of significant genes (Default = 10)
#'
#' @return A list of genes in every active subnetwork that has a score greater than
#' the `score_quan_thr`th quantile and that has at least `sig_gene_thr` significant genes.
#' @export
#'
#' @seealso See \code{\link{run_pathfindR}} for the wrapper function of the
#'   pathfindR workflow
#'
#' @examples
#' \dontshow{
#' parseActiveSnwSearch(normalizePath(system.file("extdata/resultActiveSubnetworkSearch.txt",
#' package = "pathfindR")), pathfindR::RA_input$Gene.symbol)
#' }
#' \dontrun{
#' parseActiveSnwSearch("path/to/output", significant_genes)
#' }
parseActiveSnwSearch <- function(active_snw_path, signif_genes,
                                 score_quan_thr = 0.80, sig_gene_thr = 10) {

  output <- readLines(active_snw_path)

  if (length(output) == 0)
    return(NULL)

  score_vec <- c()
  subnetworks <- list()
  for (i in 1:length(output)) {
    snw <- output[[i]]

    snw <- unlist(strsplit(snw, " "))

    score_vec <- c(score_vec, as.numeric(snw[1]))
    subnetworks[[i]] <- snw[-1]
  }

  # keep subnetworks with score over the "score_quan_thr"th quantile
  score_thr <- stats::quantile(score_vec, score_quan_thr)
  subnetworks <- subnetworks[as.numeric(score_vec) > as.numeric(score_thr)]
  # select subnetworks with at least sig_gene_thr significant genes
  snw_sig_counts <- sapply(subnetworks, function(snw) sum(snw %in% signif_genes))
  cond <- (snw_sig_counts >= sig_gene_thr)
  subnetworks <- subnetworks[cond]

  return(subnetworks)
}
