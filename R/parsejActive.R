#' Parse jActive Output File
#'
#' @param jactive_output the output of a jActive run.
#' @param signif_genes the vector of significant genes.
#'
#' @return A list of genes in every active subnetwork that has a score > 0 and
#'   that has at least 2 significant genes.
#' @export
#'
#' @seealso See \code{\link{run_pathfindr}} for the wrapper function of the
#'   pathfindr workflow
#'
#' @examples
#' \dontrun{
#' parsejActive(output, significant_genes)
#' }
parsejActive <- function(jactive_output, signif_genes) {
  ends <- grep("SPOTPvaluesig", jactive_output$V1)
  starts <- c(1, ends + 1)
  starts <- starts[-length(starts)]
  ends <- ends - 1

  subnetworks <- list()
  scores <- c()
  for (i in 1:length(starts)){
    subnetworks[[i]] <- jactive_output$V1[starts[i]:ends[i]][-1]
    scores <- c(scores, as.numeric(jactive_output$V1[starts[i]:ends[i]][1]))
  }

  # keep snws with score > 3
  subnetworks <- subnetworks[scores > 3]
  # select subnetworks with at least 2 significant genes
  # cond <- sapply(subnetworks, function(x) sum(x %in% signif_genes)) >= 2
  # subnetworks <- subnetworks[cond]

  return(subnetworks)
}
