#' Parse Active Subnetwork Search Output File and Filter the Subnetworks
#'
#' @param active_snws active subnetwork search results.
#' A list containing input \code{subnetworks} (nodes) and \code{scores} (score).
#' @param sig_genes_vec vector of significant gene symbols. In the scope of this
#'   package, these are the input genes that were used for active subnetwork search
#' @param score_quan_thr active subnetwork score quantile threshold. Must be
#' between 0 and 1 or set to -1 for not filtering. (Default = 0.8)
#' @param sig_gene_thr threshold for the minimum proportion of significant genes in
#' the subnetwork (Default = 0.02) If the number of genes to use as threshold is
#' calculated to be < 2 (e.g. 50 signif. genes x 0.01 = 0.5), the threshold number
#' is set to 2
#'
#' @return A list containing \code{subnetworks}: a list of of genes in every
#' active subnetwork that has a score greater than the \code{score_quan_thr}th
#' quantile and that contains at least \code{sig_gene_thr} of significant genes
#' and \code{scores} the score of each filtered active subnetwork
#' @export
#'
#' @seealso See \code{\link{run_pathfindR}} for the wrapper function of the
#'   pathfindR enrichment workflow
#'
#' @examples
#' filtered <- filter_active_subnetworks(
#'   active_snws = example_unfiltered_snws,
#'   sig_genes_vec = example_pathfindR_input$Gene.symbol
#' )
filter_active_subnetworks <- function(active_snws, sig_genes_vec, score_quan_thr = 0.8,
                                      sig_gene_thr = 0.02) {
  ## Arg. checks
  if (!is.atomic(sig_genes_vec)) {
    stop("`sig_genes_vec` should be a vector")
  }

  if (!is.numeric(score_quan_thr)) {
    stop("`score_quan_thr` should be numeric")
  }
  if (score_quan_thr != -1 & (score_quan_thr > 1 | score_quan_thr < 0)) {
    stop("`score_quan_thr` should be in [0, 1] or -1 (if not filtering)")
  }

  if (!is.numeric(sig_gene_thr)) {
    stop("`sig_gene_thr` should be numeric")
  }
  if (sig_gene_thr < 0 | sig_gene_thr > 1) {
    stop("`sig_gene_thr` should be in [0, 1]")
  }

  if (length(active_snws) == 0) {
    return(NULL)
  }

  score_vec <- c()
  subnetworks <- list()
  for (i in seq_along(active_snws)) {
    snw <- active_snws[[i]]

    score_vec <- c(score_vec, snw$score)
    subnetworks[[i]] <- snw$nodes
  }

  # keep subnetworks with score over the 'score_quan_thr'th quantile
  if (score_quan_thr == -1) {
    score_thr <- min(score_vec) - 1
  } else {
    score_thr <- stats::quantile(score_vec, score_quan_thr)
  }
  cond <- as.numeric(score_vec) > as.numeric(score_thr)
  subnetworks <- subnetworks[cond]
  score_vec <- as.numeric(score_vec)[cond]

  # select subnetworks containing at least 'sig_gene_thr' of significant
  # genes
  snw_sig_counts <- vapply(subnetworks, function(snw_genes) {
    sum(base::toupper(snw_genes) %in% base::toupper(sig_genes_vec))
  }, 1)
  sig_gene_num_thr <- sig_gene_thr * length(sig_genes_vec)
  sig_gene_num_thr <- max(2, sig_gene_num_thr)
  cond <- (snw_sig_counts >= sig_gene_num_thr)
  subnetworks <- subnetworks[cond]
  score_vec <- score_vec[cond]

  return(list(subnetworks = subnetworks, scores = score_vec))
}
