#' Calculate Pairwise Distances Between Given Pathways
#'
#' @param pathway_ids Vector of IDs of pathways selected to be clustered.
#' @param agg_method the agglomeration method to be used if plotting heatmap
#'   (see next argument, default: average).
#' @param plot_heatmap boolean value indicating whether or not to plot the heat
#'   map of pathway clustering (default: FALSE).
#'
#' @details See "Chen, Y. A. et al. Integrated pathway clusters with coherent
#'   biological themes for target prioritisation. PLoS One 9, e99030,
#'   doi:10.1371/journal.pone.0099030 (2014)." for details on the method of
#'   pathway clustering.
#'
#' @return Pairwise distance matrix. Optionally plots a heatmap of pathway
#'   clustering.
#' @seealso \code{\link[stats]{hclust}} for hierarchical clustering,
#'   \code{\link[stats]{heatmap}} for drawing a heatmap.
#' @export
#'
#' @examples
#' cluster_pathways(RA_output$ID)
cluster_pathways <- function(pathway_ids, agg_method = "average",
                             plot_heatmap = FALSE) {
  ## Get genes for selected pathways
  pathway_genes <- pathfindR::genes_by_pathway[pathway_ids]

  ## Sort according to number of genes
  num_genes <- sapply(pathway_genes, length)
  pathway_genes <- pathway_genes[order(num_genes, decreasing = TRUE)]

  # calculate overlap matrix ------------------------------------------------
  N <- length(pathway_genes)

  overlap_mat <- matrix(NA, nrow = N, ncol = N)
  diag(overlap_mat) <- 1
  for (i in 1:(N - 1)) {
    pw1 <- pathway_genes[[i]]
    for (j in (i + 1):N) {
      pw2 <- pathway_genes[[j]]
      tmp <- length(intersect(pw1, pw2))

      if (length(pw1) > length(pw2))
        tmp <- tmp / length(pw2)
      else
        tmp <- tmp / length(pw1)

      overlap_mat[i, j] <- overlap_mat[j, i] <- tmp
    }
  }

  # Calculate Pairwise Distances and Cluster -------------------------------
  cor_mat <- stats::cor(overlap_mat)
  PWD_mat <- 1 - cor_mat

  pw_names <- names(pathway_genes)
  rownames(PWD_mat) <- colnames(PWD_mat) <- pw_names

  ### Heatmap
  if (plot_heatmap)
    stats::heatmap(PWD_mat, symm = TRUE,
                   distfun = function(x) stats::as.dist(x),
            hclustfun = function(x) stats::hclust(x, method = agg_method))

  return(PWD_mat)
}
