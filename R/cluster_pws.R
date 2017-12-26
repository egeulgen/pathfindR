#' Cluster Pathways
#'
#' @param pathway_genes List that contains genes for selected pathways that are
#'   to be clustered.
#' @param agg_method the agglomeration method to be used.
#'
#' @return Pairwise distance matrix. See "Chen, Y. A. et al. Integrated pathway
#'   clusters with coherent biological themes for target prioritisation. PLoS
#'   One 9, e99030, doi:10.1371/journal.pone.0099030 (2014)." for more details
#'   regarding the method of pathway clustering.
#' @seealso \code{\link[stats]{hclust}} for hierchical clustering,
#'   \code{\link[stats]{heatmap}} for drawing a heatmap.
#' @export
#'
#' @examples
cluster_pws <- function(pathway_genes, agg_method = "average") {
  ## Sort according to number of genes
  num_genes <- sapply(pathway_genes, length)
  pathway_genes <- pathway_genes[order(num_genes, decreasing = T)]

  # calculate overlap matrix ------------------------------------------------
  N <- length(pathway_genes)

  overlap_mat <- matrix(NA, nrow = N, ncol = N)
  diag(overlap_mat) <- 1
  for( i in 1:(N-1) )
  {
    pw1 <- pathway_genes[[i]]
    for( j in (i+1):N )
    {
      pw2 <- pathway_genes[[j]]
      tmp <- length(intersect(pw1, pw2))

      if(length(pw1) > length(pw2))
        tmp <- tmp/length(pw2)
      else
        tmp <- tmp/length(pw1)

      overlap_mat[i,j] <- overlap_mat[j,i] <- tmp
    }
  }

  # Calculate Pairwise Distances and Cluster --------------------------------------------
  cor_mat <- cor(overlap_mat)
  PWD_mat <- 1 - cor_mat

  ### Heatmap
  pw_names <- names(pathway_genes)
  rownames(PWD_mat) <- colnames(PWD_mat) <- pw_names
  pdf("pw_clu_heatmap.pdf", height = 15, width = 15)
  heatmap(PWD_mat, symm = T, distfun = function(x) as.dist(x), hclustfun = function(x) hclust(x, method = agg_method))
  dev.off()

  return(PWD_mat)
}

