#' Calculate Pairwise Distances Between Given Pathways
#'
#' @param pathway_ids Vector of IDs of pathways to be used for calculation of pairwise distances.
#' @param pathway_names Vector of names of pathways to be used for calculation of pairwise distances.
#' @param agg_method the agglomeration method to be used if plotting heatmap
#'   (see next argument, default: average).
#' @param plot_heatmap boolean value indicating whether or not to plot the heat
#'   map of pathway clustering (default: FALSE).
#' @param use_names boolean value indicating whether to use gene set names instead of gene set ids (default: FALSE)
#' @param custom_genes a list containing the genes involved in each custom pathway. Each element
#' is a vector of gene symbols located in the given pathway. Names correspond to
#' the ID of the pathway. Must be provided if pathways were obtained using custom
#' gene sets.
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
#' calculate_pwd(RA_output$ID)
calculate_pwd <- function(pathway_ids, pathway_names, agg_method = "average",
                          plot_heatmap = FALSE,
                          use_names = FALSE,
                          custom_genes = NULL) {
  ## Get genes for selected pathways
  if (grepl("^hsa", pathway_ids[1])) {
    pathway_genes <- pathfindR::kegg_genes[pathway_ids]
  } else if (grepl("^R-HSA-", pathway_ids[1])) {
    pathway_genes <- pathfindR::reactome_genes[pathway_ids]
  } else if (grepl("^BIOCARTA_", pathway_ids[1])) {
    pathway_genes <- pathfindR::biocarta_genes[pathway_ids]
  } else {
    if (all(pathway_ids %in% names(pathfindR::go_bp_genes))) {
      pathway_genes <- pathfindR::go_bp_genes[pathway_ids]
    } else if (all(pathway_ids %in% names(pathfindR::go_cc_genes))) {
      pathway_genes <- pathfindR::go_cc_genes[pathway_ids]
    } else if (all(pathway_ids %in% names(pathfindR::go_mf_genes))) {
      pathway_genes <- pathfindR::go_mf_genes[pathway_ids]
    } else {
      warning("Could not recognize IDs.\nAssuming the user is using custom gene sets")
      if (is.null(custom_genes))
        stop("`custom_genes` must be provided when using custom gene sets")
      pathway_genes <- custom_genes[pathway_ids]
    }
  }

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
  cor_mat[is.na(cor_mat)] <- 1
  PWD_mat <- 1 - cor_mat

  pw_ids <- names(pathway_genes)
  rownames(PWD_mat) <- colnames(PWD_mat) <- pw_ids

  ### Heatmap
  if (plot_heatmap) {
    message("Plotting the heatmap\n\n")
    to_plot <- PWD_mat
    if (use_names)
      rownames(to_plot) <- colnames(to_plot) <- pathway_names[match(colnames(to_plot), pathway_ids)]
    stats::heatmap(to_plot, symm = TRUE,
                   distfun = function(x) stats::as.dist(x),
                   hclustfun = function(x) stats::hclust(x, method = agg_method))
  }

  return(PWD_mat)
}
