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
    } else if (all(pathway_ids %in% names(pathfindR::go_all_genes)) {
      pathway_genes <- pathfindR::go_all_genes[pathway_ids]
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

#' Cluster Pathways and Partition the Dendrogram
#'
#' This function first calculates the pairwise distances between the
#' pathways in the \code{result_df} data frame. Next, using this distance
#' matrix, the pathways are clustered via hierarchical clustering. By default,
#' the average silhouette width for each possible number of clusters is
#' calculated. The optimal number of clusters is selected as the one with the
#' highest average silhouette width. The dendrogram is cut into this optimal
#' number of clusters, and the pathways with the lowest p value within each
#' cluster are chosen as representative pathways. If 'auto == FALSE", the user
#' can manually select at which height to cut the dendrogram via a shiny application.
#' See "Chen, Y. A. et al. Integrated pathway clusters with coherent biological
#' themes for target prioritisation. PLoS One 9, e99030,
#' doi:10.1371/journal.pone.0099030 (2014)." for details on the method of
#' pathway clustering.
#'
#' @param result_df data frame of enriched pathways. Must-have columns are: \enumerate{
#'   \item{ID}{KEGG ID of the enriched pathway}
#'   \item{lowest_p}{the lowest adjusted-p value of the given pathway over all iterations}
#'   \item{highest_p}{the highest adjusted-p value of the given pathway over all iterations}
#'   }
#' @param p_val_threshold p value threshold for filtering the pathways prior to clustering (default: 0.05)
#' @param auto boolean value indicating whether to select the optimal number of clusters
#' automatically. If FALSE, a shiny application is displayed, where the user can manually
#' partition the clustering dendrogram (default: TRUE).
#' @param agg_method the agglomeration method to be used if plotting heatmap. Must be one of "ward.D", "ward.D2",
#' "single", "complete", "average", "mcquitty", "median" or "centroid" (default: "average").
#' @param plot_heatmap boolean value indicating whether or not to plot the heat
#'   map of pathway clustering (default: FALSE).
#' @param plot_dend boolean value indicating whether or not to plot the dendrogram
#'   partitioned into the optimal number of clusters, shown by red rectangles (default: FALSE)
#' @param use_names boolean value indicating whether to use gene set names instead of gene set ids (default: FALSE)
#' @param custom_genes a list containing the genes involved in each custom pathway. Each element
#' is a vector of gene symbols located in the given pathway. Names correspond to
#' the ID of the pathway. Must be provided if `result_df` was generated using custom
#' gene sets.
#'
#' @return  If 'auto' is FALSE, manual partitioning can be performed. Via a shiny HTML document, the
#'   hierarchical clustering dendrogram is visualized. In this HTML document,
#'   the user can select the agglomeration method and the distance value at
#'   which to cut the tree. The resulting cluster assignments of the pathways
#'   along with annotation of representative pathways (chosen by smallest lowest
#'   p value) are presented as a table and this table can be saved as a csv
#'   file.
#'   If 'auto' is TRUE, automatic partitioning of clusters is performed. The function
#'   adds 2 additional columns to the input data frame and returns it: \describe{
#'   \item{Cluster}{the cluster to which the pathway is assigned}
#'   \item{Status}{whether the pathway is the "Representative" pathway in its cluster or only a "Member"}
#' }
#'
#' @import fpc
#' @import knitr
#' @import shiny
#' @import rmarkdown
#' @import stats
#' @export
#' @seealso See \code{\link{calculate_pwd}} for calculation of pairwise
#'   distances between enriched pathways. See \code{\link[stats]{hclust}}
#'   for more information on hierarchical clustering. See \code{\link{run_pathfindR}}
#'   for the wrapper function of the pathfindR enrichment workflow.
#'
#' @examples
#' ## Cluster pathways with p <= 0.01
#' choose_clusters(RA_output, p_val_threshold = 0.01)
choose_clusters <- function(result_df, p_val_threshold = 0.05, auto = TRUE, agg_method = "average",
                            plot_heatmap = FALSE, plot_dend = FALSE, use_names = FALSE, custom_genes = NULL) {
  ## argument checks
  if (!is.logical(auto))
    stop("The argument `auto` must be either TRUE or FALSE!")

  if (!is.logical(plot_heatmap))
    stop("The argument `plot_heatmap` must be either TRUE or FALSE!")

  valid <- c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")
  if (!agg_method %in% valid)
    stop("`agg_method` must be one of ward.D, ward.D2, single, complete, average, mcquitty, median or centroid!")

  if (!is.logical(plot_heatmap))
    stop("The argument `plot_dend` must be either TRUE or FALSE!")

  if (!is.numeric(p_val_threshold))
    stop("`p_val_threshold` must be a numeric value between 0 and 1")

  if (p_val_threshold > 1 | p_val_threshold < 0)
    stop("`p_val_threshold` must be between 0 and 1")


  ## Check if clustering should be performed
  if (nrow(result_df) < 3) {
    warning("There are less than 3 pathways in result_df so clustering is not performed!")
    result_df$Cluster <- 1:nrow(result_df)
    result_df$Status <- "Representative"
    return(result_df)
  }

  ## Filter for p value
  result_df <- result_df[result_df$highest_p <= p_val_threshold, ]

  ## Create PWD matrix
  message("Calculating pairwise distances between pathways\n\n")
  PWD_mat <- pathfindR::calculate_pwd(result_df$ID, result_df$Pathway,
                                      agg_method = agg_method,
                                      plot_heatmap = plot_heatmap,
                                      use_names,
                                      custom_genes)
  if (!auto) {
    message("Creating the shiny app\n\n")
    parameters <- list(df = result_df, mat = PWD_mat, use_names = use_names)
    rmarkdown::run(system.file("rmd/clustering.Rmd", package = "pathfindR"),
                   render_args = list(output_dir = ".", params = parameters))
  } else {
    ### Calculate PWDs and Cluster
    message("Clustering pathways\n\n")
    hclu <- stats::hclust(as.dist(PWD_mat), method = agg_method)

    ### Optimal k
    message("Calculating the optimal number of clusters (based on average silhouette width)\n\n")
    kmax <- nrow(PWD_mat) - 1
    avg_sils <- c()
    for (i in 2:kmax)
      avg_sils <- c(avg_sils, fpc::cluster.stats(stats::as.dist(PWD_mat),
                                                 stats::cutree(hclu, k = i),
                                                 silhouette = TRUE)$avg.silwidth)

    k_opt <- (2:kmax)[which.max(avg_sils)]
    if (plot_dend) {
      to_plot <- hclu
      if (use_names)
        to_plot$labels <- result_df$Pathway[match(to_plot$labels, result_df$ID)]
      graphics::plot(to_plot, hang = -1)
      stats::rect.hclust(to_plot, k = k_opt)
    }
    message(paste("The maximum average silhouette width was", round(max(avg_sils), 2),
                  "for k =", k_opt, "\n\n"))

    ### Return Optimal Clusters
    clusters <- cutree(hclu, k = k_opt)

    result_df$Cluster <- clusters[match(result_df$ID, names(clusters))]
    tmp <- result_df$lowest_p
    names(tmp) <- result_df$ID
    tmp <- tapply(tmp, result_df$Cluster, function(x) names(x)[which.min(x)])
    result_df$Status <- ifelse(result_df$ID %in% tmp, "Representative", "Member")

    result_df <- result_df[order(result_df$Status, decreasing = TRUE), ]
    result_df <- result_df[order(result_df$Cluster), ]

    message("Returning the resulting data frame\n\n")
    return(result_df)
  }
}
