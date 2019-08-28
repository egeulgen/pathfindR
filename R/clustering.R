#' Create Kappa Statistics Matrix
#'
#' @param enrichment_res data frame of pathway enrichment results. Must-have
#' columns are "Pathway" or "ID", "Down_regulated", and "Up_regulated". If
#' `use_active_snw_genes = TRUE`, "non_DEG_Active_Snw_Genes" must also be
#' provided.
#' @param use_names boolean to indicate whether to use pathway names instead of
#' IDs (default = FALSE, i.e. use IDs)
#' @param use_active_snw_genes boolean to indicate whether or not to use
#' non-input active subnetwork genes
#' in the calculation of kappa statistics (default = FALSE,
#' i.e. use only affected genes)
#'
#' @return a matrix of kappa statistics between each term in the
#' enrichment results.
#'
#' @export
#'
#' @examples
#' sub_df <- RA_output[1:3,]
#' create_kappa_matrix(sub_df)
create_kappa_matrix <- function(enrichment_res,
                                use_names = FALSE,
                                use_active_snw_genes = FALSE) {

  ### Initial steps
  # Column to use for gene set names
  pw_id <- ifelse(use_names,
                  which(colnames(enrichment_res) == "Pathway"),
                  which(colnames(enrichment_res) == "ID"))

  # list of genes
  down_idx <- which(colnames(enrichment_res) == "Down_regulated")
  up_idx <- which(colnames(enrichment_res) == "Up_regulated")

  genes_lists <- apply(enrichment_res, 1, function(x)
    c(unlist(strsplit(as.character(x[up_idx]), ", ")),
      unlist(strsplit(as.character(x[down_idx]), ", "))))

  if (use_active_snw_genes) {

    if (!"non_DEG_Active_Snw_Genes" %in% colnames(enrichment_res))
      stop("No column named `non_DEG_Active_Snw_Genes`,
           please execute `run_pathfindR` with
           `list_active_snw_genes = TRUE`!")

    active_idx <- which(colnames(enrichment_res) == "non_DEG_Active_Snw_Genes")

    genes_lists <- mapply(function(x, y)
      c(x, unlist(strsplit(as.character(y), ", "))),
                          genes_lists, enrichment_res[, active_idx])
  }

  # Exclude zero-length gene sets
  excluded_idx <- which(vapply(genes_lists, length, 1) == 0)
  if (length(excluded_idx) != 0) {
    genes_lists <- genes_lists[-excluded_idx]
    enrichment_res <- enrichment_res[-excluded_idx, ]
  }

  ### Create Kappa Matrix
  all_genes <- unique(unlist(genes_lists, use.names = FALSE))
  N <- nrow(enrichment_res)
  pw_names <- enrichment_res[, pw_id]

  kappa_mat <- matrix(0, nrow = N, ncol = N,
                      dimnames = list(pw_names, pw_names))
  diag(kappa_mat) <- 1

  for (i in 1:(N - 1)) {
    for (j in (i+1):N) {
      genes_i <- genes_lists[[i]]
      genes_j <- genes_lists[[j]]

      C1_1 <- length(intersect(genes_i, genes_j))
      C0_0 <- sum(!all_genes %in% genes_i & !all_genes %in% genes_j)
      C0_1 <- sum(all_genes %in% genes_i & !all_genes %in% genes_j)
      C1_0 <- sum(!all_genes %in% genes_i & all_genes %in% genes_j)

      tot <- sum(C0_0, C0_1, C1_0, C1_1)

      observed <- (C1_1 + C0_0) / tot
      chance <- (C1_1 + C0_1) * (C1_1 + C1_0) + (C1_0 + C0_0) * (C0_1 + C0_0)
      chance <- chance / tot ^2
      kappa_mat[j, i] <- kappa_mat[i, j] <- (observed - chance) / (1 - chance)
    }
  }

  return(kappa_mat)
}


#' Hierarchical Clustering of Pathways
#'
#' @param kappa_mat matrix of kappa statistics (output of `create_kappa_matrix`)
#' @param enrichment_res data frame of pathway enrichment results
#' @param use_names boolean to indicate whether to use pathway names instead of
#' IDs (default = FALSE, i.e. use IDs)
#' @param clu_method the agglomeration method to be used
#' (default = "average", see `?hclust`)
#' @param plot_hmap boolean to indicate whether to plot the kappa statistics
#' heatmap or not (default = FALSE)
#' @param plot_dend boolean to indicate whether to plot the clustering
#' dendrogram partitioned into the optimal number of clusters (default = TRUE)
#'
#' @details The function initially performs hierarchical clustering
#' of the terms in `enrichment_res` using the kappa statistics
#' (defining the distance as `1 - kappa_statistic`). Next,
#' the clustering dendrogram is cut into k = 2, 3, ..., n - 1 clusters
#' (where n is the number of terms). The optimal number of clusters is
#' determined as the k value which yields the highest average silhouette width.
#'
#' @return a vector of clusters for each term in the enrichment results.
#' @export
#'
#' @examples
#' \dontrun{
#'hierarchical_pw_clustering(kappa_mat, enrichment_res)
#'hierarchical_pw_clustering(kappa_mat, enrichment_res, method = "complete")
#' }
hierarchical_pw_clustering <- function(kappa_mat, enrichment_res,
                                       use_names = FALSE,
                                       clu_method = "average",
                                       plot_hmap = FALSE, plot_dend = TRUE) {
  ### Set ID/Name index
  pw_id <- ifelse(use_names,
                  which(colnames(enrichment_res) == "Pathway"),
                  which(colnames(enrichment_res) == "ID"))

  ### Add excluded (zero-length) genes
  kappa_mat2 <- kappa_mat
  cond <- !enrichment_res[, pw_id] %in% rownames(kappa_mat2)
  outliers <- enrichment_res[cond, pw_id]
  outliers_mat <- matrix(-1, nrow = nrow(kappa_mat2), ncol = length(outliers),
                         dimnames = list(rownames(kappa_mat2), outliers))
  kappa_mat2 <- cbind(kappa_mat2, outliers_mat)
  outliers_mat <- matrix(-1, nrow = length(outliers), ncol = ncol(kappa_mat2),
                         dimnames = list(outliers, colnames(kappa_mat2)))
  kappa_mat2 <- rbind(kappa_mat2, outliers_mat)

  ### Hierarchical clustering
  clu <- stats::hclust(stats::as.dist(1 - kappa_mat2), method = clu_method)

  if (plot_hmap) {
    stats::heatmap(kappa_mat2,
                   distfun = function(x) stats::as.dist(1 - x),
                   hclustfun = function(x) stats::hclust(x, method = clu_method)
                   )
  }

  ### Choose optimal k
  kmax <- nrow(kappa_mat2) - 1
  avg_sils <- c()
  for (i in 2:kmax)
    avg_sils <- c(avg_sils, fpc::cluster.stats(stats::as.dist(1 - kappa_mat2),
                                               stats::cutree(clu, k = i),
                                               silhouette = TRUE)$avg.silwidth)

  k_opt <- (2:kmax)[which.max(avg_sils)]

  message(paste("The maximum average silhouette width was",
                round(max(avg_sils), 2),
                "for k =", k_opt, "\n\n"))

  if (plot_dend) {
    plot(clu)
    stats::rect.hclust(clu, k = k_opt)
  }

  ### Return clusters
  clusters <- stats::cutree(clu, k = k_opt)

  return(clusters)
}

#' Heuristic Fuzzy Multiple-linkage Partitioning of Pathways
#'
#' @param kappa_mat matrix of kappa statistics (output of `create_kappa_matrix`)
#' @param enrichment_res data frame of pathway enrichment results
#' @param kappa_threshold threshold for kappa statistics, defining strong
#' relation (default = 0.35)
#' @param use_names boolean to indicate whether to use pathway names instead of
#' IDs (default = FALSE, i.e. use IDs)
#'
#' @details The fuzzy clustering algorithm was implemented based on:
#' Huang DW, Sherman BT, Tan Q, et al. The DAVID Gene Functional
#' Classification Tool: a novel biological module-centric algorithm to
#' functionally analyze large gene lists. Genome Biol. 2007;8(9):R183.
#'
#' @return a boolean matrix of cluster assignments. Each row corresponds to a
#' term, each column corresponds to a cluster.
#' @export
#'
#' @examples
#' \dontrun{
#' fuzzy_pw_clustering(kappa_mat, enrichment_res)
#' fuzzy_pw_clustering(kappa_mat, enrichment_res, kappa_threshold = 0.45)
#' }
fuzzy_pw_clustering <- function(kappa_mat, enrichment_res,
                                kappa_threshold = 0.35, use_names = FALSE) {

  ### Check that the kappa threshold is numeric
  if (!is.numeric(kappa_threshold))
    stop("`kappa_threshold` must be numeric!")

  ### Set ID/Name index
  pw_id <- ifelse(use_names,
                  which(colnames(enrichment_res) == "Pathway"),
                  which(colnames(enrichment_res) == "ID"))

  ### Find Qualified Seeds
  qualified_seeds <- list()
  j <- 1
  for (i in base::seq_len(nrow(kappa_mat))) {
    current_term <- rownames(kappa_mat)[i]
    current_term_kappa <- kappa_mat[i, ]


    init_membership_cond <- current_term_kappa >= kappa_threshold
    if (sum(init_membership_cond) > 3) {
      related_terms <- names(current_term_kappa)[init_membership_cond]
      terms <- c(current_term, related_terms)
      related_kappa <- kappa_mat[rownames(kappa_mat) %in% terms,
                                 colnames(kappa_mat) %in% terms]
      diag(related_kappa) <- 0
      tight_relationship_cond <- sum(related_kappa >= kappa_threshold) /
        (nrow(related_kappa)^2) >= 0.5

      if (tight_relationship_cond) {
        qualified_seeds[[j]] <- related_terms
        names(qualified_seeds)[j] <- current_term
        j <- j + 1
      }
    }
  }

  ### Fuzzy Clustering
  clusters <- unique(qualified_seeds)
  i <- 1
  j <- i + 1
  while (i < length(clusters)) {
    common_terms <- intersect(clusters[[i]], clusters[[j]])
    all_terms <- union(clusters[[i]], clusters[[j]])

    if (length(common_terms) / length(all_terms) > 0.5 & i != j) {
      clusters[[i]] <- all_terms
      clusters[[j]] <- NULL
      i <- 1
      j <- i + 1
    } else if (j < length(clusters)) {
      j <- j + 1
    } else {
      i <- i + 1
      j <- 1
    }
  }

  ### Find Outliers
  cond <- !enrichment_res[, pw_id] %in% c(names(clusters), unlist(clusters))
  outliers <- enrichment_res[cond, pw_id]
  for (outlier in outliers) {
    clusters[[outlier]] <- outlier
  }
  ### Return Cluster Matrix
  names(clusters) <- base::seq_len(length(clusters))

  cluster_mat <- matrix(FALSE,
                        nrow = nrow(enrichment_res),
                        ncol = length(clusters),
                        dimnames = list(enrichment_res[, pw_id],
                                        names(clusters)))
  for(clu in names(clusters)) {
    clu_terms <- clusters[[clu]]
    cluster_mat[clu_terms, clu] <- TRUE
  }

  return(cluster_mat)
}


#' Graph Visualization of Pathway Clustering
#'
#' @param clu_obj clustering result (either a matrix obtained via
#' `fuzzy_pw_clustering` or a vector obtained via `hierarchical_pw_clustering`)
#' @param kappa_mat matrix of kappa statistics (output of `create_kappa_matrix`)
#' @param enrichment_res data frame of pathway enrichment results
#' @param kappa_threshold threshold for kappa statistics, defining strong
#' relation (default = 0.35)
#' @param use_names boolean to indicate whether to use pathway names instead of
#' IDs (default = FALSE, i.e. use IDs)
#'
#' @return Plots a graph diagram of clustering results. Each node is a term
#' from `enrichment_res`. Size of node corresponds to -log(lowest_p). Thickness
#' of the edges between nodes correspond to the kappa statistic between the two
#' terms. Color of each node corresponds to distinct clusters. For fuzzy
#' clustering, if a term is in multiple clusters, multiple colors are utilized.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' cluster_graph_vis(clu_obj, kappa_mat, enrichment_res)
#' }
cluster_graph_vis <- function(clu_obj, kappa_mat, enrichment_res,
                              kappa_threshold = 0.35, use_names = FALSE) {

  ### Init checks
  if (class(kappa_mat) != "matrix")
    stop("`kappa_mat` must be a matrix!")

  if (!isSymmetric.matrix(kappa_mat))
    stop("`kappa_mat` must be symmetric!")

  if (class(enrichment_res) != "data.frame")
    stop("`enrichment_res` must be a data.frame!")

  if (nrow(kappa_mat) != nrow(enrichment_res))
    stop("`kappa_mat` and `enrichment_res` must contain the same # of terms")


  ### Set ID/Name index
  pw_id <- ifelse(use_names,
                  which(colnames(enrichment_res) == "Pathway"),
                  which(colnames(enrichment_res) == "ID"))

  ## other checks
  if (!all(rownames(kappa_mat) %in% enrichment_res[, pw_id]))
    stop("Not all terms in `kappa_mat` and `enrichment_res` match!")

  ### For coloring nodes
  all_cols <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
                "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5",
                "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F",
                "#E5C494", "#B3B3B3", "#8DD3C7", "#FFFFB3", "#BEBADA",
                "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
                "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#A6CEE3",
                "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C",
                "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99",
                "#B15928")

  if (class(clu_obj) == "matrix") {
    ### Init checks
    if (!all(rownames(clu_obj) %in% colnames(kappa_mat)))
      stop("Not all terms in `clu_obj` present in `kappa_mat`!")

    ### Prep data
    # Remove weak links
    kappa_mat2 <- kappa_mat
    diag(kappa_mat2) <- 0
    kappa_mat2 <- ifelse(kappa_mat2 < kappa_threshold, 0, kappa_mat2)

    # Add missing pathways
    missing <- rownames(clu_obj)[!rownames(clu_obj) %in% colnames(kappa_mat2)]
    missing_mat <- matrix(0, nrow = nrow(kappa_mat2), ncol = length(missing),
                          dimnames = list(rownames(kappa_mat2), missing))
    kappa_mat2 <- cbind(kappa_mat2, missing_mat)
    missing <- rownames(clu_obj)[!rownames(clu_obj) %in% rownames(kappa_mat2)]
    missing_mat <- matrix(0, nrow = length(missing), ncol = ncol(kappa_mat2),
                          dimnames = list(missing, colnames(kappa_mat2)))
    kappa_mat2 <- rbind(kappa_mat2, missing_mat)

    ### Create Graph, Set Color, Size and Percentages
    values <- apply(clu_obj, 1, function(x) which(x))
    percs <- list()
    for (i in base::seq_len(length(values))) {
      percs[[i]] <- rep(1/length(values[[i]]), length(values[[i]]))
    }

    g <- igraph::graph_from_adjacency_matrix(kappa_mat2, weighted = TRUE)

    if (length(all_cols) < max(as.integer(colnames(clu_obj)))) {
      num_extra <- max(as.integer(colnames(clu_obj))) - length(all_cols)
      extra_colors <- grDevices::rainbow(num_extra)
      all_cols <- c(all_cols, extra_colors)
    }

    # Node shapes are either circle (single cluster) or pie (multiple clusters)
    igraph::V(g)$shape <- ifelse(vapply(percs, length, 1) > 1, "pie", "circle")

    # Node colors are cluster memberships
    cols <- lapply(values, function(x) all_cols[x])
    igraph::V(g)$color <- vapply(cols, function(x) x[1], "")

    # Node sizes are -log(lowest_p)
    p_idx <- match(names(igraph::V(g)), enrichment_res[, pw_id])
    transformed_p <- -log10(enrichment_res$lowest_p[p_idx])
    igraph::V(g)$size <- transformed_p * 2.5

    ### Plot Graph
    igraph::plot.igraph(g,
                        vertex.pie = percs,
                        vertex.pie.color = cols,
                        layout = igraph::layout_nicely(g),
                        edge.curved = FALSE,
                        vertex.label.dist = 0,
                        vertex.label.color = "black",
                        asp = 0,
                        vertex.label.cex = 0.7,
                        edge.width = igraph::E(g)$weight,
                        edge.arrow.mode = 0)

  } else if (class(clu_obj) == "integer") {
    ### Init checks
    if (!all(names(clu_obj) %in% colnames(kappa_mat)))
      stop("Not all terms in `clu_obj` present in `kappa_mat`!")

    ### Prep data
    # Remove weak links
    kappa_mat2 <- kappa_mat
    diag(kappa_mat2) <- 0
    kappa_mat2 <- ifelse(kappa_mat2 > kappa_threshold, kappa_mat2, 0)

    # Add missing pathways
    missing <- names(clu_obj)[!names(clu_obj) %in% colnames(kappa_mat2)]
    missing_mat <- matrix(0, nrow = nrow(kappa_mat2), ncol = length(missing),
                          dimnames = list(rownames(kappa_mat2), missing))
    kappa_mat2 <- cbind(kappa_mat2, missing_mat)
    missing <- names(clu_obj)[!names(clu_obj) %in% rownames(kappa_mat2)]
    missing_mat <- matrix(0, nrow = length(missing), ncol = ncol(kappa_mat2),
                          dimnames = list(missing, colnames(kappa_mat2)))
    kappa_mat2 <- rbind(kappa_mat2, missing_mat)

    ### Create Graph, Set Colors and Sizes
    g <- igraph::graph_from_adjacency_matrix(kappa_mat2, weighted = TRUE)

    igraph::V(g)$Clu <- clu_obj[match(igraph::V(g)$name, names(clu_obj))]

    if (length(all_cols) < max(as.integer(igraph::V(g)$Clu))) {
      num_extra <- max(clu_obj) - length(all_cols)
      extra_colors <- grDevices::rainbow(num_extra)
      all_cols <- c(all_cols, extra_colors)
    }

    # Node colors are cluster memberships
    igraph::V(g)$color <- all_cols[as.integer(igraph::V(g)$Clu)]

    # Node sizes are -log(lowest_p)
    p_idx <- match(names(igraph::V(g)), enrichment_res[, pw_id])
    transformed_p <- -log10(enrichment_res$lowest_p[p_idx])
    igraph::V(g)$size <- transformed_p * 2.5

    ### Plot graph
    igraph::plot.igraph(g,
                        layout = igraph::layout_nicely(g),
                        edge.curved = FALSE,
                        vertex.label.dist = 0,
                        vertex.label.color = "black",
                        asp = 0,
                        vertex.label.cex = 0.7,
                        edge.width = igraph::E(g)$weight,
                        edge.arrow.mode= 0)

  } else {
    stop("Invalid class for `clu_obj`!")
  }
}

#' Cluster Pathways
#'
#' @param enrichment_res data frame of pathway enrichment results
#' (result of `run_pathfindR`)
#' @param method Either "hierarchical" or "fuzzy". Details of clustering are
#' provided in the corresponding functions.
#' @param kappa_threshold threshold for kappa statistics, defining strong
#'relation (default = 0.35)
#' @param plot_clusters_graph boolean value indicate whether or not to plot
#' the graph diagram of clustering results (default = TRUE)
#' @param use_names boolean to indicate whether to use pathway names instead of
#' IDs (default = FALSE, i.e. use IDs)
#' @param use_active_snw_genes boolean to indicate whether or not to use
#' non-input active subnetwork genes in the calculation of kappa statistics
#' (default = FALSE, i.e. use only affected genes)
#' @param hclu_method the agglomeration method to be used
#' (default = "average", see `?hclust`)
#' @param plot_hmap boolean to indicate whether to plot the kappa statistics
#' heatmap or not (default = FALSE)
#' @param plot_dend boolean to indicate whether to plot the clustering
#' dendrogram partitioned into the optimal number of clusters (default = TRUE)
#'
#' @return a data frame of clustering results. For "hierarchical", the cluster
#' assignments (Cluster) and whether the term is representative of its cluster
#' (Status) is added as columns. For "fuzzy", terms that are in multiple
#' clusters are provided for each cluster. The cluster assignments (Cluster)
#' and whether the term is representative of its cluster (Status) is
#' added as columns.
#'
#' @export
#'
#' @examples
#' example_clustered <- cluster_pathways(RA_output[1:3,],
#' plot_clusters_graph = FALSE)
#' example_clustered <- cluster_pathways(RA_output[1:3,],
#' method = "fuzzy", plot_clusters_graph = FALSE)
#'
#' @seealso See \code{\link{hierarchical_pw_clustering}} for hierarchical
#' clustering of enriched terms.
#' See \code{\link{fuzzy_pw_clustering}} for fuzzy clustering of enriched terms.
cluster_pathways <- function(enrichment_res, method = "hierarchical",
                             kappa_threshold = 0.35,
                             plot_clusters_graph = TRUE,
                             use_names = FALSE, use_active_snw_genes = FALSE,
                             hclu_method = "average",
                             plot_hmap = FALSE, plot_dend = FALSE) {
  ### Argument Check
  if (!method %in% c("hierarchical", "fuzzy"))
    stop("the clustering `method` must either be \"hierarchical\" or \"fuzzy\"")

  ### Create Kappa Matrix
  kappa_mat <- create_kappa_matrix(enrichment_res = enrichment_res,
                                   use_names = use_names,
                                   use_active_snw_genes = use_active_snw_genes)

  ### Cluster Pathways
  if (method == "hierarchical") {
    clu_obj <- hierarchical_pw_clustering(kappa_mat = kappa_mat,
                                          enrichment_res = enrichment_res,
                                          use_names = use_names,
                                          clu_method = hclu_method,
                                          plot_hmap = plot_hmap,
                                          plot_dend = plot_dend)

  } else {
    clu_obj <- fuzzy_pw_clustering(kappa_mat = kappa_mat,
                                   enrichment_res = enrichment_res,
                                   kappa_threshold = kappa_threshold,
                                   use_names = use_names)
  }

  ### Graph Visualization of Clusters
  if (plot_clusters_graph)
    cluster_graph_vis(clu_obj = clu_obj, kappa_mat = kappa_mat,
                      enrichment_res = enrichment_res,
                      kappa_threshold = kappa_threshold, use_names = use_names)

  ### Returned Data Frame with Cluster Information
  clustered_df <- enrichment_res

  ### Set ID/Name index
  pw_id <- ifelse(use_names,
                  which(colnames(enrichment_res) == "Pathway"),
                  which(colnames(enrichment_res) == "ID"))

  if (method == "hierarchical") {
    ### Assign Clusters
    clu_idx <- match(clustered_df[, pw_id], names(clu_obj))
    clustered_df$Cluster <- clu_obj[clu_idx]
    clustered_df <- clustered_df[order(clustered_df$Cluster,
                                       clustered_df$lowest_p,
                                       decreasing = FALSE), ]

    tmp <- tapply(clustered_df$Pathway, clustered_df$Cluster, function(x) x[1])
    stat_cond <- clustered_df$Pathway %in% tmp
    clustered_df$Status <- ifelse(stat_cond, "Representative", "Member")
  } else {

    pws_list <- list()
    for (pw in rownames(clu_obj)) {
      pws_list[[pw]] <- which(clu_obj[pw, ])
    }
    ### Assign Clusters
    clustered_df2 <- c()
    for (i in base::seq_len(nrow(clustered_df))) {
      current_row <- clustered_df[i, ]
      current_clusters <- pws_list[[current_row[, pw_id]]]
      for (clu in current_clusters) {
        clustered_df2 <- rbind(clustered_df2,
                               data.frame(current_row, Cluster = clu))
      }
    }

    clustered_df <- clustered_df2
    clustered_df <- clustered_df[order(clustered_df$Cluster,
                                       clustered_df$lowest_p,
                                       decreasing = FALSE), ]

    tmp <- tapply(clustered_df$Pathway, clustered_df$Cluster, function(x) x[1])
    stat_cond <- clustered_df$Pathway %in% tmp
    clustered_df$Status <- ifelse(stat_cond, "Representative", "Member")
  }

  return(clustered_df)
}
