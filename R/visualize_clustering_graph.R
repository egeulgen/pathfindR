#' Graph Visualization of Clustered Enriched Terms
#'
#' @param clu_obj clustering result (either a matrix obtained via
#' \code{\link{hierarchical_term_clustering}} or \code{\link{fuzzy_term_clustering}}
#' `fuzzy_term_clustering` or a vector obtained via `hierarchical_term_clustering`)
#' @inheritParams fuzzy_term_clustering
#' @param vertex.label.cex font size for vertex labels; it is interpreted as a multiplication factor of some device-dependent base font size (default = 0.7)
#' @param vertex.size.scaling scaling factor for the node size (default = 2.5)
#'
#' @return Plots a graph diagram of clustering results. Each node is an enriched term
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
cluster_graph_vis <- function(clu_obj, kappa_mat, enrichment_res, kappa_threshold = 0.35,
                              use_description = FALSE, vertex.label.cex = 0.7, vertex.size.scaling = 2.5) {
  ### Set ID/Name index
  chosen_id <- ifelse(use_description, which(colnames(enrichment_res) == "Term_Description"),
    which(colnames(enrichment_res) == "ID")
  )

  ### For coloring nodes
  all_cols <- c(
    "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33",
    "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3",
    "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", "#8DD3C7", "#FFFFB3", "#BEBADA",
    "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD",
    "#CCEBC5", "#FFED6F", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
    "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928"
  )

  if (is.matrix(clu_obj)) {
    ### Argument checks
    if (!all(rownames(clu_obj) %in% colnames(kappa_mat))) {
      stop("Not all terms in `clu_obj` present in `kappa_mat`!")
    }

    ### Prep data Remove weak links
    kappa_mat2 <- kappa_mat
    diag(kappa_mat2) <- 0
    kappa_mat2 <- ifelse(kappa_mat2 < kappa_threshold, 0, kappa_mat2)

    # Add missing terms
    missing <- rownames(clu_obj)[!rownames(clu_obj) %in% colnames(kappa_mat2)]
    missing_mat <- matrix(0,
      nrow = nrow(kappa_mat2), ncol = length(missing),
      dimnames = list(rownames(kappa_mat2), missing)
    )
    kappa_mat2 <- cbind(kappa_mat2, missing_mat)
    missing <- rownames(clu_obj)[!rownames(clu_obj) %in% rownames(kappa_mat2)]
    missing_mat <- matrix(0,
      nrow = length(missing), ncol = ncol(kappa_mat2),
      dimnames = list(missing, colnames(kappa_mat2))
    )
    kappa_mat2 <- rbind(kappa_mat2, missing_mat)

    ### Create Graph, Set Color, Size and Percentages
    values <- apply(clu_obj, 1, function(x) which(x))
    percs <- list()
    for (i in base::seq_len(length(values))) {
      percs[[i]] <- rep(1 / length(values[[i]]), length(values[[i]]))
    }

    g <- igraph::graph_from_adjacency_matrix(kappa_mat2, weighted = TRUE)

    if (length(all_cols) < max(as.integer(colnames(clu_obj)))) {
      num_extra <- max(as.integer(colnames(clu_obj))) - length(all_cols)
      extra_colors <- grDevices::rainbow(num_extra)
      all_cols <- c(all_cols, extra_colors)
    }

    # Node shapes are either circle (single cluster) or pie (multiple
    # clusters)
    igraph::V(g)$shape <- ifelse(vapply(percs, length, 1) > 1, "pie", "circle")

    # Node colors are cluster memberships
    cols <- lapply(values, function(x) all_cols[x])
    igraph::V(g)$color <- vapply(cols, function(x) x[1], "")

    # Node sizes are -log(lowest_p)
    p_idx <- match(names(igraph::V(g)), enrichment_res[, chosen_id])
    transformed_p <- -log10(enrichment_res$lowest_p[p_idx])
    igraph::V(g)$size <- transformed_p * vertex.size.scaling

    ### Plot Graph
    igraph::plot.igraph(g,
      vertex.pie = percs, vertex.pie.color = cols, layout = igraph::layout_nicely(g),
      edge.curved = FALSE, vertex.label.dist = 0, vertex.label.color = "black",
      asp = 1, vertex.label.cex = vertex.label.cex, edge.width = igraph::E(g)$weight,
      edge.arrow.mode = 0
    )
  } else if (is.integer(clu_obj)) {
    ### Argument checks
    if (!all(names(clu_obj) %in% colnames(kappa_mat))) {
      stop("Not all terms in `clu_obj` present in `kappa_mat`!")
    }

    ### Prep data Remove weak links
    kappa_mat2 <- kappa_mat
    diag(kappa_mat2) <- 0
    kappa_mat2 <- ifelse(kappa_mat2 > kappa_threshold, kappa_mat2, 0)

    # Add missing terms
    missing <- names(clu_obj)[!names(clu_obj) %in% colnames(kappa_mat2)]
    missing_mat <- matrix(0,
      nrow = nrow(kappa_mat2), ncol = length(missing),
      dimnames = list(rownames(kappa_mat2), missing)
    )
    kappa_mat2 <- cbind(kappa_mat2, missing_mat)
    missing <- names(clu_obj)[!names(clu_obj) %in% rownames(kappa_mat2)]
    missing_mat <- matrix(0,
      nrow = length(missing), ncol = ncol(kappa_mat2),
      dimnames = list(missing, colnames(kappa_mat2))
    )
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
    p_idx <- match(names(igraph::V(g)), enrichment_res[, chosen_id])
    transformed_p <- -log10(enrichment_res$lowest_p[p_idx])
    igraph::V(g)$size <- transformed_p * vertex.size.scaling

    ### Plot graph
    igraph::plot.igraph(g,
      layout = igraph::layout_nicely(g), edge.curved = FALSE,
      vertex.label.dist = 0, vertex.label.color = "black", asp = 0, vertex.label.cex = vertex.label.cex,
      edge.width = igraph::E(g)$weight, edge.arrow.mode = 0
    )
  } else {
    stop("Invalid class for `clu_obj`!")
  }
}
