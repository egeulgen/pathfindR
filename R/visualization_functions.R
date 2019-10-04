#' Create Diagrams for Enriched Terms
#'
#' @param result_df Data frame of enrichment results. Must-have columns for
#'  KEGG human pathway diagrams are: "ID" and "Term_Description".
#'  Must-have columns for the rest are: "Term_Description", "Up_regulated" and
#' "Down_regulated"
#' @param input_processed input data processed via \code{\link{input_processing}},
#'  not necessary for visualizations other than KEGG human pathway diagrams.
#' @param hsa_KEGG boolean to indicate whether human KEGG gene sets were used for
#'  enrichment analysis or not (default = \code{TRUE})
#' @param pin_name_path Name of the chosen PIN or path/to/PIN.sif. If PIN name,
#'  must be one of c("Biogrid", "GeneMania", "IntAct", "KEGG"). If
#'  path/to/PIN.sif, the file must comply with the PIN specifications. (default
#'  is "Biogrid")
#'
#' @return Depending on the argument \code{hsa_KEGG}, creates visualization of
#'  interactions of genes involved in the list of enriched terms in
#'  \code{result_df} and saves them in the folder "term_visualizations" under
#'  the current working directory.
#'
#'
#' @details For \code{hsa_KEGG = TRUE}, KEGG human pathway diagrams are created,
#' affected nodes colored by up/down regulation status.
#' For other gene sets, interactions of affected genes are determined (via a shortest-path
#' algorithm) and are visualized (colored by change status) using igraph.
#'
#'
#' @export
#'
#' @seealso See \code{\link{visualize_hsa_KEGG}} for the visualization function
#' of human KEGG diagrams. See \code{\link{visualize_term_interactions}} for the
#' visualization function that generates diagrams showing the interactions of
#' input genes in the PIN. See \code{\link{run_pathfindR}} for the wrapper
#' function of the pathfindR workflow. \code{\link[pathview]{pathview}} for
#' KEGG pathway-based data integration and visualization.
#'
#' @examples
#' \dontrun{
#' visualize_terms(result_df, input_processed)
#' visualize_terms(result_df, hsa_KEGG = FALSE, pin_name_path = "IntAct")
#' }
visualize_terms <- function(result_df, input_processed = NULL,
                            hsa_KEGG = TRUE, pin_name_path = "Biogrid") {

  ############ Argument Checks
  if (!is.logical(hsa_KEGG)) {
    stop("the argument `hsa_KEGG` must be either TRUE or FALSE")
  }

  if (hsa_KEGG & is.null(input_processed)) {
    stop("`input_processed` must be specified when `hsa_KEGG` is `TRUE`")
  }

  ############ Generate Diagrams
  if (hsa_KEGG) {
    ## Prepare input
    genes_df <- input_processed[, c("GENE", "CHANGE")]
    rownames(genes_df) <- genes_df$GENE
    genes_df <- genes_df[, -1, drop = FALSE]
    visualize_hsa_KEGG(pw_table = result_df, gene_data = genes_df)
  } else {
    visualize_term_interactions(result_df = result_df, pin_name_path = pin_name_path)
  }
}

#' Visualize Interactions of Genes Involved in the Given Enriched Terms
#'
#' @param result_df Data frame of enrichment results. Must-have columns
#' are: "Term_Description", "Up_regulated" and "Down_regulated"
#' @param pin_name_path Name of the chosen PIN or path/to/PIN.sif. If PIN name,
#'   must be one of c("Biogrid", "GeneMania", "IntAct", "KEGG"). If
#'   path/to/PIN.sif, the file must comply with the PIN specifications. Defaults
#'   to Biogrid.
#'
#' @return Creates PNG files visualizing the interactions of genes involved
#' in the given enriched terms (annotated in the \code{result_df}) in the PIN used
#' for enrichment analysis (specified by \code{pin_name_path}). The PNG files are
#' saved in the folder "term_visualizations" under the current working directory.
#'
#' @details The following steps are performed for the visualization of interactions
#' of genes involved for each enriched term: \enumerate{
#'   \item shortest paths between all affected genes are determined (via \code{\link[igraph]{igraph}})
#'   \item the nodes of all shortest paths are merged
#'   \item the PIN is subsetted using the merged nodes (genes)
#'   \item using the PIN subset, the graph showing the interactions is generated
#'   \item the final graph is visualized using \code{\link[igraph]{igraph}}, colored by changed
#'   status (if provided), and is saved as a PNG file.
#' }
#'
#' @export
#'
#' @seealso See \code{\link{visualize_terms}} for the wrapper function
#'   for creating enriched term diagrams. See \code{\link{run_pathfindR}} for the
#'   wrapper function of the pathfindR enrichment workflow.
#'
#' @examples
#' \dontrun{
#' visualize_term_interactions(result_df, pin_name_path = "IntAct")
#' }
visualize_term_interactions <- function(result_df, pin_name_path) {
  ############ Initial Steps
  ## fix naming issue
  result_df$Term_Description <- gsub("\\/", "-", result_df$Term_Description)

  ## load PIN
  pin_path <- return_pin_path(pin_name_path)
  pin <- utils::read.delim(
    file = pin_path,
    header = FALSE, stringsAsFactors = FALSE
  )
  pin$V2 <- NULL

  ## pin graph
  pin_g <- igraph::graph_from_data_frame(pin, directed = FALSE)

  ## Create visualization output directory
  if (!dir.exists("term_visualizations")) {
    dir.create("term_visualizations")
  }

  ############ Visualize interactions by enriched term
  for (i in base::seq_len(nrow(result_df))) {
    current_row <- result_df[i, ]

    up_genes <- unlist(strsplit(current_row$Up_regulated, ", "))
    down_genes <- unlist(strsplit(current_row$Down_regulated, ", "))
    current_genes <- c(down_genes, up_genes)

    ## Add active snw genes if listed
    if (!is.null(result_df$non_DEG_Active_Snw_Genes)) {
      snw_genes <- unlist(strsplit(current_row$non_DEG_Active_Snw_Genes, ", "))
      current_genes <- c(current_genes, snw_genes)
    } else {
      snw_genes <- NULL
    }

    if (length(current_genes) < 2) {
      message(paste0(
        "< 2 genes, skipping visualization of ",
        current_row$Term_Description
      ))
    } else {
      cat(paste0(
        "Visualizing: ", current_row$Term_Description,
        paste(rep(" ", 200), collapse = ""), "\r"
      ))

      ## Find genes without direct interaction
      cond1 <- pin$V1 %in% current_genes
      cond2 <- pin$V3 %in% current_genes
      direct_interactions <- pin[cond1 & cond2, ]
      tmp <- c(direct_interactions$V1, direct_interactions$V3)
      missing_genes <- current_genes[!current_genes %in% tmp]

      ## Find shortest path between genes without direct interaction
      # and other current_genes
      s_path_genes <- c()
      for (gene in missing_genes) {
        tmp <- suppressWarnings(igraph::shortest_paths(
          pin_g,
          from = which(names(igraph::V(pin_g)) == gene),
          to = which(names(igraph::V(pin_g)) %in% current_genes),
          output = "vpath"
        ))

        tmp <- unique(unlist(lapply(tmp$vpath, function(x) names(x))))
        s_path_genes <- unique(c(s_path_genes, tmp))
      }

      final_genes <- unique(c(current_genes, s_path_genes))
      cond1 <- pin$V1 %in% final_genes
      cond2 <- pin$V3 %in% final_genes
      final_interactions <- pin[cond1 & cond2, ]
      g <- igraph::graph_from_data_frame(final_interactions, directed = FALSE)

      cond1 <- names(igraph::V(g)) %in% up_genes
      cond2 <- names(igraph::V(g)) %in% down_genes
      cond3 <- names(igraph::V(g)) %in% snw_genes
      igraph::V(g)$color <- ifelse(cond1, "red",
        ifelse(cond2, "green",
          ifelse(cond3, "blue", "gray60")
        )
      )

      path_to_png <- file.path("term_visualizations",
                               paste0(current_row$Term_Description, ".png"))
      grDevices::png(path_to_png, width = 1039, height = 831)
      # Plot the tree object
      igraph::plot.igraph(
        g,
        layout = igraph::layout.fruchterman.reingold,
        edge.curved = FALSE,
        vertex.size = 10,
        vertex.label.dist = 0,
        vertex.label.color = "black",
        asp = FALSE,
        vertex.label.cex = 0.8,
        edge.width = 1.2,
        edge.arrow.mode = 0,
        main = paste(current_row$Term_Description,
                     "\n Involved Gene Interactions in", pin_name_path)
      )

      if (is.null(snw_genes))
        graphics::legend("topleft",
                         legend = c(
                           "Upregulated Input Genes",
                           "Downregulated Input Genes",
                           "Other"
                         ),
                         col = c("red", "green", "gray60"),
                         pch = 19, cex = 1.5, bty = "n"
        ) else {
          graphics::legend("topleft",
                           legend = c(
                             "Non-input Active Snw. Genes",
                             "Upregulated Input Genes",
                             "Downregulated Input Genes",
                             "Other"
                           ),
                           col = c("blue", "red", "green", "gray60"),
                           pch = 19, cex = 1.5, bty = "n"
          )
        }
      grDevices::dev.off()
    }
  }
}

#' Visualize Human KEGG Pathways
#'
#' @param pw_table Data frame of enrichment results. Must-have columns are: "ID" and
#'   "Term_Description".
#' @param gene_data Single column data frame containing change values (e.g.
#'   log(fold change) values) for significant genes. Row names are gene symbols.
#'
#' @return Creates visualizations of the enriched human KEGG pathways with
#'  the package \code{pathview} and saves them in the folder
#'  "term_visualizations" under the current working directory.
#'
#' @import pathview
#' @export
#' @seealso \code{\link[pathview]{pathview}} for KEGG pathway-based data integration
#'   and visualization. See \code{\link{visualize_terms}} for the wrapper function
#'   for creating enriched term diagrams. See \code{\link{run_pathfindR}} for the
#'   wrapper function of the pathfindR enrichment workflow.
#' @examples
#' \dontrun{
#' visualize_hsa_KEGG(pathway_table, gene_data)
#' }
visualize_hsa_KEGG <- function(pw_table, gene_data) {

  ## fix KEGG names such as "Glycolysis / Gluconeogenesis"
  pw_table$Term_Description <- gsub("\\/", "-", pw_table$Term_Description)

  ## Create visualization output directory
  dir.create("term_visualizations")
  setwd("term_visualizations")
  on.exit(setwd(".."))

  for (i in base::seq_len(nrow(pw_table))) {
    ## If all change values are 100 (if no change values supplied in input)
    if (all(gene_data[, 1] == 100)) {
      pathview::pathview(
        gene.data = gene_data,
        gene.idtype = "SYMBOL",
        pathway.id = pw_table$ID[i],
        species = "hsa",
        out.suffix = pw_table$Term_Description[i],
        keys.align = "y", kegg.native = TRUE,
        key.pos = "topright", same.layer = FALSE,
        discrete = list(gene = TRUE, cpd = TRUE),
        limit = list(gene = 100, cpd = 100),
        high = list("#65909A", "red"),
        node.sum = "mean",
        both.dirs = list(gene = FALSE, cpd = FALSE),
        new.signature = FALSE, plot.col.key = FALSE
      )
      ## If binary values supplied (must be ordered)
    } else if (length(unique(gene_data[, 1])) == 2) {
      ## change to -1 to 1
      vals <- sort(unique(gene_data[, 1]))
      gene_data[, 1] <- ifelse(gene_data[, 1] == vals[1], -1, 1)

      pathview::pathview(
        gene.data = gene_data,
        gene.idtype = "SYMBOL",
        pathway.id = pw_table$ID[i],
        species = "hsa",
        out.suffix = pw_table$Term_Description[i],
        keys.align = "y", kegg.native = TRUE,
        key.pos = "topright", same.layer = FALSE,
        discrete = list(gene = TRUE, cpd = TRUE),
        bins = list(gene = 2, cpd = 2),
        new.signature = FALSE
      )
      ## If continuous change values supplied (eg. logFC)
    } else {
      pathview::pathview(
        gene.data = gene_data,
        gene.idtype = "SYMBOL",
        pathway.id = pw_table$ID[i],
        species = "hsa",
        out.suffix = pw_table$Term_Description[i],
        keys.align = "y", kegg.native = TRUE,
        key.pos = "topright", same.layer = FALSE,
        new.signature = FALSE
      )
    }
  }
}


#' Plot the Bubble Chart of Enrichment Results
#'
#' This function is used to plot a bubble chart displaying the enrichment
#' results.
#'
#' @param result_df a data frame that must contain the following columns: \describe{
#'   \item{Term_Description}{Description of the enriched term}
#'   \item{Fold_Enrichment}{Fold enrichment value for the enriched term}
#'   \item{lowest_p}{the lowest adjusted-p value of the given term over all iterations}
#'   \item{Up_regulated}{the up-regulated genes in the input involved in the given term's gene set, comma-separated}
#'   \item{Down_regulated}{the down-regulated genes in the input involved in the given term's gene set, comma-separated}
#'   \item{Cluster(OPTIONAL)}{the cluster to which the enriched term is assigned}
#' }
#' @param plot_by_cluster boolean value indicating whether or not to group the
#'  enriched terms by cluster (works if \code{result_df} contains a
#'  "Cluster" column).
#' @param num_bubbles number of sizes displayed in the legend \code{# of DEGs}
#'  (Default = 4)
#' @param even_breaks whether or not to set even breaks for the number of sizes
#'  displayed in the legend \code{# of DEGs}. If \code{TRUE} (default), sets
#'  equal breaks and the number of displayed bubbles may be different than the
#'  number set by \code{num_bubbles}. If the exact number set by
#'  \code{num_bubbles} is required, set this argument to \code{FALSE}
#'
#' @return a \code{\link[ggplot2]{ggplot2}} object containing the bubble chart.
#' The x-axis corresponds to fold enrichment values while the y-axis indicates
#' the enriched terms. Size of the bubble indicates the number of DEGs in the
#' given enriched term. Color indicates the -log10(lowest-p) value. The closer
#' the color is to red, the more significant the enrichment is. Optionally, if
#' "Cluster" is a column of \code{result_df} and \code{plot_by_cluster == TRUE},
#' the enriched terms are grouped by clusters.
#'
#' @import ggplot2
#' @export
#'
#' @examples
#' g <- enrichment_chart(RA_output)
enrichment_chart <- function(result_df, plot_by_cluster = FALSE,
                             num_bubbles = 4, even_breaks = TRUE) {
  necessary <- c(
    "Term_Description", "Fold_Enrichment", "lowest_p",
    "Up_regulated", "Down_regulated"
  )
  if (!all(necessary %in% colnames(result_df))) {
    stop(
      "The input data frame must have the columns:\n",
      paste(necessary, collapse = ", ")
    )
  }

  if (!is.logical(plot_by_cluster)) {
    stop("`plot_by_cluster` must be either TRUE or FALSE")
  }

  # sort by lowest adj.p
  result_df <- result_df[order(result_df$lowest_p), ]

  num_genes <- vapply(
    result_df$Up_regulated,
    function(x) length(unlist(strsplit(x, ", "))), 1
  )
  num_genes <- num_genes + vapply(
    result_df$Down_regulated,
    function(x) length(unlist(strsplit(x, ", "))), 1
  )

  result_df$Term_Description <- factor(result_df$Term_Description,
    levels = rev(unique(result_df$Term_Description))
  )

  g <- ggplot2::ggplot(result_df, ggplot2::aes_(~Fold_Enrichment, ~Term_Description))
  g <- g + ggplot2::geom_point(ggplot2::aes(
    color = -log10(result_df$lowest_p),
    size = num_genes
  ), na.rm = TRUE)
  g <- g + ggplot2::theme_bw()
  g <- g + ggplot2::theme(
    axis.text.x = ggplot2::element_text(size = 10),
    axis.text.y = ggplot2::element_text(size = 10),
    plot.title = ggplot2::element_blank()
  )
  g <- g + ggplot2::xlab("Fold Enrichment")
  g <- g + ggplot2::theme(axis.title.y = ggplot2::element_blank())
  g <- g + ggplot2::labs(size = "# of DEGs", color = "-log10(lowest-p)")

  ## breaks for # of DEGs
  if (max(num_genes) < num_bubbles) {
    g <- g + ggplot2::scale_size_continuous(breaks = seq(0, max(num_genes)))
  } else {
    tmp1 <- base::seq(
      0, max(num_genes),
      round(max(num_genes) / (num_bubbles + 1))
    )
    tmp2 <- base::round(base::seq(0, max(num_genes),
      length.out = num_bubbles + 1
    ))
    brks <- ifelse(even_breaks, tmp1, tmp2)

    g <- g + ggplot2::scale_size_continuous(breaks = brks)
  }

  g <- g + ggplot2::scale_color_continuous(low = "#f5efef", high = "red")

  if (plot_by_cluster & "Cluster" %in% colnames(result_df)) {
    g <- g + ggplot2::facet_grid(result_df$Cluster ~ .,
      scales = "free_y", space = "free", drop = TRUE
    )
  } else if (plot_by_cluster) {
    message("For plotting by cluster, there must a column named `Cluster` in the input data frame!")
  }

  return(g)
}


#' Plot Term-Gene Graph
#'
#' @param result_df A dataframe of pathfindR results that must contain the following columns: \describe{
#'   \item{Term_Description}{Description of the enriched term (necessary if \code{use_description = TRUE})}
#'   \item{ID}{ID of the enriched term (necessary if \code{use_description = FALSE})}
#'   \item{lowest_p}{the lowest adjusted-p value of the given term over all iterations}
#'   \item{Up_regulated}{the up-regulated genes in the input involved in the given term's gene set, comma-separated}
#'   \item{Down_regulated}{the down-regulated genes in the input involved in the given term's gene set, comma-separated}
#' }
#' @param num_terms Number of top enriched terms to use while creating the graph. Set to \code{NULL} to use
#'  all enriched terms (default = 10, i.e. top 10 terms)
#' @param layout The type of layout to create (see \code{\link[ggraph]{ggraph}} for details. Default = "auto")
#' @param use_description Boolean argument to indicate whether term descriptions
#'  (in the "Term_Description" column) should be used. (default = \code{FALSE})
#' @param node_size Argument to indicate whether to use number of DEGs ("num_DEGS")
#'  or the -log10(lowest p value) ("p_val") for adjusting the node sizes (default = "num_DEGS")
#'
#' @return a  \code{\link[ggraph]{ggraph}} object containing the term-gene graph.
#'  Each node corresponds to an enriched term (beige), an up-regulated gene (green)
#'  or a down-regulated gene (red). An edge between a term and a gene indicates
#'  that the given term involves the gene. Size of a term node is proportional
#'  to either the number of genes (if \code{node_size = "num_DEGs"}) or
#'  the -log10(lowest p value) (if \code{node_size = "p_val"}).
#'
#' @details this plotting function was created based on the Gene-Concept
#' network visualization by the R package \code{enrichplot}.
#'
#' @import ggraph
#' @export
#'
#' @examples
#' p <- term_gene_graph(RA_output)
#' p <- term_gene_graph(RA_output, num_terms = 5)
#' p <- term_gene_graph(RA_output, node_size = "p_val")
term_gene_graph <- function(result_df, num_terms = 10,
                            layout = "auto", use_description = FALSE,
                            node_size = "num_DEGs") {

  ############ Argument Checks
  ### Check num_terms is NULL or numeric
  if (!is.numeric(num_terms) & !is.null(num_terms)) {
    stop("`num_terms` must either be numeric or NULL!")
  }

  ### Check use_description is boolean
  if (!is.logical(use_description)) {
    stop("`use_description` must either be TRUE or FALSE!")
  }

  ### Set column for term labels
  ID_column <- ifelse(use_description, "Term_Description", "ID")

  ### Check node_size
  node_size <- match.arg(node_size, c("num_DEGs", "p_val"))

  ### Check necessary columnns
  necessary_cols <- c("Up_regulated", "Down_regulated", "lowest_p", ID_column)

  if (!all(necessary_cols %in% colnames(result_df))) {
    stop(paste(c(
      "All of", paste(necessary_cols, collapse = ", "),
      "must be present in `results_df`!"
    ), collapse = " "))
  }

  ############ Initial steps
  ### set num_terms to NULL if number of enriched terms is smaller than num_terms
  if (!is.null(num_terms)) {
    if (nrow(result_df) < num_terms) {
      num_terms <- NULL
    }
  }

  ### Order and filter for top N genes
  result_df <- result_df[order(result_df$lowest_p, decreasing = FALSE), ]
  if (!is.null(num_terms)) {
    result_df <- result_df[1:num_terms, ]
  }

  ### Prep data frame for graph
  graph_df <- data.frame()
  for (i in base::seq_len(nrow(result_df))) {
    up_genes <- unlist(strsplit(result_df$Up_regulated[i], ", "))
    down_genes <- unlist(strsplit(result_df$Down_regulated[i], ", "))
    genes <- c(up_genes, down_genes)

    for (gene in genes) {
      graph_df <- rbind(graph_df,
                         data.frame(Term = result_df[i, ID_column],
                                    Gene = gene))
    }
  }


  up_genes <- lapply(result_df$Up_regulated,
                     function(x) unlist(strsplit(x, ", ")))
  up_genes <- unlist(up_genes)

  ############ Create graph object and plot
  ### create igraph object
  g <- igraph::graph_from_data_frame(graph_df, directed = FALSE)
  cond_term <- names(igraph::V(g)) %in% result_df[, ID_column]
  cond_up_gene <- names(igraph::V(g)) %in% up_genes
  igraph::V(g)$type <- ifelse(cond_term, "term",
                              ifelse(cond_up_gene, "up", "down"))
  # Adjust node sizes
  if (node_size == "num_DEGs") {
    sizes <- igraph::degree(g)
    sizes <- ifelse(igraph::V(g)$type == "term", sizes, 2)
    size_label <- "# DEGs"
  } else {
    idx <- match(names(igraph::V(g)), result_df[, ID_column])
    sizes <- -log10(result_df$lowest_p[idx])
    sizes[is.na(sizes)] <- 2
    size_label <- "-log10(p)"
  }
  igraph::V(g)$size <- sizes
  igraph::V(g)$label.cex <- 0.5
  igraph::V(g)$frame.color <- "gray"
  igraph::V(g)$color <- ifelse(igraph::V(g)$type == "term", "#E5D7BF",
    ifelse(igraph::V(g)$type == "up", "green",
      "red"
    )
  )

  ### Create graph
  p <- ggraph::ggraph(g, layout = layout)
  p <- p + ggraph::geom_edge_link(alpha = .8, colour = "darkgrey")
  p <- p + ggraph::geom_node_point(ggplot2::aes_(color = ~ I(color), size = ~size))
  p <- p + ggplot2::scale_size(
    range = c(5, 10),
    breaks = round(seq(round(min(igraph::V(g)$size)),
      round(max(igraph::V(g)$size)),
      length.out = 4
    )),
    name = size_label
  )
  p <- p + ggplot2::theme_void()
  p <- p + ggraph::geom_node_text(ggplot2::aes_(label = ~name), nudge_y = .2)
  p <- p + ggplot2::scale_colour_manual(
    values = unique(igraph::V(g)$color),
    name = NULL,
    labels = c(
      "enriched term",
      "up-regulated gene",
      "down-regulated gene"
    )
  )
  if (is.null(num_terms)) {
    p <- p + ggplot2::ggtitle("Term-Gene Graph")
  } else {
    p <- p + ggplot2::ggtitle("Term-Gene Graph",
      subtitle = paste(c("Top", num_terms, "terms"),
        collapse = " "
      )
    )
  }

  p <- p + ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5),
    plot.subtitle = ggplot2::element_text(hjust = 0.5)
  )

  return(p)
}
