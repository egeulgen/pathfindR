#' Create Pathway Diagrams
#'
#' @param result_df Data frame of enrichment results. Must-have columns for
#' KEGG human pathway diagrams are: "ID" and "Pathway".
#' Must-have columns for the rest are: "Pathway", "Up_regulated" and "Down_regulated"
#' @param input_processed input data processed via `input_processing`, not necessary for
#' visualizations other than KEGG human pathway diagrams.
#' @param gene_sets the gene sets used for enrichment analysis. Possible gene sets
#'  are KEGG or non-KEGG (Default = "KEGG").
#' @param pin_name_path Name of the chosen PIN or path/to/PIN.sif. If PIN name,
#'   must be one of c("Biogrid", "GeneMania", "IntAct", "KEGG"). If
#'   path/to/PIN.sif, the file must comply with the PIN specifications. Defaults
#'   to Biogrid.
#'
#' @return Depending on `gene_sets`, creates visualization of interactions of genes
#' involved in the list of pathways in `result_df` and saves them in the folder
#' "pathway_visualizations" under the current working directory.
#'
#'
#' @details For `gene_sets == "KEGG"`, KEGG human pathway diagrams are created,
#' affected nodes colored by up/down regulation status.
#' For other `gene_sets`, interactions of affected genes are determined (via a shortest-path
#' algorithm) and are visualized (colored by change status) using igraph.
#'
#'
#' @export
#'
#' @seealso See \code{\link{visualize_hsa_KEGG}} for the visualization function
#' of human KEGG diagrams. See \code{\link{visualize_pw_interactions}} for the
#' visualization function that generates diagrams showing the interactions of
#' input genes in the PIN. See \code{\link{run_pathfindR}} for the wrapper
#' function of the pathfindR workflow. \code{\link[pathview]{pathview}} for
#' KEGG pathway-based data integration and visualization.
#'
#' @examples
#' \dontrun{
#' visualize_pws(result_df, input_processed)
#' visualize_pws(result_df, gene_sets = "GO-BP", pin_name_path = "IntAct")
#' }
visualize_pws <- function(result_df, input_processed = NULL,
                          gene_sets = "KEGG", pin_name_path = "Biogrid") {
  ############ Argument Check
  if (gene_sets == "KEGG" & is.null(input_processed))
    stop("`input_processed` must be specified when `gene_sets` is KEGG")

  ############ Generate pathway diagrams
  if (gene_sets == "KEGG") {
    ## Prepare input
    genes_df <- input_processed[, c("GENE", "CHANGE")]
    rownames(genes_df) <- genes_df$GENE
    genes_df <- genes_df[, -1, drop = FALSE]
    visualize_hsa_KEGG(result_df, genes_df)

  } else {
    visualize_pw_interactions(result_df, pin_name_path)
  }
}


#' Visualize Interactions of Genes Involved in the Given Pathways
#'
#' @param result_df Data frame of enrichment results. Must-have columns
#' are: "Pathway", "Up_regulated" and "Down_regulated"
#' @param pin_name_path Name of the chosen PIN or path/to/PIN.sif. If PIN name,
#'   must be one of c("Biogrid", "GeneMania", "IntAct", "KEGG"). If
#'   path/to/PIN.sif, the file must comply with the PIN specifications. Defaults
#'   to Biogrid.
#'
#' @return Creates PNG files visualizing the interactions of genes involved
#' in the given pathways (annotated in the `result_df`) in the PIN used for enrichment
#' analysis (specified by `pin_name_path`). The PNG files are saved in the folder
#' "pathway_visualizations" under the current working directory.
#'
#' @details The following steps are performed for the visualization of interactions
#' of genes involved in the given pathways: \enumerate{
#'   \item shortest paths between all affected genes are determined (via `igraph`)
#'   \item the nodes of all shortest pathways are merged
#'   \item the PIN is subsetted using the merged nodes (genes)
#'   \item using the PIN subset, the graph showing the interactions is generated
#'   \item the final graph is visualized using `igraph`, colored by changed status and saved as a PNG file.
#'}
#'
#' @export
#'
#' @seealso See \code{\link{visualize_pws}} for the wrapper function
#'   for creating pathway diagrams. See \code{\link{run_pathfindR}} for the
#'   wrapper function of the pathfindR workflow.
#'
#' @examples
#' \dontrun{
#' visualize_pw_interactions(result_df, pin_name_path = "IntAct")
#' }
visualize_pw_interactions <- function(result_df, pin_name_path) {
  ############ Initial Steps
  ## fix pathway naming issue
  result_df$Pathway <- gsub("\\/", "-", result_df$Pathway)

  ## load PIN
  pin_path <- return_pin_path(pin_name_path)
  pin <- utils::read.delim(file = pin_path,
                           header = FALSE, stringsAsFactors = FALSE)
  pin$V2 <- NULL

  ## pin graph
  pin_g <- igraph::graph_from_data_frame(pin, directed = FALSE)

  ## Create visualization output directory
  if(!dir.exists("pathway_visualizations"))
    dir.create("pathway_visualizations")

    ############ visualize interactions by pathway
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

    if (length(current_genes) < 2)
      message(paste0("< 2 genes, skipping visualization of ",
                     current_row$Pathway))
    else {
      cat(paste0("Visualizing: ", current_row$Pathway,
                 paste(rep(" ", 200), collapse = ""), "\r"))

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
          output = "vpath"))

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
                                          ifelse(cond3, "blue", "gray60")))

      path_to_png <- file.path("pathway_visualizations",
                               paste0(current_row$Pathway, ".png"))
      grDevices::png(path_to_png, width = 1039, height = 831)
      #Plot the tree object
      igraph::plot.igraph(
        g,
        layout=igraph::layout.fruchterman.reingold,
        edge.curved=FALSE,
        vertex.size=10,
        vertex.label.dist=0,
        vertex.label.color="black",
        asp=FALSE,
        vertex.label.cex=0.8,
        edge.width=1.2,
        edge.arrow.mode=0,
        main=paste(current_row$Pathway,
                   "\n Involved Gene Interactions in", pin_name_path)
      )

      graphics::legend("topleft", legend = c("Non-input Active Snw. Genes",
                                             "Upregulated Input Genes",
                                             "Downregulated Input Genes",
                                             "Other"),
                       col = c("blue", "red", "green", "gray60"),
                       pch = 19, cex = 1.5, bty = "n")
      grDevices::dev.off()
    }
  }
}

#' Visualize Human KEGG Pathways
#'
#' @param pw_table Data frame of enrichment results. Must-have columns are: "ID" and
#'   "Pathway".
#' @param gene_data Single column data frame containing change values (e.g.
#'   log(fold change) values) for significant genes. Row names are gene symbols.
#'
#' @return Creates visualizations of the pathways with the package \code{pathview}
#'  and saves them in the folder "pathway_visualizations" under the current working directory.
#'
#' @import pathview
#' @export
#' @seealso \code{\link[pathview]{pathview}} for KEGG pathway-based data integration
#'   and visualization. See \code{\link{visualize_pws}} for the wrapper function
#'   for creating pathway diagrams. See \code{\link{run_pathfindR}} for the
#'   wrapper function of the pathfindR workflow.
#' @examples
#' \dontrun{
#' visualize_hsa_KEGG(pathway_table, gene_data)
#' }
visualize_hsa_KEGG <- function(pw_table, gene_data) {

  ## fix KEGG names such as "Glycolysis / Gluconeogenesis"
  pw_table$Pathway <- gsub("\\/", "-", pw_table$Pathway)

  ## Create visualization output directory
  dir.create("pathway_visualizations")
  setwd("pathway_visualizations")
  on.exit(setwd(".."))

  for (i in base::seq_len(nrow(pw_table))) {
    ## If all change values are 100 (if no change values supplied in input)
    if (all(gene_data[, 1] == 100)) {
      pathview::pathview(gene.data = gene_data,
                         gene.idtype = "SYMBOL",
                         pathway.id = pw_table$ID[i],
                         species = "hsa",
                         out.suffix = pw_table$Pathway[i],
                         keys.align = "y", kegg.native = TRUE,
                         key.pos = "topright", same.layer = FALSE,
                         discrete = list(gene=TRUE, cpd=TRUE),
                         limit = list(gene=100, cpd=100),
                         high = list("#65909A", "red"),
                         node.sum = "mean",
                         both.dirs = list(gene=FALSE, cpd=FALSE),
                         new.signature = FALSE, plot.col.key = FALSE)
      ## If binary values supplied (must be ordered)
    } else if (length(unique(gene_data[, 1])) == 2) {
      ## change to -1 to 1
      vals <- sort(unique(gene_data[, 1]))
      gene_data[, 1] <- ifelse(gene_data[, 1] == vals[1], -1, 1)

      pathview::pathview(gene.data = gene_data,
                         gene.idtype = "SYMBOL",
                         pathway.id = pw_table$ID[i],
                         species = "hsa",
                         out.suffix = pw_table$Pathway[i],
                         keys.align = "y", kegg.native = TRUE,
                         key.pos = "topright", same.layer = FALSE,
                         discrete = list(gene=TRUE, cpd=TRUE),
                         bins = list(gene=2, cpd = 2),
                         new.signature = FALSE)
      ## If continuous change values supplied (eg. logFC)
    } else {
      pathview::pathview(gene.data = gene_data,
                         gene.idtype = "SYMBOL",
                         pathway.id = pw_table$ID[i],
                         species = "hsa",
                         out.suffix = pw_table$Pathway[i],
                         keys.align = "y", kegg.native = TRUE,
                         key.pos = "topright", same.layer = FALSE,
                         new.signature = FALSE)
    }
  }
}


#' Plot the Bubble Chart of Enrichment Results
#'
#' This function is used to plot a bubble chart displaying the enrichment
#' results.
#'
#' @param result_df a data frame that must contain the following columns:\describe{
#'   \item{Pathway}{Description of the enriched pathway}
#'   \item{Fold_Enrichment}{Fold enrichment value for the enriched pathway}
#'   \item{lowest_p}{the lowest adjusted-p value of the given pathway over all iterations}
#'   \item{Up_regulated}{the up-regulated genes in the input involved in the given pathway, comma-separated}
#'   \item{Down_regulated}{the down-regulated genes in the input involved in the given pathway, comma-separated}
#'   \item{Cluster(OPTIONAL)}{the cluster to which the pathway is assigned}
#' }
#' @param plot_by_cluster boolean value indicating whether or not to group the
#' pathways by cluster (works if "Cluster" is a column of `result_df`).
#' @param num_bubbles number of sizes displayed in the legend `# of DEGs` (Default = 5).
#' @param even_breaks whether or not to set even breaks for the number of sizes displayed
#' in the legend `# of DEGs`. If `TRUE` (default), sets equal breaks and the number of
#' displayed bubbles may be different than the number set by `num_bubbles`. If the exact
#' number set by `num_bubbles` is required, set this argument to `FALSE`.
#'
#' @return a `ggplot2` object containing the bubble chart. The x-axis corresponds to
#' fold enrichment values while the y-axis indicates the enriched pathways. Size of
#' the bubble indicates the number of DEGs in the given pathway. Color indicates
#' the -log10(lowest-p) value. The closer the color is to red, the more significant
#' the enrichment is. Optionally, if "Cluster" is a column of `result_df` and
#' plot_by_cluster == TRUE, the pathways are grouped by clusters.
#'
#' @import ggplot2
#' @export
#'
#' @examples
#' g <- enrichment_chart(RA_output)
enrichment_chart <- function(result_df, plot_by_cluster = FALSE,
                             num_bubbles = 4, even_breaks = TRUE) {

  necessary <- c("Pathway", "Fold_Enrichment", "lowest_p",
                 "Up_regulated", "Down_regulated")
  if (!all(necessary %in% colnames(result_df)))
    stop("The input data frame must have the columns:
         Pathway, Fold_Enrichment, lowest_p, Up_regulated, Down_regulated")

  if (!is.logical(plot_by_cluster))
    stop("plot_by_cluster must be either TRUE or FALSE")

  # sort by lowest adj.p
  result_df <- result_df[order(result_df$lowest_p), ]

  n <- vapply(result_df$Up_regulated,
              function(x) length(unlist(strsplit(x, ", "))), 1)
  n <- n + vapply(result_df$Down_regulated,
                  function(x) length(unlist(strsplit(x, ", "))), 1)

  result_df$Pathway <- factor(result_df$Pathway,
                              levels = rev(unique(result_df$Pathway)))

  g <- ggplot2::ggplot(result_df, ggplot2::aes_(~Fold_Enrichment, ~Pathway))
  g <- g + ggplot2::geom_point(ggplot2::aes(color = -log10(result_df$lowest_p),
                                            size = n), na.rm = TRUE)
  g <- g + ggplot2::theme_bw()
  g <- g + ggplot2::theme(axis.text.x = ggplot2::element_text(size = 10),
                          axis.text.y = ggplot2::element_text(size = 10),
                          plot.title = ggplot2::element_blank())
  g <- g + ggplot2::xlab("Fold Enrichment") + ggplot2::ylab("")
  g <- g + ggplot2::labs(size = "# of DEGs", color = "-log10(lowest-p)")

  ## breaks for # of DEGs
  if (max(n) < num_bubbles) {
    g <- g + ggplot2::scale_size_continuous(breaks = seq(0, max(n)))
  } else {
    tmp1 <- base::seq(0, max(n), round(max(n) / (num_bubbles + 1)))
    tmp2 <- base::round(base::seq(0, max(n), length.out = num_bubbles + 1))
    brks <- ifelse(even_breaks, tmp1, tmp2)

    g <- g + ggplot2::scale_size_continuous(breaks = brks)
  }

  g <- g + ggplot2::scale_color_continuous(low = "#f5efef", high = "red")

  if (plot_by_cluster & "Cluster" %in% colnames(result_df)) {
    g <- g + ggplot2::facet_grid(result_df$Cluster~.,
                                 scales = "free_y", space = "free", drop = TRUE)
  } else if (plot_by_cluster) {
    message("For plotting by cluster, there must a column named `Cluster`
            in the input data frame!")
  }

  return(g)
}


#' Plot Term-Gene Graph
#'
#' @param result_df A dataframe of pathfindR results that must contain the following columns:\describe{
#'   \item{Pathway}{Description of the enriched pathway (necessary if `use_names` is TRUE)}
#'   \item{ID}{ID of the enriched pathway (necessary if `use_names` is FALSE)}
#'   \item{lowest_p}{the lowest adjusted-p value of the given pathway over all iterations}
#'   \item{Up_regulated}{the up-regulated genes in the input involved in the given pathway, comma-separated}
#'   \item{Down_regulated}{the down-regulated genes in the input involved in the given pathway, comma-separated}
#' }
#' @param num_terms Number of top terms to use while creating the graph. Set to `NULL` to use
#' all terms (default = 10, i.e. top 10 terms)
#' @param layout The type of layout to create (see \code{\link[ggraph]{ggraph}} for details. Default = "auto")
#' @param use_names Boolean argument to indicate whether term descriptions (in the Pathway column) should be used. (default = `FALSE`)
#' @param node_size Argument to indicate whether to use number of DEGs ("num_DEGS")
#' or the -log10(lowest p value) ("p_val") for adjusting the node sizes (default = "num_DEGS")
#'
#' @return a `ggplot` object containing the term-gene graph. Each node corresponds to
#' an enriched term (beige), an up-regulated gene (green) or a down-regulated gene (red).
#' An edge between a term and a gene indicates that the given term involves the gene. Size
#' of a term node is proportional to either the number of genes (if `node_size`= "num_DEGs")
#' or the -log10(lowest p value) (if `node_size`= "p_val").
#'
#' @details this plotting function was created based on the Gene-Concept network visualization
#' by the R package `enrichplot`.
#'
#' @import ggraph
#' @export
#'
#' @examples
#' p <- term_gene_graph(RA_output)
#' p <- term_gene_graph(RA_output, num_terms = 5)
#' p <- term_gene_graph(RA_output, node_size = "p_val")
term_gene_graph <- function(result_df, num_terms = 10,
                            layout = "auto", use_names = FALSE,
                            node_size = "num_DEGs") {

  ############ Argument Checks
  ### Set column for term labels
  ID_column <- ifelse(use_names, "Pathway", "ID")

  ### Check num_terms is NULL or numeric
  if (!is.numeric(num_terms) & !is.null(num_terms))
    stop("`num_terms` must either be numeric or NULL!")

  ### Check use_names is boolean
  if (!is.logical(use_names))
    stop("`use_names` must either be TRUE or FALSE!")

  ### Check node_size
  if (node_size != "num_DEGs" & node_size != "p_val")
    stop("`node_size` must either be num_DEGs or p_val!")

  ### Check necessary columnns
  necessary_cols <- c("Up_regulated", "Down_regulated", "lowest_p", ID_column)

  if (!all(necessary_cols %in% colnames(result_df)))
    stop(paste(c("All of", paste(necessary_cols, collapse = ", "),
                 "must be present in `results_df`!"), collapse = " "))

  ############ Initial steps
  ### Order and filter for top N genes
  df_for_vis <- result_df[order(result_df$lowest_p, decreasing = FALSE), ]
  if (!is.null(num_terms))
    df_for_vis <- df_for_vis[1:num_terms, ]

  ### Prep data frame for graph
  for_graph <- data.frame()
  for (i in base::seq_len(nrow(df_for_vis))) {
    up_genes <- unlist(strsplit(df_for_vis$Up_regulated[i], ", "))
    down_genes <- unlist(strsplit(df_for_vis$Down_regulated[i], ", "))
    genes <- c(up_genes, down_genes)

    for (gene in genes) {
      for_graph <- rbind(for_graph,
                         data.frame(Pathway = df_for_vis[i, ID_column],
                                    Gene = gene))
    }
  }

  up_genes <- unlist(lapply(df_for_vis$Up_regulated,
                            function(x) unlist(strsplit(x, ", "))))

  ############ Create graph and plot

  ### create igraph object
  g <- igraph::graph_from_data_frame(for_graph)
  cond1 <- names(igraph::V(g)) %in% df_for_vis[, ID_column]
  cond2 <- names(igraph::V(g)) %in% up_genes
  igraph::V(g)$type <- ifelse(cond1, "pathway",
                              ifelse(cond2, "up", "down"))
  # Adjust node sizes
  if (node_size == "num_DEGs") {
    sizes <- igraph::degree(g)
    sizes <- ifelse(igraph::V(g)$type == "pathway", sizes, 2)
    size_label <- "# DEGs"
  } else {
    idx <- match(names(igraph::V(g)), df_for_vis[, ID_column])
    sizes <- -log10(df_for_vis$lowest_p[idx])
    sizes[is.na(sizes)] <- 2
    size_label <- "-log10(p)"
  }
  igraph::V(g)$size <- sizes
  igraph::V(g)$label.cex <- 0.5
  igraph::V(g)$frame.color <-  "gray"
  igraph::V(g)$color <- ifelse(igraph::V(g)$type == "pathway", "#E5D7BF",
                               ifelse(igraph::V(g)$type == "up", "green",
                                      "red"))

  ### Create graph
  p <- ggraph::ggraph(g, layout = layout)
  p <- p + ggraph::geom_edge_link(alpha=.8, colour='darkgrey')
  p <- p + ggraph::geom_node_point(ggplot2::aes_(color=~I(color), size=~size))
  p <- p + ggplot2::scale_size(range=c(5, 10),
                               breaks=round(seq(round(min(igraph::V(g)$size)),
                                                round(max(igraph::V(g)$size)),
                                                length.out = 4)),
                               name = size_label)
  p <- p + ggplot2::theme_void()
  p <- p + ggraph::geom_node_text(ggplot2::aes_(label=~name), nudge_y = .2)
  p <- p + ggplot2::scale_colour_manual(values = unique(igraph::V(g)$color),
                                        name = NULL,
                                        labels = c("enriched term",
                                                   "up-regulated gene",
                                                   "down-regulated gene"))
  if (is.null(num_terms)) {
    p <- p + ggplot2::ggtitle("Term-Gene Graph")
  } else {
    p <- p + ggplot2::ggtitle("Term-Gene Graph",
                              subtitle = paste(c("Top", num_terms, "terms"),
                                               collapse = " "))
  }

  p <- p + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                          plot.subtitle = ggplot2::element_text(hjust = 0.5))

  return(p)

}
