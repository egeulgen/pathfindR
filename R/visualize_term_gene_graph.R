#' Create Term-Gene Graph
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
#' @param layout The type of layout to create (see \code{\link[ggraph]{ggraph}} for details. Default = 'stress')
#' @param use_description Boolean argument to indicate whether term descriptions
#'  (in the 'Term_Description' column) should be used. (default = \code{FALSE})
#' @param node_size Argument to indicate whether to use number of significant genes ('num_genes')
#'  or the -log10(lowest p value) ('p_val') for adjusting the node sizes (default = 'num_genes')
#' @param node_colors vector of 3 colors to be used for coloring nodes (colors for term nodes, up, and down, respectively)
#'
#' @return a  \code{\link[ggraph]{ggraph}} object containing the term-gene graph.
#'  Each node corresponds to an enriched term (beige), an up-regulated gene (green)
#'  or a down-regulated gene (red). An edge between a term and a gene indicates
#'  that the given term involves the gene. Size of a term node is proportional
#'  to either the number of genes (if \code{node_size = 'num_genes'}) or
#'  the -log10(lowest p value) (if \code{node_size = 'p_val'}).
#'
#' @details This function (adapted from the Gene-Concept network visualization
#' by the R package \code{enrichplot}) can be utilized to visualize which input
#' genes are involved in the enriched terms as a graph. The term-gene graph
#' shows the links between genes and biological terms and allows for the
#' investigation of multiple terms to which significant genes are related. The
#' graph also enables determination of the overlap between the enriched terms
#' by identifying shared and distinct significant term-related genes.
#'
#' @import ggraph
#' @export
#'
#' @examples
#' p <- term_gene_graph(example_pathfindR_output)
#' p <- term_gene_graph(example_pathfindR_output, num_terms = 5)
#' p <- term_gene_graph(example_pathfindR_output, node_size = "p_val")
term_gene_graph <- function(result_df, num_terms = 10, layout = "stress", use_description = FALSE,
                            node_size = "num_genes", node_colors = c("#E5D7BF", "green", "red")) {
  ############ Argument Checks Check num_terms is NULL or numeric
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
  val_node_size <- c("num_genes", "p_val")
  if (!node_size %in% val_node_size) {
    stop("`node_size` should be one of ", paste(dQuote(val_node_size), collapse = ", "))
  }

  if (!is.data.frame(result_df)) {
    stop("`result_df` should be a data frame")
  }

  ### Check necessary columnns
  necessary_cols <- c(ID_column, "lowest_p", "Up_regulated", "Down_regulated")

  if (!all(necessary_cols %in% colnames(result_df))) {
    stop(paste(c("All of", paste(necessary_cols, collapse = ", "), "must be present in `results_df`!"),
      collapse = " "
    ))
  }

  if (!is.atomic(node_colors)) {
    stop("`node_colors` should be a vector of colors")
  }

  if (!all(vapply(node_colors, isColor, TRUE))) {
    stop("`node_colors` should be a vector of valid colors")
  }

  if (length(node_colors) != 3) {
    stop("`node_colors` must contain exactly 3 colors")
  }

  ############ Initial steps set num_terms to NULL if number of enriched
  ############ terms is smaller than num_terms
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
    for (gene in c(up_genes, down_genes)) {
      graph_df <- rbind(graph_df, data.frame(
        Term = result_df[i, ID_column],
        Gene = gene
      ))
    }
  }

  up_genes <- lapply(result_df$Up_regulated, function(x) unlist(strsplit(x, ", ")))
  up_genes <- unlist(up_genes)

  ############ Create graph object and plot create igraph object
  g <- igraph::graph_from_data_frame(graph_df, directed = FALSE)
  cond_term <- names(igraph::V(g)) %in% result_df[, ID_column]
  cond_up_gene <- names(igraph::V(g)) %in% up_genes

  node_type <- ifelse(cond_term, "term", ifelse(cond_up_gene, "up", "down"))
  node_type <- factor(node_type, levels = c("term", "up", "down"))
  node_type <- droplevels(node_type)
  igraph::V(g)$type <- node_type

  type_descriptions <- c(term = "enriched term", up = "up-regulated gene", down = "down-regulated gene")
  type_descriptions <- type_descriptions[levels(node_type)]

  names(node_colors) <- c("term", "up", "down")
  node_colors <- node_colors[levels(node_type)]

  # Adjust node sizes
  if (node_size == "num_genes") {
    sizes <- igraph::degree(g)
    sizes <- ifelse(igraph::V(g)$type == "term", sizes, 2)
    size_label <- "# genes"
  } else {
    idx <- match(names(igraph::V(g)), result_df[, ID_column])
    sizes <- -log10(result_df$lowest_p[idx])
    sizes[is.na(sizes)] <- 2
    size_label <- "-log10(p)"
  }
  igraph::V(g)$size <- sizes
  igraph::V(g)$label.cex <- 0.5
  igraph::V(g)$frame.color <- "gray"

  ### Create graph
  p <- ggraph::ggraph(g, layout = layout)
  p <- p + ggraph::geom_edge_link(alpha = 0.8, colour = "darkgrey")
  p <- p + ggraph::geom_node_point(ggplot2::aes(color = .data$type, size = .data$size))
  p <- p + ggplot2::scale_size(range = c(5, 10), breaks = round(seq(round(min(igraph::V(g)$size)),
    round(max(igraph::V(g)$size)),
    length.out = 4
  )), name = size_label)
  p <- p + ggplot2::theme_void()
  p <- p + suppressWarnings(ggraph::geom_node_text(ggplot2::aes(label = .data$name),
    nudge_y = 0.2, repel = TRUE, max.overlaps = 20
  ))
  p <- p + ggplot2::scale_color_manual(
    values = node_colors, name = NULL,
    labels = type_descriptions
  )
  if (is.null(num_terms)) {
    p <- p + ggplot2::ggtitle("Term-Gene Graph")
  } else {
    p <- p + ggplot2::ggtitle("Term-Gene Graph", subtitle = paste(c(
      "Top", num_terms,
      "terms"
    ), collapse = " "))
  }

  p <- p + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), plot.subtitle = ggplot2::element_text(hjust = 0.5))

  return(p)
}
