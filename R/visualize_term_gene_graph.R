#' Create Term-Gene Graph
#'
#' @param result_df A dataframe of pathfindR results that must contain the following columns: \describe{
#'   \item{Term_Description}{Description of the enriched term (necessary if \code{use_description = TRUE})}
#'   \item{ID}{ID of the enriched term (necessary if \code{use_description = FALSE})}
#'   \item{lowest_p}{the lowest adjusted-p value of the given term over all iterations}
#'   \item{Up_regulated}{the up-regulated genes in the input involved in the given term's gene set, comma-separated}
#'   \item{Down_regulated}{the down-regulated genes in the input involved in the given term's gene set, comma-separated}
#' }
#' @param genes_df (optional) the input data that was used with \code{\link{run_pathfindR}} (default: \code{NULL}).
#'   It must be a data frame with at least 2 columns: \enumerate{
#'   \item Gene.Symbol (required)
#'   \item logFC (required)
#' }
#' @param order_by Argument to order the `result_df`, this influences the `num_terms` displayed (default: \code{'lowest_p'}).
#' @param term_size Argument to indicate whether to use number of significant genes ('num_genes')
#' @param term_fill Argument to indicate by what column to fill the term nodes (e.g. \code{term_fill = "Fold_Enrichment"}) (default: \code{NULL}).
#' @param num_terms Number of top enriched terms to use while creating the graph. Set to \code{NULL} to use
#'  all enriched terms (default = 10, i.e. top 10 terms)
#' @param use_description Boolean argument to indicate whether term descriptions
#'  (in the 'Term_Description' column) should be used. (default: \code{FALSE})
#' @param use_edge_weights Boolean argument to indicate whether genes are weighted by their term interactions, similar to an Up-Set plot but in graph context (default = \code{FALSE}).
#'  or the -log10(lowest p value) ('p_val') for adjusting the term node sizes (default: \code{'num_genes'})
#' @return A \link[igraph]{igraph} object
#'
#' @details
#' This function constructs an \link[igraph]{igraph} object from pathfindR output, creating a network that connects enriched biological terms to their involved genes.
#' By default, the graph connects term nodes to up-regulated genes and down-regulated genes.
#' The size of term nodes can be adjusted by either the number of significant genes (`term_size = 'num_genes'`) or by the statistical significance
#' (`term_size = 'p_val'`, using -log10(lowest p value)).
#'
#' When `genes_df` is provided, gene nodes contain values and not mere up/down binary values, allowing visualization of
#' expression direction and magnitude. When `term_fill` is supplied, term nodes obtain values enabling simultaneous visualization of pathway enrichment strength.
#'
#' Setting `use_edge_weights = TRUE` highlights hub genes by weighting edges based on how many terms a gene participates in, similar to an Up-Set plot but in
#' a graph context. The `num_terms` parameter controls how many top enriched terms are included (default: top 10), and `order_by` determines the ordering
#' criterion for term selection. The resulting igraph object can be visualized using \link[pathfindR]{create_term_gene_plot}.
#'
#' @import ggraph
#' @export
#'
#' @examples
#' # Normal gene-term with up/down regulated genes
#' g <- create_term_gene_graph(
#'   result_df = example_pathfindR_output
#' )
#' g <- create_term_gene_graph(
#'   result_df = example_pathfindR_output,
#'   num_terms = 5
#' )
#' g <- create_term_gene_graph(
#'   result_df = example_pathfindR_output,
#'   term_size = "p_val"
#' )
#'
#' # Coloring the term nodes
#' g <- create_term_gene_graph(
#'   result_df = example_pathfindR_output,
#'   term_fill = "Fold_Enrichment"
#' )
#'
#' # Adding edge weights
#' g <- create_term_gene_graph(
#'   result_df = example_pathfindR_output,
#'   term_fill = "Fold_Enrichment",
#'   use_edge_weights = TRUE
#' )
create_term_gene_graph <- function(
  result_df,
  genes_df = NULL,
  order_by = "lowest_p",
  term_size = "num_genes",
  term_fill = NULL,
  num_terms = 10,
  use_description = FALSE,
  use_edge_weights = FALSE
) {
  ## Error handling
  #--------------------------------------------------------------------#

  if (!is.data.frame(result_df) | inherits(result_df, "data.table")) {
    stop("`result_df` should be a data.frame!")
  }

  ### Check necessary columnns
  ID_column <- ifelse(use_description, "Term_Description", "ID")
  necessary_cols <- c(ID_column, "lowest_p", "Up_regulated", "Down_regulated")
  if (!all(necessary_cols %in% colnames(result_df))) {
    stop(paste(c("All of", paste(necessary_cols, collapse = ", "), "must be present in `results_df`!"),
      collapse = " "
    ))
  }

  if (!is.null(genes_df)) {
    if (!is.data.frame(genes_df) | inherits(genes_df, "data.table")) {
      stop("`genes_df` should be a data.frame!")
    }
    required_cols <- c("Gene.symbol", "logFC")
    if (!all(required_cols %in% colnames(genes_df))) {
      stop(paste(c("All of", paste(required_cols, collapse = ", "), "must be present in `genes_df`!"),
        collapse = " "
      ))
    }
  }

  result_df <- order_df_by_columnn(result_df, order_by)

  val_term_size <- c("num_genes", "p_val")
  if (!term_size %in% val_term_size) {
    stop("`term_size` should be one of ", paste(dQuote(val_term_size), collapse = ", "))
  }

  if (!is.null(term_fill)) {
    if (!c(term_fill %in% colnames(result_df))) {
      stop("`term_fill` is not found in the supplied `result_df`!")
    }
  }

  if (!is.null(num_terms)) {
    if (!is.numeric(num_terms) & !is.null(num_terms)) {
      stop("`num_terms` must either be numeric or NULL!")
    }

    n_rows <- nrow(result_df)
    if (nrow(result_df) < num_terms) {
      num_terms <- n_rows
      cat("`num_terms` is less than n rows of `result_df`, setting `num_terms` to ", num_terms)
    }
    result_df <- result_df[1:num_terms, ]
  }

  if (!is.logical(use_description)) {
    stop("`use_description` must either be TRUE or FALSE!")
  }

  if (!is.logical(use_edge_weights)) {
    stop("`use_edge_weights` must either be TRUE or FALSE!")
  }

  ## MAIN
  #--------------------------------------------------------------------#

  ### Prep data frame for graph
  graph_df <- data.frame()
  for (i in base::seq_len(nrow(result_df))) {
    up_genes <- unlist(strsplit(result_df$Up_regulated[i], ", "))
    down_genes <- unlist(strsplit(result_df$Down_regulated[i], ", "))
    for (gene in c(up_genes, down_genes)) {
      graph_df <- rbind(graph_df, data.frame(
        Gene = gene,
        Term = result_df[i, ID_column]
      ))
    }
  }

  ### Merging genes and terms if `genes_df` is supplied
  if (!is.null(genes_df)) {
    graph_df <- merge(
      x = graph_df,
      y = genes_df,
      by.x = "Gene",
      by.y = "Gene.symbol",
      all = TRUE
    )
    graph_df <- graph_df[!is.na(graph_df$Term), ]
  }

  up_genes <- lapply(result_df$Up_regulated, function(x) unlist(strsplit(x, ", ")))
  up_genes <- unlist(up_genes)

  ############ Create graph object and plot create igraph object
  g <- igraph::graph_from_data_frame(graph_df, directed = FALSE)
  cond_term <- names(igraph::V(g)) %in% result_df[, ID_column]

  if (!is.null(genes_df)) {
    cond_gene <- !cond_term

    # store a node class for layout/shape if needed
    #-----------------------------------------------------
    igraph::V(g)$type <- ifelse(cond_term, "term", "gene")

    # store logFC only for gene nodes
    #-----------------------------------------------------
    gene_names <- igraph::V(g)$name[cond_gene]
    gene_logFC <- graph_df$logFC[match(gene_names, graph_df$Gene)]

    # Create a full-length vector with NA for term nodes
    #-----------------------------------------------------
    logFC_full <- rep(NA_real_, igraph::vcount(g))
    suppressWarnings(logFC_full[cond_gene] <- gene_logFC)
    igraph::V(g)$logFC <- logFC_full
  } else {
    up_genes <- lapply(result_df$Up_regulated, function(x) unlist(strsplit(x, ", ")))
    up_genes <- unlist(up_genes)

    cond_up_gene <- names(igraph::V(g)) %in% up_genes

    node_type <- ifelse(cond_term, "term", ifelse(cond_up_gene, "up", "down"))
    node_type <- factor(node_type, levels = c("term", "up", "down"))
    node_type <- droplevels(node_type)
    igraph::V(g)$type <- node_type
  }

  # Adjust node sizes
  if (term_size == "num_genes") {
    sizes <- igraph::degree(g)
    sizes <- ifelse(igraph::V(g)$type == "term", sizes, 2)
  } else {
    idx <- match(names(igraph::V(g)), result_df[, ID_column])
    sizes <- -log10(result_df$lowest_p[idx])
    sizes[is.na(sizes)] <- 2
  }
  igraph::V(g)$size <- sizes
  igraph::V(g)$label.cex <- 0.5
  igraph::V(g)$frame.color <- "gray"

  if (use_edge_weights) {
    gene_term_counts <- table(graph_df$Gene)

    edge_names <- gene_term_counts[igraph::as_data_frame(g, "edges")$from]
    igraph::E(g)$weight <- as.numeric(edge_names)
  }

  if (!is.null(term_fill)) {
    term_rows <- igraph::V(g)$type == "term"
    pathway_names <- igraph::V(g)$name[term_rows]
    matching_row_orders <- match(pathway_names, result_df[[ID_column]])
    node_fill_values <- result_df[matching_row_orders, ][[term_fill]]
    igraph::V(g)$term_fill <- ifelse(term_rows, node_fill_values, NA)
  }

  return(g)
}


#' Create Term-Gene Plot
#'
#' @param graph A \link[igraph]{igraph} returned from \link[pathfindR]{create_term_gene_graph}.
#' @param layout The type of layout to create (see \code{\link[ggraph]{ggraph}} for details (default: \code{'stress'})
#' @param gene_node_fill A character vector to customize the fill gradient colors of the gene nodes when `genes_df` is supplied, color order is in low -> mid -> high (default: \code{c("#7E2795", "white", "#27AE60")}).
#' @param term_node_fill A character vector to customize the fill gradient colors of the term nodes when `term_fill` is supplied, color order is in low -> mid -> high (default: \code{c("#CCBB44", "white", "#4477AA")}).
#' @param gene_node_color A character vector to customize the fill gradient colors of the term nodes when `genes_df` is not supplied, color order is in up -> down (default: \code{c("green", "red")}).
#' @param term_node_color A character to customize the fill color of the terms when `term_fill` is not specified (default: \code{"#E5D7BF"}).
#' @param term_fill_label A character to change the term node legend name (default: \code{NULL}).
#' @param term_size_label A character to change the term node size legend name (default: \code{NULL}).
#' @return A \link[ggraph]{ggraph} object
#'
#' @details
#' This function creates a visualization of the term-gene graph (adapted from the Gene-Concept network visualization in the \code{enrichplot} package).
#' It displays which input genes are involved in enriched biological terms, showing connections between genes and pathway/terms nodes.
#' The graph facilitates investigation of multi-term relationships and identifies shared versus distinct genes across enriched terms.
#'
#' Node coloring depends on the inputs provided to \link[pathfindR]{create_term_gene_graph}:
#' \itemize{
#'   \item If `genes_df` was NOT supplied: term nodes are beige (`term_node_color`), up-regulated genes are green, and down-regulated genes are red.
#'   \item If `genes_df` WAS supplied: gene nodes are colored by logFC using a gradient (`gene_node_fill`: default purple → white → green),
#'         and term nodes can be colored by `term_fill` values (default yellow → white → blue) if `term_fill` was provided.
#' }
#'
#' Term node size reflects either the number of associated genes (`term_size = 'num_genes'`) or statistical significance (`term_size = 'p_val'`).
#' When `use_edge_weights = TRUE` was set in `create_term_gene_graph`, edge widths represent hub gene importance (genes appearing in multiple terms).
#' The layout can be customized via the `layout` parameter (default: "stress"), and legends automatically reflect the applied coloring schemes.
#'
#' @import ggraph
#' @export
#'
#' @examples
#' # Normal gene-term with up/down regulated genes
#' g <- create_term_gene_graph(
#'   result_df = example_pathfindR_output
#' )
#' plt <- create_term_gene_plot(g)
create_term_gene_plot <- function(
  graph,
  layout = "stress",
  gene_node_fill = c("#7E2795", "white", "#27AE60"),
  term_node_fill = c("#CCBB44", "white", "#4477AA"),
  gene_node_color = c("green", "red"),
  term_node_color = "#E5D7BF",
  term_fill_label = NULL,
  term_size_label = NULL
) {
  ## Error handling
  #--------------------------------------------------------------------#

  if (!inherits(graph, "igraph")) {
    stop("`graph` needs to be of class 'igraph'!")
  }

  # Extract essential values
  term_fill <- igraph::V(graph)$term_fill
  num_terms <- sum(igraph::V(graph)$type == "term")
  gene_node_values <- igraph::V(graph)$logFC
  term_node_values <- igraph::V(graph)$term_fill
  weight_node_values <- igraph::E(graph)$weight

  if (length(gene_node_fill) == 3) {
    if (!all(sapply(X = gene_node_fill, FUN = isColor))) {
      stop("Not all elements in `gene_node_fill` are valid colors!")
    }
  } else {
    stop("`gene_node_fill` needs to be of length 3!")
  }

  if (length(term_node_fill) == 3) {
    if (!all(sapply(X = term_node_fill, FUN = isColor))) {
      stop("Not all elements in `term_node_fill` are valid colors!")
    }
  } else {
    stop("`term_node_fill` needs to be of length 3!")
  }

  if (length(gene_node_color) == 2) {
    if (!all(sapply(X = gene_node_color, FUN = isColor))) {
      stop("Not all elements in `gene_node_color` are valid colors!")
    }
  } else {
    stop("`gene_node_color` needs to be of length 2!")
  }

  if (!isColor(term_node_color)) {
    stop("`term_node_color` is not a valid color!")
  }

  ## MAIN
  #--------------------------------------------------------------------#

  if (is.null(gene_node_values) || is.null(term_fill)) {
    node_type <- igraph::V(graph)$type

    type_descriptions <- c(term = "enriched term", up = "up-regulated gene", down = "down-regulated gene")
    type_descriptions <- type_descriptions[levels(node_type)]

    node_colors <- c(term_node_color, gene_node_color)

    names(node_colors) <- names(type_descriptions)
    node_colors <- node_colors[levels(node_type)]
  }

  ### Create graph
  if (!is.null(weight_node_values)) {
    p <- ggraph::ggraph(graph, layout = layout) +
      ggraph::geom_edge_link(
        mapping = ggplot2::aes(
          width = .data$weight * 0.1
        ),
        color = "darkgrey",
        alpha = 0.8,
        show.legend = FALSE
      )
  } else {
    p <- ggraph::ggraph(graph, layout = layout) +
      ggraph::geom_edge_link(
        alpha = 0.8,
        color = "darkgrey",
        show.legend = FALSE
      )
  }

  # First layer for gene nodes, if `genes_df` is supplied
  if (!is.null(gene_node_values)) {
    p1 <- p +
      ggraph::scale_edge_width(guide = "none") +
      ggraph::geom_node_point(
        mapping = ggplot2::aes(
          fill = .data$logFC,
          size = .data$size
        ),
        shape = 21,
        colour = "black",
        show.legend = TRUE
      ) +
      ggplot2::scale_fill_gradient2(
        low = gene_node_fill[1],
        mid = gene_node_fill[2],
        high = gene_node_fill[3],
        name = "LogFC"
      )
  } else {
    # `genes_df` is absent, coloring up/down and term nodes
    p1 <- p +
      ggraph::scale_edge_width(guide = "none") +
      ggraph::geom_node_point(
        mapping = ggplot2::aes(
          fill = .data$type,
          size = .data$size
        ),
        shape = 21,
        colour = "black"
      ) +
      ggplot2::scale_fill_manual(
        values = node_colors,
        labels = type_descriptions
      )
  }


  if (!is.null(gene_node_values) || !is.null(term_fill)) {
    ## Extract plot data
    plot_data <- ggplot2::ggplot_build(p1)$data

    ## Extract graph data
    term_rows <- igraph::V(graph)$type == "term"

    ## Second layer contains the previous `geom_node_point`
    node_data <- plot_data[[2]]
    node_data$type <- igraph::V(graph)$type
    term_data <- node_data[term_rows, ]
    term_data$term_fill <- stats::na.omit(igraph::V(graph)$term_fill)
    term_data$size <- igraph::V(graph)$size[term_rows]

    gene_term <- node_data[!term_rows, ]
    gene_term$size <- igraph::V(graph)$size[!term_rows]

    if (!is.null(gene_node_values)) {
      gene_term$logFC <- stats::na.omit(igraph::V(graph)$logFC)
      p1 <- p +
        # First gene layer
        ggraph::scale_edge_width(guide = "none") +
        ggraph::geom_node_point(
          data = gene_term,
          mapping = ggplot2::aes(
            x = .data$x,
            y = .data$y,
            fill = .data$logFC,
            size = .data$size
          ),
          shape = 21,
          colour = "black",
          show.legend = TRUE
        ) +
        ggplot2::scale_fill_gradient2(
          low = gene_node_fill[1],
          mid = gene_node_fill[2],
          high = gene_node_fill[3],
          name = "LogFC"
        )

      if (!is.null(term_fill)) {
        p1 <- p1 +
          # Second Term layer
          ggnewscale::new_scale_fill() +
          ggraph::geom_node_point(
            data = term_data,
            mapping = ggplot2::aes(
              x = .data$x,
              y = .data$y,
              fill = .data$term_fill,
              size = .data$size
            ),
            shape = 21,
            colour = "black",
            show.legend = TRUE
          ) +
          ggplot2::scale_fill_gradient2(
            low = term_node_fill[1],
            mid = term_node_fill[2],
            high = term_node_fill[3],
            name = ifelse(is.null(term_fill_label), "term_fill", term_fill_label)
          )
      } else {
        p1 <- p1 +
          # Second Term layer
          ggnewscale::new_scale_fill() +
          ggraph::geom_node_point(
            data = term_data,
            mapping = ggplot2::aes(
              x = .data$x,
              y = .data$y,
              fill = .data$type,
              size = .data$size
            ),
            shape = 21,
            colour = "black",
            show.legend = TRUE
          ) +
          ggplot2::scale_fill_manual(
            values = node_colors[1],
            labels = type_descriptions[1]
          )
      }
    } else {
      p1 <- p +
        ggraph::scale_edge_width(guide = "none") +
        ggraph::geom_node_point(
          data = gene_term,
          mapping = ggplot2::aes(
            fill = .data$type,
            size = .data$size
          ),
          colour = "black",
          shape = 21
        ) +
        ggplot2::scale_fill_manual(
          values = node_colors[2:3],
          labels = type_descriptions[2:3]
        )

      if (!is.null(term_fill)) {
        p1 <- p1 +
          # Second Term layer
          ggnewscale::new_scale_fill() +
          ggraph::geom_node_point(
            data = term_data,
            mapping = ggplot2::aes(
              x = .data$x,
              y = .data$y,
              fill = .data$term_fill,
              size = .data$size
            ),
            shape = 21,
            colour = "black",
            show.legend = TRUE
          ) +
          ggplot2::scale_fill_gradient2(
            low = term_node_fill[1],
            mid = term_node_fill[2],
            high = term_node_fill[3],
            name = paste0(term_fill)
          )
      }
    }
  }

  p <- p1 +
    ggplot2::scale_size(
      range = c(5, 10),
      breaks = round(seq(
        round(min(igraph::V(graph)$size)),
        round(max(igraph::V(graph)$size)),
        length.out = 4
      )),
      name = ifelse(is.null(term_size_label), "term size", term_size_label)
    ) +
    ggplot2::theme_void() +
    suppressWarnings(
      ggraph::geom_node_text(
        mapping = ggplot2::aes(label = .data$name),
        nudge_y = 0.2,
        repel = TRUE,
        max.overlaps = 20
      )
    )
  if (!is.null(num_terms)) {
    p <- p + ggplot2::ggtitle(
      "Term-Gene Graph",
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
