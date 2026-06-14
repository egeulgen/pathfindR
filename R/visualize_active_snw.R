#' Visualize Active Subnetworks
#'
#' @inheritParams filter_active_subnetworks
#' @inheritParams term_gene_heatmap
#' @inheritParams return_pin_path
#' @param num_snws number of top subnetworks to be visualized (leave blank if
#' you want to visualize all subnetworks)
#' @param layout The type of layout to create (see \code{\link[ggraph]{ggraph}} for details (default: \code{'stress'})
#' @param ... additional arguments for \code{\link{input_processing}}
#'
#' @return a list of ggplot objects of graph visualizations of identified active
#' subnetworks. Green nodes are down-regulated genes, reds are up-regulated genes
#' and yellows are non-input genes
#' @export
#'
#' @examples
#' # visualize top 2 active subnetworks
#' g_list <- visualize_active_subnetworks(
#'   active_snws = example_unfiltered_snws,
#'   genes_df = example_pathfindR_input[1:10, ],
#'   pin_name_path = "KEGG",
#'   num_snws = 2
#' )
visualize_active_subnetworks <- function(active_snws, genes_df, pin_name_path = "Biogrid",
                                         num_snws, layout = "stress", score_quan_thr = 0.8, sig_gene_thr = 0.02, ...) {
  # process input data frame
  processed_input <- input_processing(genes_df,
    pin_name_path = pin_name_path,
    ...
  )

  # parse and filter active subnetworks
  active_snw_list <- filter_active_subnetworks(
    active_snws = active_snws, sig_genes_vec = processed_input$GENE,
    score_quan_thr = score_quan_thr, sig_gene_thr = sig_gene_thr
  )
  if (is.null(active_snw_list) | length(active_snw_list$scores) == 0) {
    return(NULL)
  }

  score_vec <- active_snw_list$scores
  subnetworks <- active_snw_list$subnetworks

  if (missing(num_snws)) {
    num_snws <- length(subnetworks)
  }

  if (num_snws > length(subnetworks)) {
    num_snws <- length(subnetworks)
  }

  # load PIN data load PIN
  pin_path <- return_pin_path(pin_name_path)
  pin <- utils::read.delim(file = pin_path, header = FALSE)
  pin$V2 <- NULL

  pin[, 1] <- base::toupper(pin[, 1])
  pin[, 2] <- base::toupper(pin[, 2])

  # create graphs
  graphs_list <- list()
  for (idx in seq_len(num_snws)) {
    snw <- subnetworks[[idx]]

    num_input_genes <- sum(processed_input$GENE %in% snw)
    perc_input_genes <- round(num_input_genes / length(processed_input$GENE) *
      100, 2)

    snw_interactions <- pin[pin[, 1] %in% snw & pin[, 2] %in% snw, ]
    g <- igraph::graph_from_data_frame(snw_interactions, directed = FALSE)
    cond_up_gene <- names(igraph::V(g)) %in% processed_input$GENE[processed_input$CHANGE >
      0]
    cond_down_gene <- names(igraph::V(g)) %in% processed_input$GENE[processed_input$CHANGE <
      0]
    igraph::V(g)$type <- ifelse(cond_up_gene, "up", ifelse(cond_down_gene, "down",
      "non-input"
    ))

    igraph::V(g)$label.cex <- 0.5
    igraph::V(g)$frame.color <- "gray"
    igraph::V(g)$color <- ifelse(igraph::V(g)$type == "non-input", "#FFD500",
      ifelse(igraph::V(g)$type == "up", "#D2222D", "#35CD35")
    )

    color_lookup <- c(
      `#35CD35` = "down-regulated gene", `#D2222D` = "up-regulated gene",
      `#FFD500` = "non-input gene"
    )


    p <- ggraph::ggraph(g, layout = layout)
    p <- p + ggraph::geom_edge_link(alpha = 0.8, colour = "darkgrey")
    p <- p + ggraph::geom_node_point(ggplot2::aes(color = .data$color), size = 2)
    p <- p + ggplot2::theme_void()
    p <- p + ggraph::geom_node_text(ggplot2::aes(label = .data$name), nudge_y = 0.2)
    p <- p + ggplot2::scale_colour_manual(
      values = unique(igraph::V(g)$color),
      name = NULL, labels = color_lookup[unique(igraph::V(g)$color)]
    )
    p <- p + ggplot2::labs(title = paste0("Active Subnetwork #", idx), subtitle = paste0(
      "Score=",
      round(score_vec[idx], 2), ", ", num_input_genes, "(", perc_input_genes,
      "%) input genes"
    ))
    p <- p + ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5), legend.position = "bottom"
    )
    graphs_list[[idx]] <- p
  }

  return(graphs_list)
}
