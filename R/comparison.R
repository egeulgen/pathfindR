#' Combine 2 pathfindR Results
#'
#' @param result_A data frame of first pathfindR enrichment results
#' @param result_B data frame of second pathfindR enrichment results
#' @param plot_common boolean to indicate whether or not to plot the term-gene
#' graph of the common terms (default=\code{TRUE})
#'
#' @return Data frame of combined pathfindR enrichment results. Columns are: \describe{
#'   \item{ID}{ID of the enriched term}
#'   \item{Term_Description}{Description of the enriched term}
#'   \item{Fold_Enrichment_A}{Fold enrichment value for the enriched term (Calculated using ONLY the input genes)}
#'   \item{occurrence_A}{the number of iterations that the given term was found to enriched over all iterations}
#'   \item{lowest_p_A}{the lowest adjusted-p value of the given term over all iterations}
#'   \item{highest_p_A}{the highest adjusted-p value of the given term over all iterations}
#'   \item{Up_regulated_A}{the up-regulated genes in the input involved in the given term's gene set, comma-separated}
#'   \item{Down_regulated_A}{the down-regulated genes in the input involved in the given term's gene set, comma-separated}
#'   \item{Fold_Enrichment_B}{Fold enrichment value for the enriched term (Calculated using ONLY the input genes)}
#'   \item{occurrence_B}{the number of iterations that the given term was found to enriched over all iterations}
#'   \item{lowest_p_B}{the lowest adjusted-p value of the given term over all iterations}
#'   \item{highest_p_B}{the highest adjusted-p value of the given term over all iterations}
#'   \item{Up_regulated_B}{the up-regulated genes in the input involved in the given term's gene set, comma-separated}
#'   \item{Down_regulated_B}{the down-regulated genes in the input involved in the given term's gene set, comma-separated}
#'   \item{combined_p}{the combined p value (via Fisher's method)}
#'   \item{status}{whether the term is found in both analyses ("common"), found only in the first ("A only") or found only in the second ("B only)}
#' }
#' By default, the function also displays the term-gene graph of the common terms
#'
#' @export
#'
#' @examples
#' combined_results <- combine_pathfindR_results(example_pathfindR_output, example_comparison_output)
combine_pathfindR_results <- function(result_A, result_B, plot_common = TRUE) {
  combined_df <- merge(
    result_A, result_B,
    by = c("ID", "Term_Description"), all = TRUE,
    suffixes = c("_A", "_B")
  )

  ### Calculate combined p values
  combined_df$combined_p <- NA
  for (i in seq_len(nrow(combined_df))) {
    p_vec <- c(combined_df$lowest_p_A[i], combined_df$lowest_p_B[i])
    p_vec <- p_vec[!is.na(p_vec)]
    combined_df$combined_p[i] <- stats::pchisq(
      q = sum(log(p_vec)) * -2,
      df = length(p_vec) * 2,
      lower.tail = FALSE
    )
  }
  ### Indicate intersection status
  combined_df$status <- ifelse(is.na(combined_df$lowest_p_A), "B only",
    ifelse(is.na(combined_df$lowest_p_B), "A only", "common")
  )

  ### Plot graph common terms
  if (plot_common) {
    graphics::plot(combined_results_graph(combined_df))
  }

  message("You may run `combined_results_graph()` to create visualizations of combined term-gene graphs of selected terms")

  return(combined_df)
}



#' Combined Results Graph
#'
#' @param combined_df Data frame of combined pathfindR enrichment results
#' @param selected_terms the vector of selected terms for creating the graph
#' (either IDs or term descriptions). If set to \code{"common"}, all of the
#' common terms are used. (default = "common")
#' @inheritParams term_gene_graph
#'
#' @return a  \code{\link[ggraph]{ggraph}} object containing the combined term-gene graph.
#'  Each node corresponds to an enriched term (orange if common, different shades of blue otherwise),
#'  an up-regulated gene (green), a down-regulated gene (red) or
#'  a conflicting (i.e. up in one analysis, down in the other or vice versa) gene
#'  (gray). An edge between a term and a gene indicates
#'  that the given term involves the gene. Size of a term node is proportional
#'  to either the number of genes (if \code{node_size = "num_genes"}) or
#'  the -log10(lowest p value) (if \code{node_size = "p_val"}).
#' @export
#'
#' @examples
#' combined_results <- combine_pathfindR_results(
#'   example_pathfindR_output,
#'   example_comparison_output,
#'   plot_common = FALSE
#' )
#' g <- combined_results_graph(combined_results, selected_terms = sample(combined_results$ID, 3))
combined_results_graph <- function(combined_df, selected_terms = "common",
                                   use_description = FALSE,
                                   layout = "stress",
                                   node_size = "num_genes") {
  ############ Argument Checks
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

  if (!is.data.frame(combined_df)) {
    stop("`combined_df` should be a data frame")
  }

  ### Check necessary columnns
  necessary_cols <- c(
    ID_column, "combined_p",
    "Up_regulated_A", "Down_regulated_A",
    "Up_regulated_B", "Down_regulated_B"
  )

  if (!all(necessary_cols %in% colnames(combined_df))) {
    stop(paste(c(
      "All of", paste(necessary_cols, collapse = ", "),
      "must be present in `results_df`!"
    ), collapse = " "))
  }

  ############ Initial steps
  ### Filter for selected terms
  if (any(selected_terms == "common")) {
    if (!any(combined_df$status == "common")) {
      stop("There are no common terms")
    }
    combined_df <- combined_df[combined_df$status == "common", ]
  } else {
    if (!any(selected_terms %in% combined_df[, ID_column])) {
      stop("None of the `selected_terms` are in the combined results!")
    }
    combined_df <- combined_df[combined_df[, ID_column] %in% selected_terms, ]
  }

  ### Prep data frame for graph
  graph_df <- data.frame()
  for (i in base::seq_len(nrow(combined_df))) {
    up_genes <- c(
      unlist(strsplit(combined_df$Up_regulated_A[i], ", ")),
      unlist(strsplit(combined_df$Up_regulated_B[i], ", "))
    )
    down_genes <- c(
      unlist(strsplit(combined_df$Down_regulated_A[i], ", ")),
      unlist(strsplit(combined_df$Down_regulated_B[i], ", "))
    )
    genes <- c(up_genes, down_genes)
    genes <- genes[!is.na(genes)]
    for (gene in genes) {
      graph_df <- rbind(
        graph_df,
        data.frame(
          Term = combined_df[i, ID_column],
          Gene = gene,
          stringsAsFactors = FALSE
        )
      )
    }
  }
  graph_df <- unique(graph_df)

  up_genes_A <- unlist(lapply(
    combined_df$Up_regulated_A,
    function(x) unlist(strsplit(x, ", "))
  ))
  down_genes_A <- unlist(lapply(
    combined_df$Down_regulated_A,
    function(x) unlist(strsplit(x, ", "))
  ))
  up_genes_B <- unlist(lapply(
    combined_df$Up_regulated_B,
    function(x) unlist(strsplit(x, ", "))
  ))
  down_genes_B <- unlist(lapply(
    combined_df$Down_regulated_B,
    function(x) unlist(strsplit(x, ", "))
  ))

  terms_A <- combined_df[!is.na(combined_df$lowest_p_A) & is.na(combined_df$lowest_p_B), ID_column]
  terms_B <- combined_df[is.na(combined_df$lowest_p_A) & !is.na(combined_df$lowest_p_B), ID_column]

  ############ Create graph object and plot
  ### create igraph object
  g <- igraph::graph_from_data_frame(graph_df, directed = FALSE)
  igraph::V(g)$type <- ifelse(names(igraph::V(g)) %in% terms_A, "A-only term",
    ifelse(names(igraph::V(g)) %in% terms_B, "B-only term",
      ifelse(names(igraph::V(g)) %in% combined_df[, ID_column], "common term", "gene")
    )
  )

  # Adjust node sizes
  if (node_size == "num_genes") {
    sizes <- igraph::degree(g)
    sizes <- ifelse(grepl("term", igraph::V(g)$type), sizes, 2)
    size_label <- "# genes"
  } else {
    idx <- match(names(igraph::V(g)), combined_df[, ID_column])
    sizes <- -log10(combined_df$combined_p[idx])
    sizes[is.na(sizes)] <- 2
    size_label <- "-log10(p)"
  }
  igraph::V(g)$size <- sizes
  igraph::V(g)$label.cex <- 0.5
  igraph::V(g)$frame.color <- "gray"

  cond_up_A <- names(igraph::V(g)) %in% up_genes_A
  cond_up_B <- names(igraph::V(g)) %in% up_genes_B
  cond_down_A <- names(igraph::V(g)) %in% down_genes_A
  cond_down_B <- names(igraph::V(g)) %in% down_genes_B
  missing_A <- !cond_up_A & !cond_down_A
  missing_B <- !cond_up_B & !cond_down_B

  up_cond <- (cond_up_A & cond_up_B) | (missing_A & cond_up_B) | (cond_up_A & missing_B)
  down_cond <- (cond_down_A & cond_down_B) | (missing_A & cond_down_B) | (cond_down_A & missing_B)

  igraph::V(g)$for_coloring <- ifelse(igraph::V(g)$type == "common term", "Common term",
    ifelse(igraph::V(g)$type == "A-only term", "A-only term",
      ifelse(igraph::V(g)$type == "B-only term", "B-only term",
        ifelse(up_cond, "Up gene",
          ifelse(down_cond, "Down gene", "Conflicting gene")
        )
      )
    )
  )

  ### Create graph
  p <- ggraph::ggraph(g, layout = layout)
  p <- p + ggraph::geom_edge_link(alpha = .8, colour = "darkgrey")
  p <- p + ggraph::geom_node_point(ggplot2::aes_(color = ~for_coloring, size = ~size))
  p <- p + ggplot2::scale_size(
    range = c(5, 10),
    breaks = round(seq(round(min(igraph::V(g)$size)),
      round(max(igraph::V(g)$size)),
      length.out = 4
    )),
    name = size_label
  )
  p <- p + ggplot2::theme_void()
  p <- p + suppressWarnings(ggraph::geom_node_text(ggplot2::aes_(label = ~name),
    nudge_y = .2,
    repel = TRUE, max.overlaps = 20
  ))

  vertex_cols <- c(
    "Common term" = "#FCCA46",
    "A-only term" = "#9FB8AD",
    "B-only term" = "#619B8A",
    "Up gene" = "green",
    "Down gene" = "red",
    "Conflicting gene" = "gray"
  )
  p <- p + ggplot2::scale_colour_manual(
    values = vertex_cols,
    name = NULL
  )
  p <- p + ggplot2::ggtitle("Combined Terms Graph")
  p <- p + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  return(p)
}
