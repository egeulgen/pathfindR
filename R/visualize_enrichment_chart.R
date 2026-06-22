#' Create Bubble Chart of Enrichment Results
#'
#' This function is used to create a ggplot2 bubble chart displaying the
#' enrichment results.
#'
#' @param result_df a data frame that must contain the following columns: \describe{
#'   \item{Term_Description}{Description of the enriched term}
#'   \item{Fold_Enrichment}{Fold enrichment value for the enriched term}
#'   \item{lowest_p}{the lowest adjusted-p value of the given term over all iterations}
#'   \item{Up_regulated}{the up-regulated genes in the input involved in the given term's gene set, comma-separated}
#'   \item{Down_regulated}{the down-regulated genes in the input involved in the given term's gene set, comma-separated}
#'   \item{Cluster(OPTIONAL)}{the cluster to which the enriched term is assigned}
#' }
#' @param top_terms number of top terms (according to the 'lowest_p' column)
#'  to plot (default = 10). If \code{plot_by_cluster = TRUE}, selects the top
#'  \code{top_terms} terms per each cluster. Set \code{top_terms = NULL} to plot
#'  for all terms.If the total number of terms is less than \code{top_terms},
#'  all terms are plotted.
#' @param plot_by_cluster boolean value indicating whether or not to group the
#'  enriched terms by cluster (works if \code{result_df} contains a
#'  'Cluster' column).
#' @param num_bubbles number of sizes displayed in the legend \code{# genes}
#'  (Default = 4)
#' @param even_breaks whether or not to set even breaks for the number of sizes
#'  displayed in the legend \code{# genes}. If \code{TRUE} (default), sets
#'  equal breaks and the number of displayed bubbles may be different than the
#'  number set by \code{num_bubbles}. If the exact number set by
#'  \code{num_bubbles} is required, set this argument to \code{FALSE}
#' @param order_by the order and coloring of the dots (default: \code{'lowest_p'}).
#'
#' @return a \code{\link[ggplot2]{ggplot2}} object containing the bubble chart.
#' The x-axis corresponds to fold enrichment values while the y-axis indicates
#' the enriched terms. Size of the bubble indicates the number of significant
#' genes in the given enriched term. Color indicates the -log10(lowest-p) value.
#' The closer the color is to red, the more significant the enrichment is.
#' Optionally, if 'Cluster' is a column of \code{result_df} and
#' \code{plot_by_cluster == TRUE}, the enriched terms are grouped by clusters.
#'
#' @import ggplot2
#' @export
#'
#' @examples
#' g <- enrichment_chart(example_pathfindR_output)
enrichment_chart <- function(result_df, top_terms = 10, plot_by_cluster = FALSE,
                             num_bubbles = 4, even_breaks = TRUE, order_by = "lowest_p") {
  message("Plotting the enrichment bubble chart")
  necessary <- c(
    "Term_Description", "Fold_Enrichment", "lowest_p", "Up_regulated",
    "Down_regulated"
  )

  if (!all(necessary %in% colnames(result_df))) {
    stop("The input data frame must have the columns:\n", paste(necessary, collapse = ", "))
  }

  if (!is.logical(plot_by_cluster)) {
    stop("`plot_by_cluster` must be either TRUE or FALSE")
  }

  if (!is.numeric(top_terms) & !is.null(top_terms)) {
    stop("`top_terms` must be either numeric or NULL")
  }

  if (!is.null(top_terms)) {
    if (top_terms < 1) {
      stop("`top_terms` must be > 1")
    }
  }

  result_df <- order_df_by_columnn(result_df, order_by)

  ## Filter for top_terms
  if (!is.null(top_terms)) {
    if (plot_by_cluster & "Cluster" %in% colnames(result_df)) {
      keep_ids <- tapply(result_df$ID, result_df$Cluster, function(x) {
        x[seq_len(min(top_terms, length(x)))]
      })
      keep_ids <- unlist(keep_ids)
      result_df <- result_df[result_df$ID %in% keep_ids, ]
    } else if (top_terms < nrow(result_df)) {
      result_df <- result_df[seq_len(top_terms), ]
    }
  }

  num_genes <- vapply(result_df$Up_regulated, function(x) {
    length(unlist(strsplit(
      x,
      ", "
    )))
  }, 1)
  num_genes <- num_genes + vapply(result_df$Down_regulated, function(x) {
    length(unlist(strsplit(
      x,
      ", "
    )))
  }, 1)

  result_df$Term_Description <- factor(result_df$Term_Description, levels = rev(unique(result_df$Term_Description)))

  g <- ggplot2::ggplot(
    data = result_df,
    mapping = ggplot2::aes(
      x = .data$Fold_Enrichment,
      y = .data$Term_Description
    )
  )

  if (order_by %in% c("lowest_p", "highest_p")) {
    log_p <- -log10(result_df[[order_by]])

    g <- g + ggplot2::geom_point(
      mapping = ggplot2::aes(
        color = log_p,
        size = num_genes
      ),
      na.rm = TRUE
    )

    color_label <- expression(-log[10](p))
  } else {
    g <- g + ggplot2::geom_point(
      mapping = ggplot2::aes(
        color = !!sym(order_by),
        size = num_genes
      ),
      na.rm = TRUE
    )
    color_label <- order_by
  }

  g <- g + ggplot2::theme_bw()
  g <- g + ggplot2::theme(
    axis.text.x = ggplot2::element_text(size = 10), axis.text.y = ggplot2::element_text(size = 10),
    plot.title = ggplot2::element_blank()
  )
  g <- g + ggplot2::xlab("Fold Enrichment")
  g <- g + ggplot2::theme(axis.title.y = ggplot2::element_blank())
  g <- g + ggplot2::labs(size = "# genes", color = color_label)

  ## breaks for # genes
  if (max(num_genes) < num_bubbles) {
    g <- g + ggplot2::scale_size_continuous(breaks = seq(0, max(num_genes)))
  } else {
    if (even_breaks) {
      brks <- base::seq(0, max(num_genes), round(max(num_genes) / (num_bubbles +
        1)))
    } else {
      brks <- base::round(base::seq(0, max(num_genes), length.out = num_bubbles +
        1))
    }
    g <- g + ggplot2::scale_size_continuous(breaks = brks)
  }

  g <- g + ggplot2::scale_color_gradient(low = "#f5efef", high = "red")

  if (plot_by_cluster & "Cluster" %in% colnames(result_df)) {
    g <- g + ggplot2::facet_grid(result_df$Cluster ~ .,
      scales = "free_y", space = "free",
      drop = TRUE
    )
  } else if (plot_by_cluster) {
    message("For plotting by cluster, there must a column named `Cluster` in the input data frame!")
  }

  return(g)
}
