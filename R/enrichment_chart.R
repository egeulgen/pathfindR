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
enrichment_chart <- function(result_df, plot_by_cluster = FALSE) {
  necessary <- c("Pathway", "Fold_Enrichment", "lowest_p", "Up_regulated", "Down_regulated")
  if (!all(necessary %in% colnames(result_df)))
    stop("The input data frame must have the columns: Pathway, Fold_Enrichment, lowest_p, Up_regulated, Down_regulated")

  if (!is.logical(plot_by_cluster))
    stop("plot_by_cluster must be either TRUE or FALSE")

  # sort by lowest adj.p
  result_df <- result_df[order(result_df$lowest_p), ]

  n <- sapply(result_df$Up_regulated, function(x) length(unlist(strsplit(x, ", "))))
  n <- n + sapply(result_df$Down_regulated, function(x) length(unlist(strsplit(x, ", "))))

  result_df$Pathway <- factor(result_df$Pathway, levels = rev(result_df$Pathway))

  g <- ggplot2::ggplot(result_df, ggplot2::aes_(x = ~Fold_Enrichment, y = ~Pathway))
  g <- g + ggplot2::geom_point(ggplot2::aes(color = -log10(result_df$lowest_p),
                                            size = n), na.rm = TRUE)
  g <- g + ggplot2::theme_bw()
  g <- g + ggplot2::theme(axis.text.x = ggplot2::element_text(size = 10),
                          axis.text.y = ggplot2::element_text(size = 10),
                          plot.title = ggplot2::element_blank())
  g <- g + ggplot2::xlab("Fold Enrichment") + ggplot2::ylab("")
  g <- g + ggplot2::labs(size = "# of DEGs", color = "-log10(lowest-p)")
  g <- g + ggplot2::scale_color_continuous(low = "#f5efef", high = "red")

  if (plot_by_cluster & "Cluster" %in% colnames(result_df)) {
    g <- g + ggplot2::facet_grid(result_df$Cluster~., scales = "free_y", space = "free", drop = TRUE)
  } else if (plot_by_cluster) {
    warning("For plotting by cluster, there must a column named `Cluster` in the input data frame!")
  }

  return(g)
}
