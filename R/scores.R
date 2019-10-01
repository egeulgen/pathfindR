#' Calculate Agglomerated Scores of Enriched Terms for Each Subject
#'
#' @param enrichment_table a data frame that must contain the 3 columns below: \describe{
#'   \item{Term_Description}{Description of the enriched term (necessary if \code{use_description = TRUE})}
#'   \item{ID}{ID of the enriched term (necessary if \code{use_description = FALSE})}
#'   \item{Up_regulated}{the up-regulated genes in the input involved in the given term's gene set, comma-separated}
#'   \item{Down_regulated}{the down-regulated genes in the input involved in the given term's gene set, comma-separated}
#' }
#' @param exp_mat the experiment (e.g., gene expression/methylation) matrix.
#' Columns are samples and rows are genes. Column names must contain sample
#' names and row names must contain the gene symbols.
#' @param cases (Optional) A vector of sample names that are cases in the
#' case/control experiment. (default = NULL)
#' @param use_description Boolean argument to indicate whether term descriptions
#'  (in the "Term_Description" column) should be used. (default = \code{FALSE})
#' @param plot_hmap Boolean value to indicate whether or not to draw the
#' heatmap plot of the scores. (default = TRUE)
#' @param ... Additional arguments for \code{\link{plot_scores}} for aesthetics
#' of the heatmap plot
#'
#' @return Matrix of agglomerated scores of each enriched term per sample.
#' Columns are samples, rows are enriched terms. Optionally, displays a heatmap
#' of this matrix.
#'
#' @export
#'
#' @examples
#' score_matrix <- calculate_scores(RA_output, RA_exp_mat, plot_hmap = FALSE)
calculate_scores <- function(enrichment_table, exp_mat, cases = NULL,
                             use_description = FALSE,
                             plot_hmap = TRUE,
                             ...) {

  ############ Argument Checks
  ### Check use_description is boolean
  if (!is.logical(use_description)) {
    stop("`use_description` must either be TRUE or FALSE!")
  }

  ### Check whether any cases (if provided) are missing from exp_mat
  if (!is.null(cases) & any(!cases %in% colnames(exp_mat))) {
    stop("Missing cases in the expression matrix!")
  }

  ### Set column for term labels
  ID_column <- ifelse(use_description, "Term_Description", "ID")

  all_scores_matrix <- c()
  for (i in base::seq_len(nrow(enrichment_table))) {
    # Get DEGs
    up_genes <- enrichment_table$Up_regulated[i]
    down_genes <- enrichment_table$Down_regulated[i]
    up_genes <- unlist(strsplit(up_genes, ", "))
    down_genes <- unlist(strsplit(down_genes, ", "))

    genes <- c(up_genes, down_genes)

    # some genes may not be in exp. matrix
    genes <- genes[genes %in% rownames(exp_mat)]

    if (length(genes) != 0) {
      # subset exp. matrix to include only DEGs
      sub_mat <- exp_mat[rownames(exp_mat) %in% genes, , drop = FALSE]

      current_term_score_matrix <- c()
      for (gene in genes) {
        gene_vec <- sub_mat[rownames(sub_mat) == gene, ]
        gene_vec <- as.numeric(gene_vec)
        names(gene_vec) <- colnames(sub_mat)

        # calculate mean and sd across samples
        gene_mean <- base::mean(gene_vec)
        gene_sd <- stats::sd(gene_vec)

        gene_scores <- vapply(
          gene_vec, function(x) (x - gene_mean) / gene_sd,
          1.2
        )
        current_term_score_matrix <- rbind(current_term_score_matrix, gene_scores)
        rownames(current_term_score_matrix)[nrow(current_term_score_matrix)] <- gene
      }

      current_term_scores <- apply(current_term_score_matrix, 2, base::mean)
      all_scores_matrix <- rbind(all_scores_matrix, current_term_scores)
      rownames(all_scores_matrix)[nrow(all_scores_matrix)] <- enrichment_table[i, ID_column]
    }
  }

  if (!is.null(cases)) {
    ## order as cases, then controls
    match1 <- match(cases, colnames(all_scores_matrix))
    match2 <- setdiff(base::seq_len(ncol(all_scores_matrix)), match1)
    all_scores_matrix <- all_scores_matrix[, c(match1, match2)]
  }

  if (plot_hmap) {
    heatmap <- plot_scores(score_matrix = all_scores_matrix, cases = cases, ...)
    graphics::plot(heatmap)
  }

  return(all_scores_matrix)
}

#' Plot the Heatmap of Score Matrix of Enriched Terms per Sample
#'
#' @param score_matrix Matrix of agglomerated enriched term scores per sample. Columns are
#' samples, rows are enriched terms
#' @param cases (Optional) A vector of sample names that are cases in the
#' case/control experiment.
#' @param label_samples Boolean value to indicate whether or not to label the
#' samples in the heatmap plot (default = TRUE)
#' @param case_control_titles A vector of length two for naming of the 'Case'
#' and 'Control' groups (in order) (default = \code{c('Case', 'Control')})
#' @param low a string indicating the color of 'low' values in the score coloring gradient (default = 'green')
#' @param high a string indicating the color of 'high' values in the score coloring gradient (default = 'red')
#'
#' @return A `ggplot2` object containing the heatmap plot. x-axis indicates
#' the samples. y-axis indicates the enriched terms. "Score" indicates the
#' score of the term in a given sample. If \code{cases} are provided, the plot is
#' divided into 2 facets, named by the \code{case_control_titles}.
#'
#' @import ggplot2
#' @export
#'
#' @examples
#' score_mat <- calculate_scores(RA_output, RA_exp_mat, plot_hmap = FALSE)
#' hmap <- plot_scores(score_mat)
plot_scores <- function(score_matrix, cases = NULL, label_samples = TRUE,
                        case_control_titles = c("Case", "Control"),
                        low = "green", high = "red") {

  ############ Argument Checks
  if (length(case_control_titles) != 2) {
    stop("\"case_control_titles\" must contain exactly two elements!")
  }
  if (any(!cases %in% colnames(score_matrix)) & !is.null(cases)) {
    stop("Missing cases in the score matrix!")
  }

  ## sort according to activity (up/down)
  if (!is.null(cases)) {
    tmp <- rowMeans(score_matrix[, cases, drop = FALSE])
    score_matrix <- score_matrix[c(which(tmp >= 0), which(tmp < 0)), ]
  }

  ## transform the matrix
  var_names <- list()
  var_names[["Term"]] <- factor(rownames(score_matrix),
    levels = rev(rownames(score_matrix))
  )
  var_names[["Sample"]] <- factor(colnames(score_matrix),
    levels = colnames(score_matrix)
  )

  score_df <- expand.grid(var_names,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  scores <- as.vector(score_matrix)
  scores <- data.frame(scores)
  score_df <- cbind(score_df, scores)
  if (!is.null(cases)) {
    score_df$Type <- factor(ifelse(score_df$Sample %in% cases,
      case_control_titles[1],
      case_control_titles[2]
    ),
    levels = case_control_titles
    )
  }

  g <- ggplot2::ggplot(score_df, ggplot2::aes_(x = ~Sample, y = ~Term))
  g <- g + ggplot2::geom_tile(ggplot2::aes_(fill = ~scores), color = "white")
  g <- g + ggplot2::scale_fill_gradient2(low = low, mid = "black", high = high)
  g <- g + ggplot2::theme(
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1
    ),
    legend.title = ggplot2::element_text(size = 10),
    legend.text = ggplot2::element_text(size = 12)
  )
  g <- g + ggplot2::labs(fill = "Score")
  if (!is.null(cases)) {
    g <- g + ggplot2::facet_grid(~Type, scales = "free_x", space = "free")
    g <- g + ggplot2::theme(strip.text.x = ggplot2::element_text(
      size = 12,
      face = "bold"
    ))
  }
  if (!label_samples) {
    g <- g + ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )
  }
  return(g)
}
