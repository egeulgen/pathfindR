#' Calculate Pathway Scores for Each Subject
#'
#' @param pw_table a data frame that must contain the 3 columns below:\describe{
#'   \item{Pathway}{Description of the enriched pathway}
#'   \item{Up_regulated}{the up-regulated genes in the input involved in the given pathway, comma-separated}
#'   \item{Down_regulated}{the down-regulated genes in the input involved in the given pathway, comma-separated}
#' }
#' @param exp_mat the gene expression/methylation matrix. Columns are samples
#' and rows are genes. Column names must contain sample names and row names must
#' contain the gene symbols.
#' @param cases (Optional) A vector of sample names that are cases in the
#' case/control experiment.
#' @param plot_hmap Boolean value to indicate whether or not to draw the
#' heatmap plot of the scores. (default = TRUE)
#' @param ... Additional arguments for `plot_scores` for aesthetics of the heatmap plot
#'
#' @return Matrix of pathway scores per sample. Columns are samples, rows are
#' pathways. Optionally, displays a heatmap of this matrix.
#'
#' @export
#'
#' @examples
#' score_matrix <- calculate_pw_scores(RA_output, RA_exp_mat, plot_hmap = FALSE)
calculate_pw_scores <- function(pw_table, exp_mat,
                                cases = NULL, plot_hmap = TRUE, ...) {
  if (any(!cases %in% colnames(exp_mat)) & !is.null(cases))
    stop("Missing cases in the expression matrix!")

  all_pws_scores <- c()
  for (i in 1:nrow(pw_table)) {
    # Get DEGs
    up_genes <- pw_table$Up_regulated[i]
    down_genes <- pw_table$Down_regulated[i]
    up_genes <- unlist(strsplit(up_genes, ", "))
    down_genes <- unlist(strsplit(down_genes, ", "))

    genes <- c(up_genes, down_genes)

    # some genes may not be in exp. matrix
    genes <- genes[genes %in% rownames(exp_mat)]

    if (length(genes) != 0) {
      # subset exp. matrix to include only DEGs
      sub_mat <- exp_mat[rownames(exp_mat) %in% genes,, drop = FALSE]

      pw_score_matrix <- c()
      for (gene in genes) {
        gene_vec <- sub_mat[rownames(sub_mat) == gene, ]
        gene_vec <- as.numeric(gene_vec)
        names(gene_vec) <- colnames(sub_mat)

        # calculate mean and sd across samples
        gene_mean <- base::mean(gene_vec)
        gene_sd <- stats::sd(gene_vec)

        gene_scores <- sapply(gene_vec, function(x) (x - gene_mean) / gene_sd)
        pw_score_matrix <- rbind(pw_score_matrix, gene_scores)
        rownames(pw_score_matrix)[nrow(pw_score_matrix)] <- gene
      }

      pw_scores <- apply(pw_score_matrix, 2, base::mean)
      all_pws_scores <- rbind(all_pws_scores, pw_scores)
      rownames(all_pws_scores)[nrow(all_pws_scores)] <- pw_table$Pathway[i]
    }
  }

  if (!is.null(cases)) {
    ## order as cases, then controls
    match1 <- match(cases, colnames(all_pws_scores))
    match2 <- setdiff(1:ncol(all_pws_scores), match1)
    all_pws_scores <- all_pws_scores[, c(match1, match2)]
    if (plot_hmap) {
      heatmap <- plot_scores(score_matrix = all_pws_scores, cases = cases, ...)
      graphics::plot(heatmap)
    }
  } else if (plot_hmap) {
    heatmap <- plot_scores(score_matrix = all_pws_scores, ...)
    graphics::plot(heatmap)
  }

  return(all_pws_scores)
}

#' Plot the Heatmap of Pathway Scores
#'
#' @param score_matrix Matrix of pathway scores per sample. Columns are
#' samples, rows are pathways.
#' @param cases (Optional) A vector of sample names that are cases in the
#' case/control experiment.
#' @param label_cases Boolean value to indicate whether or not to label the
#' cases in the heatmap plot
#' @param case_control_titles A vector of length two for naming of the 'Case'
#' and 'Control' groups (in order) (default = c('Case', 'Control'))
#' @param low a string indicating the color of 'low' values in the score coloring gradient (default = 'green')
#' @param high a string indicating the color of 'high' values in the score coloring gradient (default = 'red')
#'
#' @return A `ggplot2` object containing the heatmap plot. x-axis indicates
#' the samples. y-axis indicates the pathways. "Pathway Score" indicates the
#' pathway score of a sample. If `cases` are provided, the plot is divided
#' into 2 facets, named by the `case_control_titles`.
#'
#' @import ggplot2
#' @export
#'
#' @examples
#' score_mat <- calculate_pw_scores(RA_output, RA_exp_mat, plot_hmap = FALSE)
#' hmap <- plot_scores(score_mat)
plot_scores <- function(score_matrix, cases = NULL, label_cases = TRUE,
                        case_control_titles = c("Case", "Control"), low = "green", high = "red") {
  if (length(case_control_titles) != 2)
    stop("\"case_control_titles\" must contain two elements!")
  if (any(!cases %in% colnames(score_matrix)) & !is.null(cases))
    stop("Missing cases in the score matrix!")

  ## sort according to activity
  if (!is.null(cases)) {
    tmp <- rowMeans(score_matrix[, cases, drop = FALSE])
    score_matrix <- score_matrix[c(which(tmp >= 0), which(tmp < 0)), ]
  }

  ## transform the matrix
  var_names <- list()
  var_names[["Pathway"]] <- factor(rownames(score_matrix), levels = rev(rownames(score_matrix)))
  var_names[["Sample"]] <- factor(colnames(score_matrix), levels = colnames(score_matrix))

  score_df <- expand.grid(var_names, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  scores <- as.vector(score_matrix)
  scores <- data.frame(scores)
  score_df <- cbind(score_df, scores)
  if (!is.null(cases))
    score_df$Type <- factor(ifelse(score_df$Sample %in% cases, case_control_titles[1], case_control_titles[2]), levels = case_control_titles)

  g <- ggplot2::ggplot(score_df, ggplot2::aes_(x = ~Sample, y = ~Pathway))
  g <- g + ggplot2::geom_tile(ggplot2::aes_(fill = ~scores), color = "white")
  g <- g + ggplot2::scale_fill_gradient2(low = low, mid = "black", high = high)
  g <- g + ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                          axis.title.y = ggplot2::element_blank(),
                          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                          legend.title = ggplot2::element_text(size = 10),
                          legend.text = ggplot2::element_text(size = 12))
  g <- g + ggplot2::labs(fill = "Pathway\nscore")
  if (!is.null(cases)) {
    g <- g + ggplot2::facet_grid(~ Type, scales = "free_x", space = "free")
    g <- g + ggplot2::theme(strip.text.x = ggplot2::element_text(size = 12, face="bold"))
  }
  if (!label_cases) {
    g <- g + ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                            axis.ticks.x = ggplot2::element_blank())
  }
  return(g)
}
