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
#' @section Conceptual Background:
#' For an experiment matrix (containing expression, methylation, etc. values),
#' the rows of which are genes and the columns of which are samples,
#' we denote: \itemize{
#' \item E as a matrix of size \ifelse{html}{\out{m x n}}{\eqn{m \times n}}
#' \item G as the set of all genes in the experiment \ifelse{html}{\out{G = E<sub>i.</sub>,  i &#8712; [1, m]}}{\eqn{G = E_{i\cdot},  \ \ i \in [1, m]}}
#' \item S as the set of all samples in the experiment \ifelse{html}{\out{S = E<sub>.j</sub>,  i &#8712; [1, n]}}{\eqn{S = E_{j\cdot},  \ \ \in [1, n]}}
#' }
#'
#' We next define the gene score matrix GS (the standardized experiment matrix,
#' also of size \ifelse{html}{\out{m x n}}{\eqn{m \times n}}) as:
#'
#' \ifelse{html}{\out{GS<sub>gs</sub> = (E<sub>gs</sub> - &#x113;<sub>g</sub>) / s<sub>g</sub>}}{\eqn{GS_{gs} = \frac{E_{gs} - \bar{e_g}}{s_g}}}
#'
#' where \ifelse{html}{\out{g &#8712; G}}{\eqn{g \in G}}, \ifelse{html}{\out{s &#8712; S}}{\eqn{s \in S}},
#' \ifelse{html}{\out{&#x113;<sub>g</sub>}}{\eqn{\bar{e_g}}} is the mean of
#' all values for gene g and \ifelse{html}{\out{s<sub>g</sub>}}{\eqn{\bar{s_g}}}
#' is the standard deviation of all values for gene g.
#'
#' We next denote T to be a set of terms (where each \ifelse{html}{\out{t &#8712; T}}{\eqn{t \in T}}
#' is a set of term-related genes, i.e.,
#' \ifelse{html}{\out{t = \{g<sub>x</sub>, ..., g<sub>y</sub>\} &sub; G}}{\eqn{t = \{g_x, ..., g_y\} \subset G}})
#' and finally define the agglomerated term scores matrix TS (where rows
#' correspond to genes and columns corresponds to samples s.t. the matrix has size
#' \ifelse{html}{\out{|T| x n}}{\eqn{|T| \times n}}) as:
#'
#' \ifelse{html}{\out{TS<sub>ts</sub> = 1/|t| &#x2211; <sub>g &#8712; t</sub> GS<sub>gs</sub>}}{\eqn{TS_{ts} = \frac{1}{|t|}\sum_{g \in t} GS_{gs}}},
#' where \ifelse{html}{\out{t &#8712; T}}{\eqn{t \in T}} and \ifelse{html}{\out{s &#8712; S}}{\eqn{s \in S}}.
#'
#' @export
#'
#' @examples
#' score_matrix <- score_terms(RA_output, RA_exp_mat, plot_hmap = FALSE)
score_terms <- function(enrichment_table, exp_mat, cases = NULL,
                        use_description = FALSE,
                        plot_hmap = TRUE,
                        ...) {
  #### Argument Checks
  if (!is.logical(use_description)) {
    stop("`use_description` should either be TRUE or FALSE")
  }

  if (!is.logical(plot_hmap)) {
    stop("`plot_hmap` should either be TRUE or FALSE")
  }

  if (!is.data.frame(enrichment_table)) {
    stop("`enrichment_table` should be a data frame of enrichment results")
  }
  ID_column <- ifelse(use_description, "Term_Description", "ID")
  nec_cols <- c(ID_column, "Up_regulated", "Down_regulated")
  if (!all(nec_cols %in% colnames(enrichment_table))) {
    stop("`enrichment_table` should contain all of ",
         paste(dQuote(nec_cols), collapse = ", "))
  }

  if(!is.matrix(exp_mat)) {
    stop("`exp_mat` should be a matrix")
  }

  if (!is.null(cases)) {
    if (!is.atomic(cases)) {
      stop("`cases` should be a vector")
    }

    if (!all(cases %in% colnames(exp_mat))) {
      stop("Missing `cases` in `exp_mat`")
    }
  }

  #### Create score matrix
  all_scores_matrix <- c()
  for (i in base::seq_len(nrow(enrichment_table))) {
    # Get signif. genes
    up_genes <- enrichment_table$Up_regulated[i]
    down_genes <- enrichment_table$Down_regulated[i]
    up_genes <- unlist(strsplit(up_genes, ", "))
    down_genes <- unlist(strsplit(down_genes, ", "))

    genes <- c(up_genes, down_genes)

    # some genes may not be in exp. matrix
    genes <- genes[genes %in% rownames(exp_mat)]

    if (length(genes) != 0) {
      # subset exp. matrix to include only genes
      sub_mat <- exp_mat[rownames(exp_mat) %in% genes, , drop = FALSE]

      current_term_score_matrix <- c()
      for (gene in genes) {
        gene_vec <- sub_mat[rownames(sub_mat) == gene, ]
        gene_vec <- as.numeric(gene_vec)
        names(gene_vec) <- colnames(sub_mat)

        # calculate mean and sd across samples
        gene_mean <- base::mean(gene_vec)
        gene_sd <- stats::sd(gene_vec)

        gene_scores <- vapply(gene_vec, function(x) (x - gene_mean) / gene_sd, 1.2)
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
#' @inheritParams score_terms
#' @param label_samples Boolean value to indicate whether or not to label the
#' samples in the heatmap plot (default = TRUE)
#' @param case_title Naming of the 'Case' group (as in \code{cases}) (default = "Case")
#' @param control_title Naming of the 'Control' group (default = "Control")
#' @param low a string indicating the color of 'low' values in the score coloring gradient (default = 'green')
#' @param mid a string indicating the color of 'mid' values in the score coloring gradient (default = 'black')
#' @param high a string indicating the color of 'high' values in the score coloring gradient (default = 'red')
#'
#' @return A `ggplot2` object containing the heatmap plot. x-axis indicates
#' the samples. y-axis indicates the enriched terms. "Score" indicates the
#' score of the term in a given sample. If \code{cases} are provided, the plot is
#' divided into 2 facets, named by \code{case_title} and \code{control_title}.
#'
#' @import ggplot2
#' @export
#'
#' @examples
#' score_mat <- score_terms(RA_output, RA_exp_mat, plot_hmap = FALSE)
#' hmap <- plot_scores(score_mat)
plot_scores <- function(score_matrix, cases = NULL, label_samples = TRUE,
                        case_title = "Case",
                        control_title = "Control",
                        low = "green", mid = "black", high = "red") {

  #### Argument Checks
  if (!is.matrix(score_matrix)) {
    stop("`score_matrix` should be a matrix")
  }

  if (!is.null(cases)) {
    if (!is.atomic(cases)) {
      stop("`cases` should be a vector")
    }

    if (!all(cases %in% colnames(score_matrix))) {
      stop("Missing `cases` in `score_matrix`")
    }
  }

  if (!is.logical(label_samples)) {
    stop("`label_samples` should be TRUE or FALSE")
  }

  if (!is.character(case_title) | length(case_title) != 1) {
    stop("`case_title` should be a single character value")
  }

  if (!is.character(control_title) | length(control_title) != 1) {
    stop("`control_title` should be a single character value")
  }


  #### Create plot
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
    score_df$Type <- ifelse(score_df$Sample %in% cases, case_title, control_title)
    score_df$Type <- factor(score_df$Type,
                            levels = c(case_title, control_title))
  }

  g <- ggplot2::ggplot(score_df, ggplot2::aes_(x = ~Sample, y = ~Term))
  g <- g + ggplot2::geom_tile(ggplot2::aes_(fill = ~scores), color = "white")
  g <- g + ggplot2::scale_fill_gradient2(low = low, mid = mid, high = high)
  g <- g + ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                          axis.title.y = ggplot2::element_blank(),
                          axis.text.x = ggplot2::element_text(angle = 45,
                                                              hjust = 1),
                          legend.title = ggplot2::element_text(size = 10),
                          legend.text = ggplot2::element_text(size = 12))
  g <- g + ggplot2::labs(fill = "Score")
  if (!is.null(cases)) {
    g <- g + ggplot2::facet_grid(~Type, scales = "free_x", space = "free")
    g <- g + ggplot2::theme(strip.text.x = ggplot2::element_text(size = 12,
                                                                 face = "bold"))
  }
  if (!label_samples) {
    g <- g + ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                            axis.ticks.x = ggplot2::element_blank())
  }
  return(g)
}
