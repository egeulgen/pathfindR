#' Create Terms by Genes Heatmap
#'
#' @param result_df A dataframe of pathfindR results that must contain the following columns: \describe{
#'   \item{Term_Description}{Description of the enriched term (necessary if \code{use_description = TRUE})}
#'   \item{ID}{ID of the enriched term (necessary if \code{use_description = FALSE})}
#'   \item{lowest_p}{the highest adjusted-p value of the given term over all iterations}
#'   \item{Up_regulated}{the up-regulated genes in the input involved in the given term's gene set, comma-separated}
#'   \item{Down_regulated}{the down-regulated genes in the input involved in the given term's gene set, comma-separated}
#' }
#' @param genes_df the input data that was used with \code{\link{run_pathfindR}}.
#'   It must be a data frame with 3 columns: \enumerate{
#'   \item Gene Symbol (Gene Symbol)
#'   \item Change value, e.g. log(fold change) (optional)
#'   \item p value, e.g. adjusted p value associated with differential expression
#' } The change values in this data frame are used to color the affected genes
#' @param num_terms Number of top enriched terms to use while creating the plot. Set to \code{NULL} to use
#'  all enriched terms (default = 10)
#' @param use_description Boolean argument to indicate whether term descriptions
#'  (in the 'Term_Description' column) should be used. (default: \code{FALSE})
#' @inheritParams plot_scores
#' @param legend_title legend title (default = 'change')
#' @param sort_terms_by_p boolean to indicate whether to sort terms by 'lowest_p'
#' (\code{TRUE}) or by number of genes (\code{FALSE}) (default = \code{FALSE})
#' @param ... additional arguments for \code{\link{input_processing}} (used if
#' \code{genes_df} is provided)
#'
#' @return a ggplot2 object of a heatmap where rows are enriched terms and
#' columns are involved input genes. If \code{genes_df} is provided, colors of
#' the tiles indicate the change values.
#' @export
#'
#' @examples
#' term_gene_heatmap(example_pathfindR_output, num_terms = 3)
term_gene_heatmap <- function(result_df, genes_df, num_terms = 10, use_description = FALSE,
                              low = "red", mid = "black", high = "green", legend_title = "change", sort_terms_by_p = FALSE,
                              ...) {
  ############ Arg checks
  if (!is.logical(use_description)) {
    stop("`use_description` must either be TRUE or FALSE!")
  }

  ### Set column for term labels
  ID_column <- ifelse(use_description, "Term_Description", "ID")

  if (!is.data.frame(result_df)) {
    stop("`result_df` should be a data frame")
  }

  nec_cols <- c(ID_column, "lowest_p", "Up_regulated", "Down_regulated")
  if (!all(nec_cols %in% colnames(result_df))) {
    stop("`result_df` should have the following columns: ", paste(dQuote(nec_cols),
      collapse = ", "
    ))
  }

  if (!missing(genes_df)) {
    suppressMessages(input_testing(genes_df))
  }

  if (!is.null(num_terms)) {
    if (!is.numeric(num_terms)) {
      stop("`num_terms` should be numeric or NULL")
    }

    if (num_terms < 1) {
      stop("`num_terms` should be > 0 or NULL")
    }
  }

  if (!isColor(low)) {
    stop("`low` should be a valid color")
  }

  if (!isColor(mid)) {
    stop("`mid` should be a valid color")
  }

  if (!isColor(high)) {
    stop("`high` should be a valid color")
  }

  ############ Init prep steps
  result_df <- result_df[order(result_df$lowest_p), ]
  ### select num_terms genes
  if (!is.null(num_terms)) {
    if (num_terms < nrow(result_df)) {
      result_df <- result_df[1:num_terms, ]
    }
  }

  ### process input genes (if provided)
  if (!missing(genes_df)) {
    genes_df <- input_processing(input = genes_df, ...)
  }

  ### parse genes from enrichment results
  parse_genes <- function(vec, idx) {
    return(unname(unlist(strsplit(vec[idx], ", "))))
  }

  up_genes <- apply(result_df, 1, parse_genes, which(colnames(result_df) == "Up_regulated"))
  down_genes <- apply(result_df, 1, parse_genes, which(colnames(result_df) == "Down_regulated"))

  if (length(down_genes) == 0) {
    down_genes <- rep(NA, nrow(result_df))
  }
  if (length(up_genes) == 0) {
    up_genes <- rep(NA, nrow(result_df))
  }

  names(up_genes) <- names(down_genes) <- result_df[, ID_column]

  ############ Create terms-by-genes matrix and order
  all_genes <- unique(c(unlist(up_genes), unlist(down_genes)))
  all_genes <- all_genes[!is.na(all_genes)]
  all_terms <- result_df[, ID_column]

  term_genes_mat <- matrix(0,
    nrow = nrow(result_df), ncol = length(all_genes),
    dimnames = list(all_terms, all_genes)
  )
  for (i in seq_len(nrow(term_genes_mat))) {
    current_term <- rownames(term_genes_mat)[i]
    current_genes <- c(up_genes[[current_term]], down_genes[[current_term]])
    current_genes <- current_genes[!is.na(current_genes)]
    term_genes_mat[i, match(current_genes, colnames(term_genes_mat))] <- 1
  }

  ### Order by column
  term_genes_mat <- term_genes_mat[, order(colSums(term_genes_mat), decreasing = TRUE)]

  ### Order by row
  ordering_func <- function(row) {
    n <- length(row)
    pow <- 2^-(0:(n - 1))
    return(row %*% pow)
  }
  term_genes_mat <- term_genes_mat[order(apply(term_genes_mat, 1, ordering_func),
    decreasing = TRUE
  ), ]

  ### Transform the matrix
  var_names <- list()
  var_names[["Enriched_Term"]] <- factor(rownames(term_genes_mat), levels = rev(rownames(term_genes_mat)))
  var_names[["Symbol"]] <- factor(colnames(term_genes_mat), levels = colnames(term_genes_mat))


  term_genes_df <- expand.grid(var_names, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  value <- as.vector(term_genes_mat)
  value <- data.frame(value)
  term_genes_df <- cbind(term_genes_df, value)
  term_genes_df$value[term_genes_df$value == 0] <- NA

  bg_df <- expand.grid(Enriched_Term = all_terms, Symbol = all_genes)
  if (sort_terms_by_p) {
    bg_df$Enriched_Term <- factor(bg_df$Enriched_Term, levels = rev(result_df[
      ,
      ID_column
    ]))
  } else {
    bg_df$Enriched_Term <- factor(bg_df$Enriched_Term, levels = rev(rownames(term_genes_mat)))
  }


  bg_df$Symbol <- factor(bg_df$Symbol, levels = colnames(term_genes_mat))

  if (!missing(genes_df)) {
    for (i in seq_len(nrow(term_genes_df))) {
      if (!is.na(term_genes_df$value[i])) {
        if (all(genes_df$CHANGE == 1e+06)) {
          term_genes_df$value[i] <- ifelse(term_genes_df$Symbol[i] %in% up_genes[[as.character(term_genes_df$Enriched_Term[i])]],
            1, -1
          )
        } else {
          term_genes_df$value[i] <- genes_df$CHANGE[genes_df$GENE == term_genes_df$Symbol[i]]
        }
      }
    }

    if (all(genes_df$CHANGE == 1e+06)) {
      term_genes_df$value <- factor(term_genes_df$value, levels = c(-1, 1))
    }
  } else {
    for (i in seq_len(nrow(term_genes_df))) {
      if (!is.na(term_genes_df$value[i])) {
        term_genes_df$value[i] <- ifelse(term_genes_df$Symbol[i] %in% unlist(up_genes),
          "up", "down"
        )
      }
    }
  }

  g <- ggplot2::ggplot(bg_df, ggplot2::aes(x = .data$Symbol, y = .data$Enriched_Term))
  g <- g + ggplot2::geom_tile(fill = "white", color = "white")
  g <- g + ggplot2::theme(
    axis.ticks.y = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(
      angle = 90,
      hjust = 1
    ), axis.text.y = ggplot2::element_text(colour = "#000000"), axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank(), panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_blank(), panel.grid.minor.x = ggplot2::element_blank(),
    panel.grid.minor.y = ggplot2::element_blank(), panel.background = ggplot2::element_rect(fill = "#ffffff")
  )
  g <- g + ggplot2::geom_tile(
    data = term_genes_df, ggplot2::aes(fill = .data$value),
    color = "gray60"
  )
  if (!missing(genes_df)) {
    if (all(genes_df$CHANGE == 1e+06)) {
      g <- g + ggplot2::scale_fill_manual(
        values = c(low, high), na.value = "white",
        name = legend_title
      )
    } else {
      g <- g + ggplot2::scale_fill_gradient2(
        low = low, mid = mid, high = high,
        na.value = "white", name = legend_title
      )
    }
  } else {
    g <- g + ggplot2::scale_fill_manual(
      values = c(low, high), na.value = "white",
      name = legend_title
    )
  }
  return(g)
}
