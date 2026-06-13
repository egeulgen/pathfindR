#' Create UpSet Plot of Enriched Terms
#'
#' @inheritParams term_gene_heatmap
#' @param method the option for producing the plot. Options include 'heatmap',
#' 'boxplot' and 'barplot'. (default = 'heatmap')
#'
#' @return UpSet plots are plots of the intersections of sets as a matrix. This
#' function creates a ggplot object of an UpSet plot where the x-axis is the
#' UpSet plot of intersections of enriched terms. By default (i.e.
#' \code{method = 'heatmap'}) the main plot is a heatmap of genes at the
#' corresponding intersections, colored by up/down regulation (if
#' \code{genes_df} is provided, colored by change values). If
#' \code{method = 'barplot'}, the main plot is bar plots of the number of genes
#' at the corresponding intersections. Finally, if \code{method = 'boxplot'} and
#' if \code{genes_df} is provided, then the main plot displays the boxplots of
#' change values of the genes at the corresponding intersections.
#' @export
#'
#' @examples
#' UpSet_plot(example_pathfindR_output)
UpSet_plot <- function(
  result_df, genes_df, num_terms = 10, method = "heatmap", use_description = FALSE,
  low = "red", mid = "black", high = "green", ...
) {
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

  valid_opts <- c("heatmap", "boxplot", "barplot")
  if (!method %in% valid_opts) {
    stop("`method` should be one of` ", paste(dQuote(valid_opts), collapse = ", "))
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

  ########## Init prep steps
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
  all_terms <- result_df[, ID_column]

  term_genes_mat <- matrix(0,
    nrow = nrow(result_df), ncol = length(all_genes),
    dimnames = list(all_terms, all_genes)
  )
  for (i in seq_len(nrow(term_genes_mat))) {
    current_term <- rownames(term_genes_mat)[i]
    current_genes <- c(up_genes[[current_term]], down_genes[[current_term]])
    term_genes_mat[i, match(current_genes, colnames(term_genes_mat))] <- 1
  }

  ### Transform the matrix
  var_names <- list()
  var_names[["Enriched_Term"]] <- factor(rownames(term_genes_mat), levels = rownames(term_genes_mat))
  var_names[["Symbol"]] <- factor(colnames(term_genes_mat), levels = colnames(term_genes_mat))


  term_genes_df <- expand.grid(var_names, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  value <- as.vector(term_genes_mat)
  value <- data.frame(value)
  term_genes_df <- cbind(term_genes_df, value)
  term_genes_df <- term_genes_df[term_genes_df$value != 0, ]

  ### Order according to frequencies
  term_genes_df$Enriched_Term <- factor(term_genes_df$Enriched_Term, levels = names(sort(table(term_genes_df$Enriched_Term),
    decreasing = TRUE
  )))
  term_genes_df$Symbol <- factor(term_genes_df$Symbol, levels = rev(names(sort(table(term_genes_df$Symbol)))))

  terms_lists <- rev(split(term_genes_df$Enriched_Term, term_genes_df$Symbol))

  plot_df <- data.frame(
    Gene = names(terms_lists),
    Up_Down = ifelse(names(terms_lists) %in% unlist(up_genes), "up", "down"),
    stringsAsFactors = FALSE
  )

  plot_df$Term <- terms_lists

  bg_df <- expand.grid(Gene = unique(plot_df$Gene), Term = unique(plot_df$Term))

  if (method == "heatmap") {
    g <- ggplot2::ggplot(bg_df, ggplot2::aes(x = .data$Term, y = .data$Gene))
    g <- g + ggplot2::geom_tile(fill = "white", color = "gray60")

    if (missing(genes_df)) {
      g <- g + ggplot2::geom_tile(data = plot_df, ggplot2::aes(
        x = .data$Term,
        y = .data$Gene, fill = .data$Up_Down
      ), color = "gray60")
      g <- g + ggplot2::scale_fill_manual(values = c(low, high))
    } else {
      plot_df$Value <- genes_df$CHANGE[match(names(plot_df$Term), genes_df$GENE)]
      g <- g + ggplot2::geom_tile(data = plot_df, ggplot2::aes(
        x = .data$Term,
        y = .data$Gene, fill = .data$Value
      ), color = "gray60")
      g <- g + ggplot2::scale_fill_gradient2(low = low, mid = mid, high = high)
    }
    g <- g + ggplot2::theme_minimal()
    g <- g + ggplot2::theme(
      axis.title = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(), legend.title = ggplot2::element_blank()
    )
  } else if (method == "boxplot") {
    if (missing(genes_df)) {
      stop("For `method = boxplot`, you must provide `genes_df`")
    }

    plot_df$Value <- genes_df$CHANGE[match(names(plot_df$Term), genes_df$GENE)]
    g <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$Term, y = .data$Value))
    g <- g + ggplot2::geom_boxplot()
    g <- g + ggplot2::geom_jitter(width = 0.1)
  } else {
    g <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$Term))
    g <- g + ggplot2::geom_bar()
  }

  g <- g + ggupset::scale_x_upset(order_by = ifelse(missing(genes_df), "freq",
    "degree"
  ), reverse = !missing(genes_df))
  return(g)
}
