#' Check if value is a valid color
#'
#' @param x value
#'
#' @return TRUE if x is a valid color, otherwise FALSE
isColor <- function(x) {
  if (!is.character(x) | length(x) != 1) {
    return(FALSE)
  }
  tryCatch(is.matrix(grDevices::col2rgb(x)), error = function(e) FALSE)
}

#' Order input data frame by provided columnn
#'
#' @param df the input data frame to be ordered
#' @param order_by A column name
#'
#' @return The ordered data frame or raises error
order_df_by_columnn <- function(df, order_by) {
  if (!c(order_by %in% colnames(df))) {
    stop("`order_by` column doesn't exist in `df`")
  }
  col_values <- df[[order_by]]
  if (anyNA(col_values)) {
    stop("Column values of `order_by` cannot have NAs!")
  }
  result_df <- tryCatch(
    {
      df[order(df[[order_by]], decreasing = FALSE), ]
    },
    error = function(e) {
      stop(
        sprintf(
          "`order_by` (%s) cannot be used to order the `df`: %s",
          order_by,
          e$message
        ),
        call. = FALSE
      )
    }
  )
  return(result_df)
}

#' Color hsa KEGG pathway
#'
#' @param pw_id hsa KEGG pathway id (e.g. hsa05012)
#' @param change_vec vector of change values, names should be hsa KEGG gene ids
#' @param scale_vals should change values be scaled? (default = \code{TRUE})
#' @param node_cols low, middle and high color values for coloring the pathway nodes
#' (default = \code{NULL}). If \code{node_cols=NULL}, the low, middle and high color
#' are set as 'green', 'gray' and 'red'. If all change values are 1e6 (in case no
#' changes are supplied, this dummy value is assigned by
#' \code{\link{input_processing}}), only one color ('#F38F18' if NULL) is used.
#' @inheritParams ggplot2::theme
#'
#' @return a ggplot object containing the colored KEGG pathway diagram visualization
#'
#' @examples
#' \dontrun{
#' pw_id <- "hsa00010"
#' change_vec <- c(-2, 4, 6)
#' names(change_vec) <- c("hsa:2821", "hsa:226", "hsa:229")
#' result <- pathfindR:::color_kegg_pathway(pw_id, change_vec)
#' }
color_kegg_pathway <- function(pw_id, change_vec, scale_vals = TRUE, node_cols = NULL, legend.position = "top") {
  ############ Arg checks
  if (!is.logical(scale_vals)) {
    stop("`scale_vals` should be logical")
  }

  ## check node_cols
  if (!is.null(node_cols)) {
    if (!is.atomic(node_cols)) {
      stop("`node_cols` should be a vector of colors")
    }

    if (!all(change_vec == 1e+06) & length(node_cols) != 3) {
      stop("the length of `node_cols` should be 3")
    }

    if (!all(vapply(node_cols, isColor, TRUE))) {
      stop("`node_cols` should be a vector of valid colors")
    }
  }
  ############ Set node palette if node_cols not supplied, use default
  ############ color(s)
  if (!is.null(node_cols)) {
    if (all(change_vec == 1e+06)) {
      message("all `change_vec` values are 1e6, using the first color in `node_cols`")
      low_col <- mid_col <- high_col <- node_cols[1]
    } else {
      low_col <- node_cols[1]
      mid_col <- node_cols[2]
      high_col <- node_cols[3]
    }
  } else if (all(change_vec == 1e+06)) {
    ## NO CHANGES SUPPLIED
    low_col <- mid_col <- high_col <- "#F38F18"
  } else {
    low_col <- "red"
    mid_col <- "gray"
    high_col <- "green"
  }

  ############ Assign the input change values to any corresponding pathway gene nodes
  # create pathway graph object and collect all pathway genes
  ggkegg_temp_dir <- file.path(tempdir(check = TRUE), "ggkegg")
  dir.create(ggkegg_temp_dir, showWarnings = FALSE)

  g <- tryCatch(
    {
      ggkegg::pathway(pid = pw_id, directory = ggkegg_temp_dir)
    },
    error = function(e) {
      message(paste("Cannot parse KEGG pathway for:", pw_id))
      message("Here's the original error message:")
      message(e$message)
      return(NULL)
    },
    warning = function(w) {
      message(paste("Cannot parse KEGG pathway for:", pw_id))
      message("Here's the original error message:")
      message(w$message)
      return(NULL)
    }
  )

  if (is.null(g)) {
    return(NULL)
  }

  gene_nodes <- names(igraph::V(g))[igraph::V(g)$type == "gene"]

  ## aggregate change values over all pathway gene nodes
  pw_vis_changes <- c()
  for (i in seq_len(length(gene_nodes))) {
    node_name <- gene_nodes[i]
    node <- unlist(strsplit(node_name, " "))
    cond <- names(change_vec) %in% node

    if (any(cond)) {
      node_val <- mean(change_vec[cond])
      names(node_val) <- node_name
      pw_vis_changes <- c(pw_vis_changes, node_val)
    }
  }
  ## if no input genes present in chosen pathway
  if (all(is.na(pw_vis_changes))) {
    return(NULL)
  }

  ############ Determine node colors
  ### scaling
  if (!all(pw_vis_changes == 1e+06) & scale_vals) {
    common_limit <- max(abs(pw_vis_changes))
    pw_vis_changes <- ifelse(pw_vis_changes < 0,
      -abs(pw_vis_changes) / common_limit,
      pw_vis_changes / common_limit
    )
  }


  ############ Create pathway diagram visualisation
  igraph::V(g)$change_value <- NA
  igraph::V(g)$change_value[match(names(pw_vis_changes), names(igraph::V(g)))] <- pw_vis_changes

  p <- ggraph::ggraph(g, layout = "manual", x = igraph::V(g)$x, y = igraph::V(g)$y)
  p <- p + ggkegg::geom_node_rect(ggplot2::aes(filter = !is.na(.data$change_value), fill = .data$change_value))
  p <- p + ggkegg::overlay_raw_map(pw_id)
  p <- p + ggplot2::scale_fill_gradient2(low = low_col, mid = mid_col, high = high_col)
  p <- p + ggplot2::theme_void()
  p <- p + ggplot2::theme(
    legend.title = ggplot2::element_blank(),
    legend.position = legend.position
  )

  return(p)
}
