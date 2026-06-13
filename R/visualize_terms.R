#' Create Diagrams for Enriched Terms
#'
#' @param result_df Data frame of enrichment results. Must-have columns for
#'  KEGG human pathway diagrams (\code{is_KEGG_result = TRUE}) are: 'ID' and 'Term_Description'.
#'  Must-have columns for the rest are: 'Term_Description', 'Up_regulated' and
#' 'Down_regulated'
#' @param input_processed input data processed via \code{\link{input_processing}},
#'  not necessary when \code{is_KEGG_result = FALSE}
#' @param is_KEGG_result boolean to indicate whether KEGG gene sets were used for
#'  enrichment analysis or not (default = \code{TRUE})
#' @inheritParams return_pin_path
#' @param ... additional arguments for \code{\link{visualize_KEGG_diagram}} (used
#' when \code{is_KEGG_result = TRUE}) or \code{\link{visualize_term_interactions}}
#' (used when \code{is_KEGG_result = FALSE})
#'
#' @return Depending on the argument \code{is_KEGG_result}, creates visualization of
#'  interactions of genes involved in the list of enriched terms in
#'  \code{result_df}. Returns a list of ggplot objects named by Term ID.
#'
#'
#' @details For \code{is_KEGG_result = TRUE}, KEGG pathway diagrams are created,
#' affected nodes colored by up/down regulation status.
#' For other gene sets, interactions of affected genes are determined (via a shortest-path
#' algorithm) and are visualized (colored by change status) using igraph.
#'
#'
#' @export
#'
#' @seealso See \code{\link{visualize_KEGG_diagram}} for the visualization function
#' of KEGG diagrams. See \code{\link{visualize_term_interactions}} for the
#' visualization function that generates diagrams showing the interactions of
#' input genes in the PIN. See \code{\link{run_pathfindR}} for the wrapper
#' function of the pathfindR workflow.
#'
#' @examples
#' \dontrun{
#' input_processed <- data.frame(
#'   GENE = c("PARP1", "NDUFA1", "STX6", "SNAP23"),
#'   CHANGE = c(1.5, -2, 3, 5)
#' )
#' result_df <- example_pathfindR_output[1:2, ]
#'
#' gg_list <- visualize_terms(result_df, input_processed)
#' gg_list2 <- visualize_terms(result_df, is_KEGG_result = FALSE, pin_name_path = "IntAct")
#' }
visualize_terms <- function(
  result_df, input_processed = NULL, is_KEGG_result = TRUE, pin_name_path = "Biogrid", ...
) {
  ############ Argument Checks
  if (!is.data.frame(result_df)) {
    stop("`result_df` should be a data frame")
  }

  if (!is.logical(is_KEGG_result)) {
    stop("the argument `is_KEGG_result` should be either TRUE or FALSE")
  }

  if (is_KEGG_result) {
    nec_cols <- "ID"
  } else {
    nec_cols <- c("Term_Description", "Up_regulated", "Down_regulated")
  }
  if (!all(nec_cols %in% colnames(result_df))) {
    stop("`result_df` should contain the following columns: ", paste(dQuote(nec_cols),
      collapse = ", "
    ))
  }

  if (is_KEGG_result) {
    if (is.null(input_processed)) {
      stop("`input_processed` should be specified when `is_KEGG_result = TRUE`")
    }
  }

  ############ Generate Diagrams
  if (is_KEGG_result) {
    visualize_KEGG_diagram(
      kegg_pw_ids = result_df$ID, input_processed = input_processed, ...
    )
  } else {
    visualize_term_interactions(
      result_df = result_df, pin_name_path = pin_name_path, ...
    )
  }
}

#' Visualize Interactions of Genes Involved in the Given Enriched Terms
#'
#' @param result_df Data frame of enrichment results. Must-have columns
#' are: 'Term_Description', 'Up_regulated' and 'Down_regulated'
#' @inheritParams return_pin_path
#' @param show_legend Boolean to indicate whether to display the legend (\code{TRUE})
#' or not (\code{FALSE}) (default: \code{TRUE})
#'
#' @return list of ggplot objects (named by Term ID) visualizing the interactions of genes involved
#' in the given enriched terms (annotated in the \code{result_df}) in the PIN used
#' for enrichment analysis (specified by \code{pin_name_path}).
#'
#' @details The following steps are performed for the visualization of interactions
#' of genes involved for each enriched term: \enumerate{
#'   \item shortest paths between all affected genes are determined (via \code{\link[igraph]{igraph}})
#'   \item the nodes of all shortest paths are merged
#'   \item the PIN is subsetted using the merged nodes (genes)
#'   \item using the PIN subset, the graph showing the interactions is generated
#'   \item the final graph is visualized using \code{\link[igraph]{igraph}}, colored by changed
#'   status (if provided)
#' }
#'
#' @export
#'
#' @seealso See \code{\link{visualize_terms}} for the wrapper function
#'   for creating enriched term diagrams. See \code{\link{run_pathfindR}} for the
#'   wrapper function of the pathfindR enrichment workflow.
#'
#' @examples
#' \dontrun{
#' result_df <- example_pathfindR_output[1:2, ]
#' gg_list <- visualize_term_interactions(result_df, pin_name_path = "IntAct")
#' }
visualize_term_interactions <- function(result_df, pin_name_path, show_legend = TRUE) {
  ############ Initial Steps fix naming issue
  result_df$Term_Description <- gsub("\\/", "-", result_df$Term_Description)

  ## load PIN
  pin_path <- return_pin_path(pin_name_path)
  pin <- utils::read.delim(file = pin_path, header = FALSE)
  pin$V2 <- NULL

  pin[, 1] <- base::toupper(pin[, 1])
  pin[, 2] <- base::toupper(pin[, 2])

  ## pin graph
  pin_g <- igraph::graph_from_data_frame(pin, directed = FALSE)

  ############ Visualize interactions by enriched term
  pw_vis_list <- list()
  for (i in base::seq_len(nrow(result_df))) {
    current_row <- result_df[i, ]

    up_genes <- base::toupper(unlist(strsplit(current_row$Up_regulated, ", ")))
    down_genes <- base::toupper(unlist(strsplit(current_row$Down_regulated, ", ")))
    current_genes <- c(down_genes, up_genes)

    ## Add active snw genes if listed
    if (!is.null(result_df$non_Signif_Snw_Genes)) {
      snw_genes <- unlist(strsplit(current_row$non_Signif_Snw_Genes, ", "))
      snw_genes <- base::toupper(snw_genes)
      current_genes <- c(current_genes, snw_genes)
    } else {
      snw_genes <- NULL
    }

    if (length(current_genes) < 2) {
      message(paste0("< 2 genes, skipping visualization of ", current_row$Term_Description))
    } else {
      cat("Visualizing:", paste0("(", i, ")"), current_row$Term_Description, paste(rep(" ", 200),
        collapse = ""
      ), "\r")

      ## Find genes without direct interaction
      cond1 <- pin$V1 %in% current_genes
      cond2 <- pin$V3 %in% current_genes
      direct_interactions <- pin[cond1 & cond2, ]
      tmp <- c(direct_interactions$V1, direct_interactions$V3)
      missing_genes <- current_genes[!current_genes %in% tmp]

      ## Find shortest path between genes without direct interaction and
      ## other current_genes
      s_path_genes <- c()
      for (gene in missing_genes) {
        tmp <- suppressWarnings(igraph::shortest_paths(pin_g,
          from = which(names(igraph::V(pin_g)) ==
            gene), to = which(names(igraph::V(pin_g)) %in% current_genes),
          output = "vpath"
        ))
        tmp <- unique(unlist(lapply(tmp$vpath, function(x) names(x))))
        s_path_genes <- unique(c(s_path_genes, tmp))
      }

      final_genes <- unique(c(current_genes, s_path_genes))
      cond1 <- pin$V1 %in% final_genes
      cond2 <- pin$V3 %in% final_genes
      final_interactions <- pin[cond1 & cond2, ]
      g <- igraph::graph_from_data_frame(final_interactions, directed = FALSE)

      cond1 <- names(igraph::V(g)) %in% up_genes
      cond2 <- names(igraph::V(g)) %in% down_genes
      cond3 <- names(igraph::V(g)) %in% snw_genes
      node_type <- as.factor(ifelse(cond1, "up",
        ifelse(cond2, "down",
          ifelse(cond3,
            "interactor", "none"
          )
        )
      ))
      igraph::V(g)$type <- node_type

      node_colors <- c("green", "red", "blue", "gray")
      names(node_colors) <- c("up", "down", "interactor", "none")
      node_colors <- node_colors[levels(node_type)]

      type_descriptions <- c(
        none = "other", up = "up-regulated gene", down = "down-regulated gene", interactor = "interacting non-input gene"
      )
      type_descriptions <- type_descriptions[levels(node_type)]

      p <- ggraph::ggraph(g, layout = "stress")
      p <- p + ggraph::geom_edge_link(alpha = 0.8, colour = "darkgrey", linewidth = 0.5)
      p <- p + ggraph::geom_node_point(ggplot2::aes(color = .data$type), size = 5)
      p <- p + ggplot2::theme_void()
      p <- p + suppressWarnings(ggraph::geom_node_text(ggplot2::aes(label = .data$name),
        nudge_y = 0.2, repel = TRUE, max.overlaps = 20
      ))
      p <- p + ggplot2::scale_color_manual(
        values = node_colors, name = NULL,
        labels = type_descriptions
      )
      p <- p + ggplot2::ggtitle(
        paste(current_row$Term_Description, "\n Involved Gene Interactions in", pin_name_path)
      )
      pw_vis_list[[current_row$ID]] <- p
    }
  }
  return(pw_vis_list)
}

#' Visualize Human KEGG Pathways
#'
#' @param kegg_pw_ids KEGG ids of pathways to be colored and visualized
#' @param input_processed input data processed via \code{\link{input_processing}}
#' @inheritParams color_kegg_pathway
#'
#' @return Creates colored visualizations of the enriched human KEGG pathways
#' and returns them as a list of ggplot objects, named by Term ID.
#'
#' @seealso See \code{\link{visualize_terms}} for the wrapper function for
#' creating enriched term diagrams. See \code{\link{run_pathfindR}} for the
#' wrapper function of the pathfindR enrichment workflow.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' input_processed <- data.frame(
#'   GENE = c("PKLR", "GPI", "CREB1", "INS"),
#'   CHANGE = c(1.5, -2, 3, 5)
#' )
#' gg_list <- visualize_KEGG_diagram(c("hsa00010", "hsa04911"), input_processed)
#' }
visualize_KEGG_diagram <- function(
  kegg_pw_ids,
  input_processed,
  scale_vals = TRUE,
  node_cols = NULL,
  legend.position = "top"
) {
  message("This function utilises one functionality of `ggkegg`. For more options, visit https://github.com/noriakis/ggkegg")

  ############ Arg checks

  ### kegg_pw_ids
  if (!is.atomic(kegg_pw_ids)) {
    stop("`kegg_pw_ids` should be a vector of KEGG IDs")
  }
  if (!all(grepl("^[a-z]{3}[0-9]{5}$", kegg_pw_ids))) {
    stop("`kegg_pw_ids` should be a vector of valid hsa KEGG IDs")
  }

  ### input_processed
  if (!is.data.frame(input_processed)) {
    stop("`input_processed` should be a data frame")
  }

  nec_cols <- c("GENE", "CHANGE")
  if (!all(nec_cols %in% colnames(input_processed))) {
    stop("`input_processed` should contain the following columns: ", paste(dQuote(nec_cols),
      collapse = ", "
    ))
  }

  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    message(
      "Package 'org.Hs.eg.db' is not installed; returning empty list.\n",
      "Install it with:\n",
      "  if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')\n",
      "  BiocManager::install('org.Hs.eg.db')"
    )
    return(list())
  }

  ############ Create change vector Convert gene symbols into NCBI gene IDs
  tmp <- AnnotationDbi::mget(input_processed$GENE, AnnotationDbi::revmap(org.Hs.eg.db::org.Hs.egSYMBOL),
    ifnotfound = NA
  )
  input_processed$EG_ID <- vapply(tmp, function(x) as.character(x[1]), "EGID")
  input_processed <- input_processed[!is.na(input_processed$EG_ID), ]

  ### A rule of thumb for the 'kegg' ID is entrezgene ID for eukaryote
  ### species
  input_processed$KEGG_ID <- paste0("hsa:", input_processed$EG_ID)

  ############ Fetch all pathway genes, create vector of change values and
  ############ Generate colored pathway diagrams for each pathway
  change_vec <- input_processed$CHANGE
  names(change_vec) <- input_processed$KEGG_ID

  cat("Generating pathway diagrams of", length(kegg_pw_ids), "KEGG pathways\n\n")
  pw_vis_list <- lapply(
    kegg_pw_ids,
    color_kegg_pathway,
    change_vec = change_vec,
    scale_vals = scale_vals,
    node_cols = node_cols,
    legend.position = legend.position
  )
  names(pw_vis_list) <- kegg_pw_ids

  return(pw_vis_list)
}
