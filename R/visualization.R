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
#' gg_list2 <- visualize_terms(result_df, is_KEGG_result = FALSE, pin_name_path = 'IntAct')
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
            collapse = ", "))
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
#' gg_list <- visualize_term_interactions(result_df, pin_name_path = 'IntAct')
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
            cat("Visualizing:", paste0("(", i, ")") , current_row$Term_Description, paste(rep(" ", 200),
                collapse = ""), "\r")

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
                tmp <- suppressWarnings(igraph::shortest_paths(pin_g, from = which(names(igraph::V(pin_g)) ==
                  gene), to = which(names(igraph::V(pin_g)) %in% current_genes),
                  output = "vpath"))
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
                                                        "interactor", "none"))))
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
                                                             nudge_y = 0.2, repel = TRUE, max.overlaps = 20))
            p <- p + ggplot2::scale_color_manual(values = node_colors, name = NULL,
                                                 labels = type_descriptions)
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
            collapse = ", "))
    }

    ############ Create change vector Convert gene symbols into NCBI gene IDs
    tmp <- AnnotationDbi::mget(input_processed$GENE, AnnotationDbi::revmap(org.Hs.eg.db::org.Hs.egSYMBOL),
        ifnotfound = NA)
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
      change_vec=change_vec,
      scale_vals = scale_vals,
      node_cols = node_cols,
      legend.position = legend.position
    )
    names(pw_vis_list) <- kegg_pw_ids

    return(pw_vis_list)
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
#' pw_id <- 'hsa00010'
#' change_vec <- c(-2, 4, 6)
#' names(change_vec) <- c('hsa:2821', 'hsa:226', 'hsa:229')
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

    g <- tryCatch({
      ggkegg::pathway(pid = pw_id, directory = tempdir(), use_cache = FALSE)
    }, error = function(e) {
      message(paste("Cannot parse KEGG pathway for:", pw_id))
      message("Here's the original error message:")
      message(e$message)
      return(NULL)
    }, warning = function(w) {
      message(paste("Cannot parse KEGG pathway for:", pw_id))
      message("Here's the original error message:")
      message(w$message)
      return(NULL)
    })

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
                               pw_vis_changes / common_limit)
    }


    ############ Create pathway diagram visualisation
    igraph::V(g)$change_value <- NA
    igraph::V(g)$change_value[match(names(pw_vis_changes), names(igraph::V(g)))] <- pw_vis_changes

    p <- ggraph::ggraph(g, layout="manual", x=igraph::V(g)$x, y=igraph::V(g)$y)
    p <- p + ggkegg::geom_node_rect(ggplot2::aes(filter = !is.na(.data$change_value), fill = .data$change_value))
    p <- p + ggkegg::overlay_raw_map(pw_id, directory = tempdir(), use_cache = FALSE)
    p <- p + ggplot2::scale_fill_gradient2(low = low_col, mid = mid_col, high = high_col)
    p <- p + ggplot2::theme_void()
    p <- p + ggplot2::theme(
      legend.title = ggplot2::element_blank(),
      legend.position = legend.position
    )

    return(p)
}

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
    num_bubbles = 4, even_breaks = TRUE) {
    message("Plotting the enrichment bubble chart")
    necessary <- c("Term_Description", "Fold_Enrichment", "lowest_p", "Up_regulated",
        "Down_regulated")

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

    # sort by lowest adj.p
    result_df <- result_df[order(result_df$lowest_p), ]

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

    num_genes <- vapply(result_df$Up_regulated, function(x) length(unlist(strsplit(x,
        ", "))), 1)
    num_genes <- num_genes + vapply(result_df$Down_regulated, function(x) length(unlist(strsplit(x,
        ", "))), 1)

    result_df$Term_Description <- factor(result_df$Term_Description, levels = rev(unique(result_df$Term_Description)))

    log_p <- -log10(result_df$lowest_p)

    g <- ggplot2::ggplot(result_df, ggplot2::aes(.data$Fold_Enrichment, .data$Term_Description))
    g <- g + ggplot2::geom_point(ggplot2::aes(color = log_p, size = num_genes), na.rm = TRUE)
    g <- g + ggplot2::theme_bw()
    g <- g + ggplot2::theme(axis.text.x = ggplot2::element_text(size = 10), axis.text.y = ggplot2::element_text(size = 10),
        plot.title = ggplot2::element_blank())
    g <- g + ggplot2::xlab("Fold Enrichment")
    g <- g + ggplot2::theme(axis.title.y = ggplot2::element_blank())
    g <- g + ggplot2::labs(size = "# genes", color = expression(-log[10](p)))

    ## breaks for # genes
    if (max(num_genes) < num_bubbles) {
        g <- g + ggplot2::scale_size_continuous(breaks = seq(0, max(num_genes)))
    } else {
        if (even_breaks) {
            brks <- base::seq(0, max(num_genes), round(max(num_genes)/(num_bubbles +
                1)))
        } else {
            brks <- base::round(base::seq(0, max(num_genes), length.out = num_bubbles +
                1))
        }
        g <- g + ggplot2::scale_size_continuous(breaks = brks)
    }

    g <- g + ggplot2::scale_color_gradient(low = "#f5efef", high = "red")

    if (plot_by_cluster & "Cluster" %in% colnames(result_df)) {
        g <- g + ggplot2::facet_grid(result_df$Cluster ~ ., scales = "free_y", space = "free",
            drop = TRUE)
    } else if (plot_by_cluster) {
        message("For plotting by cluster, there must a column named `Cluster` in the input data frame!")
    }

    return(g)
}


#' Create Term-Gene Graph
#'
#' @param result_df A dataframe of pathfindR results that must contain the following columns: \describe{
#'   \item{Term_Description}{Description of the enriched term (necessary if \code{use_description = TRUE})}
#'   \item{ID}{ID of the enriched term (necessary if \code{use_description = FALSE})}
#'   \item{lowest_p}{the lowest adjusted-p value of the given term over all iterations}
#'   \item{Up_regulated}{the up-regulated genes in the input involved in the given term's gene set, comma-separated}
#'   \item{Down_regulated}{the down-regulated genes in the input involved in the given term's gene set, comma-separated}
#' }
#' @param num_terms Number of top enriched terms to use while creating the graph. Set to \code{NULL} to use
#'  all enriched terms (default = 10, i.e. top 10 terms)
#' @param layout The type of layout to create (see \code{\link[ggraph]{ggraph}} for details. Default = 'stress')
#' @param use_description Boolean argument to indicate whether term descriptions
#'  (in the 'Term_Description' column) should be used. (default = \code{FALSE})
#' @param node_size Argument to indicate whether to use number of significant genes ('num_genes')
#'  or the -log10(lowest p value) ('p_val') for adjusting the node sizes (default = 'num_genes')
#' @param node_colors vector of 3 colors to be used for coloring nodes (colors for term nodes, up, and down, respectively)
#'
#' @return a  \code{\link[ggraph]{ggraph}} object containing the term-gene graph.
#'  Each node corresponds to an enriched term (beige), an up-regulated gene (green)
#'  or a down-regulated gene (red). An edge between a term and a gene indicates
#'  that the given term involves the gene. Size of a term node is proportional
#'  to either the number of genes (if \code{node_size = 'num_genes'}) or
#'  the -log10(lowest p value) (if \code{node_size = 'p_val'}).
#'
#' @details This function (adapted from the Gene-Concept network visualization
#' by the R package \code{enrichplot}) can be utilized to visualize which input
#' genes are involved in the enriched terms as a graph. The term-gene graph
#' shows the links between genes and biological terms and allows for the
#' investigation of multiple terms to which significant genes are related. The
#' graph also enables determination of the overlap between the enriched terms
#' by identifying shared and distinct significant term-related genes.
#'
#' @import ggraph
#' @export
#'
#' @examples
#' p <- term_gene_graph(example_pathfindR_output)
#' p <- term_gene_graph(example_pathfindR_output, num_terms = 5)
#' p <- term_gene_graph(example_pathfindR_output, node_size = 'p_val')
term_gene_graph <- function(result_df, num_terms = 10, layout = "stress", use_description = FALSE,
    node_size = "num_genes", node_colors = c("#E5D7BF", "green", "red")) {
    ############ Argument Checks Check num_terms is NULL or numeric
    if (!is.numeric(num_terms) & !is.null(num_terms)) {
        stop("`num_terms` must either be numeric or NULL!")
    }

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

    if (!is.data.frame(result_df)) {
        stop("`result_df` should be a data frame")
    }

    ### Check necessary columnns
    necessary_cols <- c(ID_column, "lowest_p", "Up_regulated", "Down_regulated")

    if (!all(necessary_cols %in% colnames(result_df))) {
        stop(paste(c("All of", paste(necessary_cols, collapse = ", "), "must be present in `results_df`!"),
            collapse = " "))
    }

    if (!is.atomic(node_colors)) {
      stop("`node_colors` should be a vector of colors")
    }

    if (!all(vapply(node_colors, isColor, TRUE))) {
      stop("`node_colors` should be a vector of valid colors")
    }

    if (length(node_colors) != 3) {
      stop("`node_colors` must contain exactly 3 colors")
    }

    ############ Initial steps set num_terms to NULL if number of enriched
    ############ terms is smaller than num_terms
    if (!is.null(num_terms)) {
        if (nrow(result_df) < num_terms) {
            num_terms <- NULL
        }
    }

    ### Order and filter for top N genes
    result_df <- result_df[order(result_df$lowest_p, decreasing = FALSE), ]
    if (!is.null(num_terms)) {
        result_df <- result_df[1:num_terms, ]
    }

    ### Prep data frame for graph
    graph_df <- data.frame()
    for (i in base::seq_len(nrow(result_df))) {
        up_genes <- unlist(strsplit(result_df$Up_regulated[i], ", "))
        down_genes <- unlist(strsplit(result_df$Down_regulated[i], ", "))
        for (gene in c(up_genes, down_genes)) {
            graph_df <- rbind(graph_df, data.frame(Term = result_df[i, ID_column],
                Gene = gene))
        }
    }

    up_genes <- lapply(result_df$Up_regulated, function(x) unlist(strsplit(x, ", ")))
    up_genes <- unlist(up_genes)

    ############ Create graph object and plot create igraph object
    g <- igraph::graph_from_data_frame(graph_df, directed = FALSE)
    cond_term <- names(igraph::V(g)) %in% result_df[, ID_column]
    cond_up_gene <- names(igraph::V(g)) %in% up_genes

    node_type <-  ifelse(cond_term, "term", ifelse(cond_up_gene, "up", "down"))
    node_type <- factor(node_type, levels = c("term", "up", "down"))
    node_type <- droplevels(node_type)
    igraph::V(g)$type <- node_type

    type_descriptions <- c(term="enriched term", up="up-regulated gene", down="down-regulated gene")
    type_descriptions <- type_descriptions[levels(node_type)]

    names(node_colors) <- c("term", "up", "down")
    node_colors <- node_colors[levels(node_type)]

    # Adjust node sizes
    if (node_size == "num_genes") {
        sizes <- igraph::degree(g)
        sizes <- ifelse(igraph::V(g)$type == "term", sizes, 2)
        size_label <- "# genes"
    } else {
        idx <- match(names(igraph::V(g)), result_df[, ID_column])
        sizes <- -log10(result_df$lowest_p[idx])
        sizes[is.na(sizes)] <- 2
        size_label <- "-log10(p)"
    }
    igraph::V(g)$size <- sizes
    igraph::V(g)$label.cex <- 0.5
    igraph::V(g)$frame.color <- "gray"

    ### Create graph
    p <- ggraph::ggraph(g, layout = layout)
    p <- p + ggraph::geom_edge_link(alpha = 0.8, colour = "darkgrey")
    p <- p + ggraph::geom_node_point(ggplot2::aes(color = .data$type, size = .data$size))
    p <- p + ggplot2::scale_size(range = c(5, 10), breaks = round(seq(round(min(igraph::V(g)$size)),
        round(max(igraph::V(g)$size)), length.out = 4)), name = size_label)
    p <- p + ggplot2::theme_void()
    p <- p + suppressWarnings(ggraph::geom_node_text(ggplot2::aes(label = .data$name),
        nudge_y = 0.2, repel = TRUE, max.overlaps = 20))
    p <- p + ggplot2::scale_color_manual(values = node_colors, name = NULL,
        labels = type_descriptions)
    if (is.null(num_terms)) {
        p <- p + ggplot2::ggtitle("Term-Gene Graph")
    } else {
        p <- p + ggplot2::ggtitle("Term-Gene Graph", subtitle = paste(c("Top", num_terms,
            "terms"), collapse = " "))
    }

    p <- p + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), plot.subtitle = ggplot2::element_text(hjust = 0.5))

    return(p)
}


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
#' @inheritParams term_gene_graph
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
            collapse = ", "))
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

    term_genes_mat <- matrix(0, nrow = nrow(result_df), ncol = length(all_genes),
        dimnames = list(all_terms, all_genes))
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
        decreasing = TRUE), ]

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
        bg_df$Enriched_Term <- factor(bg_df$Enriched_Term, levels = rev(result_df[,
            ID_column]))
    } else {
        bg_df$Enriched_Term <- factor(bg_df$Enriched_Term, levels = rev(rownames(term_genes_mat)))
    }


    bg_df$Symbol <- factor(bg_df$Symbol, levels = colnames(term_genes_mat))

    if (!missing(genes_df)) {
        for (i in seq_len(nrow(term_genes_df))) {
            if (!is.na(term_genes_df$value[i])) {
                if (all(genes_df$CHANGE == 1e+06)) {
                  term_genes_df$value[i] <- ifelse(term_genes_df$Symbol[i] %in% up_genes[[as.character(term_genes_df$Enriched_Term[i])]],
                    1, -1)
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
                  "up", "down")
            }
        }
    }

    g <- ggplot2::ggplot(bg_df, ggplot2::aes(x = .data$Symbol, y = .data$Enriched_Term))
    g <- g + ggplot2::geom_tile(fill = "white", color = "white")
    g <- g + ggplot2::theme(axis.ticks.y = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(angle = 90,
        hjust = 1), axis.text.y = ggplot2::element_text(colour = "#000000"), axis.title.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(), panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_blank(), panel.grid.minor.x = ggplot2::element_blank(),
        panel.grid.minor.y = ggplot2::element_blank(), panel.background = ggplot2::element_rect(fill = "#ffffff"))
    g <- g + ggplot2::geom_tile(data = term_genes_df, ggplot2::aes(fill = .data$value),
        color = "gray60")
    if (!missing(genes_df)) {
        if (all(genes_df$CHANGE == 1e+06)) {
            g <- g + ggplot2::scale_fill_manual(values = c(low, high), na.value = "white",
                name = legend_title)
        } else {
            g <- g + ggplot2::scale_fill_gradient2(low = low, mid = mid, high = high,
                na.value = "white", name = legend_title)
        }
    } else {
        g <- g + ggplot2::scale_fill_manual(values = c(low, high), na.value = "white",
            name = legend_title)
    }
    return(g)
}


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
UpSet_plot <- function(result_df, genes_df, num_terms = 10, method = "heatmap", use_description = FALSE,
    low = "red", mid = "black", high = "green", ...) {
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
            collapse = ", "))
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

    term_genes_mat <- matrix(0, nrow = nrow(result_df), ncol = length(all_genes),
        dimnames = list(all_terms, all_genes))
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
        decreasing = TRUE)))
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
            g <- g + ggplot2::geom_tile(data = plot_df, ggplot2::aes(x = .data$Term,
                y = .data$Gene, fill = .data$Up_Down), color = "gray60")
            g <- g + ggplot2::scale_fill_manual(values = c(low, high))
        } else {
            plot_df$Value <- genes_df$CHANGE[match(names(plot_df$Term), genes_df$GENE)]
            g <- g + ggplot2::geom_tile(data = plot_df, ggplot2::aes(x = .data$Term,
                y = .data$Gene, fill = .data$Value), color = "gray60")
            g <- g + ggplot2::scale_fill_gradient2(low = low, mid = mid, high = high)
        }
        g <- g + ggplot2::theme_minimal()
        g <- g + ggplot2::theme(axis.title = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(), legend.title = ggplot2::element_blank())
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
        "degree"), reverse = !missing(genes_df))
    return(g)
}
