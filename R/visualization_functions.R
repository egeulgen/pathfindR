#' Create Diagrams for Enriched Terms
#'
#' @param result_df Data frame of enrichment results. Must-have columns for
#'  KEGG human pathway diagrams (\code{hsa_kegg = TRUE}) are: "ID" and "Term_Description".
#'  Must-have columns for the rest are: "Term_Description", "Up_regulated" and
#' "Down_regulated"
#' @param input_processed input data processed via \code{\link{input_processing}},
#'  not necessary when \code{hsa_KEGG = FALSE}
#' @param hsa_KEGG boolean to indicate whether human KEGG gene sets were used for
#'  enrichment analysis or not (default = \code{TRUE})
#' @inheritParams return_pin_path
#' @param ... additional arguments for \code{\link{visualize_hsa_KEGG}} (used
#' when \code{hsa_kegg = TRUE})
#'
#' @return Depending on the argument \code{hsa_KEGG}, creates visualization of
#'  interactions of genes involved in the list of enriched terms in
#'  \code{result_df} and saves them in the folder "term_visualizations" under
#'  the current working directory.
#'
#'
#' @details For \code{hsa_KEGG = TRUE}, KEGG human pathway diagrams are created,
#' affected nodes colored by up/down regulation status.
#' For other gene sets, interactions of affected genes are determined (via a shortest-path
#' algorithm) and are visualized (colored by change status) using igraph.
#'
#'
#' @export
#'
#' @seealso See \code{\link{visualize_hsa_KEGG}} for the visualization function
#' of human KEGG diagrams. See \code{\link{visualize_term_interactions}} for the
#' visualization function that generates diagrams showing the interactions of
#' input genes in the PIN. See \code{\link{run_pathfindR}} for the wrapper
#' function of the pathfindR workflow.
#'
#' @examples
#' \dontrun{
#' visualize_terms(result_df, input_processed)
#' visualize_terms(result_df, hsa_KEGG = FALSE, pin_name_path = "IntAct")
#' }
visualize_terms <- function(result_df, input_processed = NULL,
                            hsa_KEGG = TRUE, pin_name_path = "Biogrid",
                            ...) {

  ############ Argument Checks
  if (!is.data.frame(result_df)) {
    stop("`result_df` should be a data frame")
  }

  if (!is.logical(hsa_KEGG)) {
    stop("the argument `hsa_KEGG` should be either TRUE or FALSE")
  }

  if (hsa_KEGG) {
    nec_cols <- "ID"
  } else {
    nec_cols <- c("Term_Description", "Up_regulated", "Down_regulated")
  }
  if (!all(nec_cols %in% colnames(result_df))) {
    stop("`result_df` should contain the following columns: ",
         paste(dQuote(nec_cols), collapse = ", "))
  }

  if (hsa_KEGG) {
    if (is.null(input_processed)) {
      stop("`input_processed` should be specified when `hsa_KEGG = TRUE`")
    }
  }

  ############ Generate Diagrams
  if (hsa_KEGG) {
    visualize_hsa_KEGG(hsa_kegg_ids = result_df$ID,
                       input_processed = input_processed,
                       ...)
  } else {
    visualize_term_interactions(result_df = result_df,
                                pin_name_path = pin_name_path)
  }
}

#' Visualize Interactions of Genes Involved in the Given Enriched Terms
#'
#' @param result_df Data frame of enrichment results. Must-have columns
#' are: "Term_Description", "Up_regulated" and "Down_regulated"
#' @inheritParams return_pin_path
#'
#' @return Creates PNG files visualizing the interactions of genes involved
#' in the given enriched terms (annotated in the \code{result_df}) in the PIN used
#' for enrichment analysis (specified by \code{pin_name_path}). The PNG files are
#' saved in the folder "term_visualizations" under the current working directory.
#'
#' @details The following steps are performed for the visualization of interactions
#' of genes involved for each enriched term: \enumerate{
#'   \item shortest paths between all affected genes are determined (via \code{\link[igraph]{igraph}})
#'   \item the nodes of all shortest paths are merged
#'   \item the PIN is subsetted using the merged nodes (genes)
#'   \item using the PIN subset, the graph showing the interactions is generated
#'   \item the final graph is visualized using \code{\link[igraph]{igraph}}, colored by changed
#'   status (if provided), and is saved as a PNG file.
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
#' visualize_term_interactions(result_df, pin_name_path = "IntAct")
#' }
visualize_term_interactions <- function(result_df, pin_name_path) {
  ############ Initial Steps
  ## fix naming issue
  result_df$Term_Description <- gsub("\\/", "-", result_df$Term_Description)

  ## load PIN
  pin_path <- return_pin_path(pin_name_path)
  pin <- utils::read.delim(file = pin_path,
                           header = FALSE, stringsAsFactors = FALSE)
  pin$V2 <- NULL

  pin[, 1] <- base::toupper(pin[, 1])
  pin[, 2] <- base::toupper(pin[, 2])

  ## pin graph
  pin_g <- igraph::graph_from_data_frame(pin, directed = FALSE)

  ## Create visualization output directory
  dir.create("term_visualizations", showWarnings = FALSE)

  ############ Visualize interactions by enriched term
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
      message(paste0("< 2 genes, skipping visualization of ",
                     current_row$Term_Description))
    } else {
      cat("Visualizing:", current_row$Term_Description,
          paste(rep(" ", 200), collapse = ""), "\r")

      ## Find genes without direct interaction
      cond1 <- pin$V1 %in% current_genes
      cond2 <- pin$V3 %in% current_genes
      direct_interactions <- pin[cond1 & cond2, ]
      tmp <- c(direct_interactions$V1, direct_interactions$V3)
      missing_genes <- current_genes[!current_genes %in% tmp]

      ## Find shortest path between genes without direct interaction
      # and other current_genes
      s_path_genes <- c()
      for (gene in missing_genes) {
        tmp <- suppressWarnings(igraph::shortest_paths(pin_g,
                                                       from = which(names(igraph::V(pin_g)) == gene),
                                                       to = which(names(igraph::V(pin_g)) %in% current_genes),
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
      igraph::V(g)$color <- ifelse(cond1, "red",
                                   ifelse(cond2, "green",
                                          ifelse(cond3, "blue", "gray60")))

      path_to_png <- file.path("term_visualizations",
                               paste0(current_row$Term_Description, ".png"))

      #### Generate diagram
      term_diagram <- magick::image_graph(width = 1200, height = 900, res = 100)
      # Plot the tree object
      igraph::plot.igraph(g,
                          layout = igraph::layout.fruchterman.reingold,
                          edge.curved = FALSE,
                          vertex.size = 10,
                          vertex.label.dist = 0,
                          vertex.label.color = "black",
                          asp = FALSE,
                          vertex.label.cex = 0.8,
                          edge.width = 1.2,
                          edge.arrow.mode = 0,
                          main = paste(current_row$Term_Description,
                                       "\n Involved Gene Interactions in",
                                       pin_name_path))

      if (is.null(snw_genes)) {
        graphics::legend("topleft",
                         legend = c("Upregulated Input Genes",
                                    "Downregulated Input Genes",
                                    "Other"),
                         col = c("red", "green", "gray60"),
                         pch = 19, cex = 1.5, bty = "n")
      } else {
        graphics::legend("topleft",
                         legend = c("Non-input Active Snw. Genes",
                                    "Upregulated Input Genes",
                                    "Downregulated Input Genes",
                                    "Other"),
                         col = c("blue", "red", "green", "gray60"),
                         pch = 19, cex = 1.5, bty = "n")
      }
      grDevices::dev.off()

      #### Add logo to diagram
      ### Read logo
      path_logo <- system.file("extdata", "logo.png", package = "pathfindR")
      logo_img <- magick::image_read(path_logo)

      ### Add logo
      term_diagram <- magick::image_composite(term_diagram,
                                              magick::image_scale(logo_img, "x100"),
                                              gravity = "northeast",
                                              offset = "+10+10")

      #### Save file
      magick::image_write(term_diagram, path = path_to_png, format = "png")
    }
  }
}

#' Visualize Human KEGG Pathways
#'
#' @param hsa_kegg_ids hsa KEGG ids of pathways to be colored and visualized
#' @param input_processed input data processed via \code{\link{input_processing}}
#' @inheritParams color_kegg_pathway
#' @param key_gravity gravity value (character) for the color key legend placement
#' (see \code{\link[magick]{gravity_types}})
#' @param logo_gravity gravity value (character) for the logo placement
#' (see \code{\link[magick]{gravity_types}})
#'
#' @return Creates colored visualizations of the enriched human KEGG pathways
#' and saves them in the folder "term_visualizations" under the current working
#' directory.
#'
#' @seealso See \code{\link{visualize_terms}} for the wrapper function for
#' creating enriched term diagrams. See \code{\link{run_pathfindR}} for the
#' wrapper function of the pathfindR enrichment workflow.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' visualize_hsa_KEGG(hsa_kegg_ids, input_processed)
#' }
visualize_hsa_KEGG <- function(hsa_kegg_ids, input_processed,
                               normalize_vals = TRUE, node_cols = NULL,
                               quiet = TRUE,
                               key_gravity = "northeast",
                               logo_gravity = "southeast") {
  ############ Arg checks

  ### hsa_kegg_ids
  if (!is.atomic(hsa_kegg_ids)) {
    stop("`hsa_kegg_ids` should be a vector of hsa KEGG IDs")
  }
  if (!all(grepl("^hsa[0-9]{5}$", hsa_kegg_ids))) {
    stop("`hsa_kegg_ids` should be a vector of valid hsa KEGG IDs")
  }

  ### input_processed
  if (!is.data.frame(input_processed)) {
    stop("`input_processed` should be a data frame")
  }

  nec_cols <- c("GENE", "CHANGE")
  if (!all(nec_cols %in% colnames(input_processed))) {
    stop("`input_processed` should contain the following columns: ",
         paste(dQuote(nec_cols), collapse = ", "))
  }

  ############ Create change vector
  ### Convert gene symbols into NCBI gene IDs
  tmp <- AnnotationDbi::mget(input_processed$GENE,
                             AnnotationDbi::revmap(org.Hs.eg.db::org.Hs.egSYMBOL),
                             ifnotfound = NA)
  input_processed$EG_ID <- unlist(tmp)

  ### A rule of thumb for the 'kegg' ID is entrezgene ID for eukaryote species
  input_processed$KEGG_ID  <- paste0("hsa:", input_processed$EG_ID)

  ############ Fetch all pathway genes, create vector of change values and
  ############ Generate colored pathway diagrams for each pathway
  change_vec <- input_processed$CHANGE
  names(change_vec) <- input_processed$KEGG_ID

  cat("Downloading pathway diagrams of", length(hsa_kegg_ids), "KEGG pathways\n\n")
  pw_vis_list <- list()
  pb <- utils::txtProgressBar(min = 0, max = length(hsa_kegg_ids), style = 3)
  for (i in seq_len(length(hsa_kegg_ids))) {
    pw_id <- hsa_kegg_ids[i]

    pw_vis_list[[pw_id]] <- color_kegg_pathway(pw_id = pw_id,
                                               change_vec = change_vec,
                                               normalize_vals = normalize_vals,
                                               node_cols = node_cols,
                                               quiet = quiet)

    utils::setTxtProgressBar(pb, i)
  }
  close(pb)

  ############ Add logo and color key legend per each diagram

  ### Read logo
  path_logo <- system.file("extdata", "logo.png", package = "pathfindR")
  logo_img <- magick::image_read(path_logo)

  dir.create("term_visualizations", showWarnings = FALSE)

  cat("Saving colored pathway diagrams of", length(pw_vis_list), "KEGG pathways\n\n")
  pb <- utils::txtProgressBar(min = 0, max = length(pw_vis_list), style = 3)
  for (i in seq_len(length(pw_vis_list))) {
    ### Read image
    f_path <- pw_vis_list[[i]]$file_path
    pw_diag <- magick::image_read(f_path)

    ### Add logo
    pw_diag <- magick::image_composite(pw_diag,
                                       magick::image_scale(logo_img, "x90"),
                                       gravity = logo_gravity,
                                       offset = "+10+10")

    ### Prep for color keys
    key_col_df <- data.frame(bin_val = seq_along(pw_vis_list[[i]]$all_key_cols),
                             color = pw_vis_list[[i]]$all_key_cols,
                             y_val = 1)

    key_breaks <- pw_vis_list[[i]]$all_brks
    names(key_breaks) <- seq_along(key_breaks)
    key_breaks <- c(key_breaks[1], mean(key_breaks[1:6]), key_breaks[6], mean(key_breaks[6:11]), key_breaks[11])
    brks <- c(.5, 3, 5.5, 8, 10.5)

    ### Generate color legend image
    col_key_legend <- magick::image_graph(width = 200, height = 90, res = 100)
    g <- ggplot2::ggplot(key_col_df, ggplot2::aes_(~bin_val, ~y_val))
    g <- g + ggplot2::geom_tile(fill = key_col_df$color,
                                colour = "black")
    g <- g + ggplot2::scale_x_continuous(expand = c(0, 0),
                                         breaks = brks,
                                         labels = base::format(key_breaks,
                                                               digits = 2))
    g <- g + ggplot2::scale_y_discrete(expand = c(0, 0))
    g <- g + ggplot2::theme_bw()
    g <- g + ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                            panel.grid.minor = ggplot2::element_blank(),
                            axis.title.x = ggplot2::element_blank(),
                            axis.title.y = ggplot2::element_blank(),
                            axis.ticks = ggplot2::element_line(colour = "black",
                                                               size = .6),
                            axis.ticks.length = ggplot2::unit(.2, "cm"),
                            axis.text.x = ggplot2::element_text(size = 14,
                                                                face = "bold"),
                            panel.border = ggplot2::element_rect(colour = "black",
                                                                 fill = NA,
                                                                 size = .5),
                            plot.margin = ggplot2::unit(c(0, .7 , 0, .7), "cm"))
    print(g)
    grDevices::dev.off()

    ### Add color key legend
    if (all(change_vec == 1e6)) {
      final_pw_img <- pw_diag
    } else {
      final_pw_img <- magick::image_composite(pw_diag,
                                              magick::image_scale(col_key_legend, "x45"),
                                              gravity  = key_gravity,
                                              offset = "+10+10")
    }

    final_path <- file.path("term_visualizations", basename(f_path))
    magick::image_write(final_pw_img, path = final_path, format = "png")

    utils::setTxtProgressBar(pb, i)
  }
  close(pb)
}

#' Color hsa KEGG pathway
#'
#' @param pw_id hsa KEGG pathway id (e.g. hsa05012)
#' @param change_vec vector of change values, names should be hsa KEGG gene ids
#' @param normalize_vals should change values be normalized (default = \code{TRUE})
#' @param node_cols low, middle and high color values for coloring the pathway nodes
#' (default = \code{NULL}). If \code{node_cols=NULL}, the low, middle and high color
#' are set as "green", "gray" and "red". If all change values are 1e6 (in case no
#' changes are supplied, this dummy value is assigned by
#' \code{\link{input_processing}}), only one color ("#F38F18" if NULL) is used.
#' @param quiet If \code{TRUE}, suppress status messages (if any), and the
#' progress bar while downloading file(s)
#'
#' @return list containing: \enumerate{
#'    \item file_path: path to colored hsa KEGG pathway diagram
#'    \item all_key_cols: colors used for each change value bin
#'    \item all_brks: breaks used for separating change values into bins
#' }
#'
#' @examples
#' \dontrun{
#'    pw_id <- "hsa00010"
#'    change_vec <- c(-2, 4, 6)
#'    names(change_vec) <- c("hsa:2821", "hsa:226", "hsa:229")
#'    result <- pathfindR:::color_kegg_pathway(pw_id, change_vec)
#' }
color_kegg_pathway <- function(pw_id, change_vec, normalize_vals = TRUE,
                               node_cols = NULL, quiet = TRUE) {
  ############ Arg checks
  if (!is.logical(normalize_vals)) {
    stop("`normalize_vals` should be logical")
  }

  ## check node_cols
  if (!is.null(node_cols)) {
    if (!is.atomic(node_cols)) {
      stop("`node_cols` should be a vector of colors")
    }

    if (!all(change_vec == 1e6) & length(node_cols) != 3) {
      stop("the length of `node_cols` should be 3")
    }

    isColor <- function(x) {
      tryCatch(is.matrix(grDevices::col2rgb(x)),
               error = function(e) FALSE)
    }

    if (!all(vapply(node_cols, isColor, TRUE))) {
      stop("`node_cols` should be a vector of valid colors")
    }
  }
  ############ Set node palette
  ### if node_cols not supplied, use default color(s)
  if (!is.null(node_cols)) {
    if (all(change_vec == 1e6)) {
      message("all `change_vec` values are 1e6, using the first color in `node_cols`")
      low_col <- mid_col <- high_col <- node_cols[1]
    } else {
      low_col <- node_cols[1]
      mid_col <- node_cols[2]
      high_col <- node_cols[3]
    }
  } else if (all(change_vec == 1e6)) { ## NO CHANGES SUPPLIED
    low_col <- mid_col <- high_col <- "#F38F18"
  } else {
    low_col <- "green"
    mid_col <- "gray"
    high_col <- "red"
  }

  ############ Summarization for genes that are in the same node
  ############ and handling of non-input pathway genes
  ## download KGML to determine gene nodes
  pwKGML <- tempfile()

  tryCatch({status <- KEGGgraph::retrieveKGML(sub("hsa", "", pw_id),
                                              organism = "hsa",
                                              destfile = pwKGML,
                                              quiet = quiet)},
           warning = function(w) {
             message(paste("Warning for :", pw_id))
             message("Here's the original warning message:")
             message(w)
             return(NULL)
           }, error = function(e) {
             message(paste("Cannot reach URL, please check your connection:", pw_id))
             message("Here's the original error message:")
             message(e)
             return(NA)
           }, finally = {
             invisible(status)
           })

  current_pw <- KEGGgraph::parseKGML(pwKGML)

  all_nodes <- KEGGgraph::nodes(current_pw)
  gene_nodes <- list()
  for (node in all_nodes) {
    if (KEGGgraph::getType(node) == "gene") {
      gene_nodes[[KEGGgraph::getEntryID(node)]] <- KEGGgraph::getName(node)
    }
  }
  gene_nodes <- unique(gene_nodes)
  gene_nodes <- gene_nodes[order(vapply(gene_nodes, length, 1L))]

  ## summarize over all pathway gene nodes
  ## keeping one unique gene per all nodes
  input_genes <- names(change_vec)
  pw_vis_changes <- c()
  for (i in seq_len(length(gene_nodes))) {
    node <- gene_nodes[[i]]
    cond <- input_genes %in% node

    tmp_val <- mean(change_vec[input_genes[cond]])

    if (length(node) != 1) {
      rest_nodes <- unlist(gene_nodes[-i])
      chosen_nm <- setdiff(node, rest_nodes)
    } else {
      chosen_nm <- node
    }

    names(tmp_val) <- chosen_nm[!chosen_nm %in% names(pw_vis_changes)][1]
    pw_vis_changes <- c(pw_vis_changes, tmp_val)
  }

  ############ Determine node colors
  vals <- pw_vis_changes[!is.na(pw_vis_changes)]
  ### Normalization
  if (!all(vals == 1e6) & normalize_vals) {
    vals <- 2 * (vals - min(vals)) / diff(range(vals)) - 1
  }

  ### determine limit
  lim <- round(max(abs(vals)), 2)

  ### generate low colors
  low_vals <- vals[vals < 0]
  low_brks <- seq(from = -lim, to = 0, length.out = 6)
  low_bins <- cut(low_vals, breaks = low_brks,
                  ordered_result = TRUE, include.lowest = TRUE)
  key_col_fn <- grDevices::colorRampPalette(c(low_col, mid_col))
  low_key_cols <- key_col_fn(5)
  names(low_key_cols) <- levels(low_bins)
  names(low_bins) <- names(low_vals)

  ### generate high colors
  high_vals <- vals[vals > 0]
  high_brks <- seq(from = 0, to = lim, length.out = 6)
  high_bins <- cut(high_vals, breaks = high_brks,
                   ordered_result = TRUE, include.lowest = TRUE)
  key_col_fn <- grDevices::colorRampPalette(c(mid_col, high_col))
  high_key_cols <- key_col_fn(5)
  names(high_key_cols) <- levels(high_bins)
  names(high_bins) <- names(high_vals)

  all_brks <- c(low_brks, high_brks[-1])
  all_key_cols <- c(low_key_cols, high_key_cols)

  ############ Label each pw gene with the appropriate color
  fg_cols <- ifelse(names(pw_vis_changes) %in% names(low_bins),
                    all_key_cols[low_bins[names(pw_vis_changes)]],
                    ifelse(names(pw_vis_changes) %in% names(high_bins),
                           all_key_cols[high_bins[names(pw_vis_changes)]],
                           "white"))

  bg_cols <- rep("black", length(pw_vis_changes))

  ############ Download colored KEGG pathway diagram
  fname <- paste0(pw_id, "_pathfindR.png")
  f_path <- file.path(tempdir(check = TRUE), fname)
  pw_url <- KEGGREST::color.pathway.by.objects(pw_id, names(pw_vis_changes),
                                               fg.color.list = fg_cols,
                                               bg.color.list = bg_cols)

  tryCatch({status <- utils::download.file(url = pw_url,
                                           destfile = f_path,
                                           quiet = quiet)},
           warning = function(w) {
             message(paste("URL caused a warning:", url))
             message("Here's the original warning message:")
             message(w)
             return(NULL)
           }, error = function(e) {
             message(paste("Cannot reach URL, please check your connection:", url))
             message("Here's the original error message:")
             message(e)
             return(NA)
           }, finally = {
             invisible(status)
           })

  return(list(file_path = f_path,
              all_key_cols = all_key_cols,
              all_brks = all_brks))
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
#' @param top_terms number of top terms (according to the "lowest_p" column)
#'  to plot (default = 10). If \code{plot_by_cluster = TRUE}, selects the top
#'  \code{top_terms} terms per each cluster. Set \code{top_terms = NULL} to plot
#'  for all terms.If the total number of terms is less than \code{top_terms},
#'  all terms are plotted.
#' @param plot_by_cluster boolean value indicating whether or not to group the
#'  enriched terms by cluster (works if \code{result_df} contains a
#'  "Cluster" column).
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
#' Optionally, if "Cluster" is a column of \code{result_df} and
#' \code{plot_by_cluster == TRUE}, the enriched terms are grouped by clusters.
#'
#' @import ggplot2
#' @export
#'
#' @examples
#' g <- enrichment_chart(RA_output)
enrichment_chart <- function(result_df,
                             top_terms = 10,
                             plot_by_cluster = FALSE,
                             num_bubbles = 4,
                             even_breaks = TRUE) {
  necessary <- c("Term_Description", "Fold_Enrichment", "lowest_p",
                 "Up_regulated", "Down_regulated")

  if (!all(necessary %in% colnames(result_df))) {
    stop("The input data frame must have the columns:\n",
         paste(necessary, collapse = ", "))
  }

  if (!is.logical(plot_by_cluster)) {
    stop("`plot_by_cluster` must be either TRUE or FALSE")
  }

  if (!is.numeric(top_terms) & !is.null(top_terms)) {
    stop("`top_terms` must be either numeric or NULL")
  }

  if (!is.null(top_terms)){
    if (top_terms < 1) {
      stop("`top_terms` must be > 1")
    }
  }

  # sort by lowest adj.p
  result_df <- result_df[order(result_df$lowest_p), ]

  ## Filter for top_terms
  if (!is.null(top_terms)) {
    if (plot_by_cluster & "Cluster" %in% colnames(result_df)) {
      keep_ids <- tapply(result_df$ID, result_df$Cluster, function(x)
                                          x[seq_len(min(top_terms, length(x)))])
      keep_ids <- unlist(keep_ids)
      result_df <- result_df[result_df$ID %in% keep_ids, ]
    } else if (top_terms < nrow(result_df)){
      result_df <- result_df[seq_len(top_terms), ]
    }
  }

  num_genes <- vapply(result_df$Up_regulated,
    function(x) length(unlist(strsplit(x, ", "))), 1)
  num_genes <- num_genes + vapply(result_df$Down_regulated,
    function(x) length(unlist(strsplit(x, ", "))), 1)

  result_df$Term_Description <- factor(result_df$Term_Description,
          levels = rev(unique(result_df$Term_Description)))

  g <- ggplot2::ggplot(result_df, ggplot2::aes_(~Fold_Enrichment, ~Term_Description))
  g <- g + ggplot2::geom_point(ggplot2::aes(color = -log10(result_df$lowest_p),
                                            size = num_genes), na.rm = TRUE)
  g <- g + ggplot2::theme_bw()
  g <- g + ggplot2::theme(axis.text.x = ggplot2::element_text(size = 10),
                          axis.text.y = ggplot2::element_text(size = 10),
                          plot.title = ggplot2::element_blank())
  g <- g + ggplot2::xlab("Fold Enrichment")
  g <- g + ggplot2::theme(axis.title.y = ggplot2::element_blank())
  g <- g + ggplot2::labs(size = "# genes", color = expression(-log[10](p)))

  ## breaks for # genes
  if (max(num_genes) < num_bubbles) {
    g <- g + ggplot2::scale_size_continuous(breaks = seq(0, max(num_genes)))
  } else {

    if (even_breaks) {
      brks <- base::seq(0, max(num_genes),
                        round(max(num_genes) / (num_bubbles + 1)))
    } else {
      brks <- base::round(base::seq(0, max(num_genes),
                                    length.out = num_bubbles + 1))
    }
    g <- g + ggplot2::scale_size_continuous(breaks = brks)
  }

  g <- g + ggplot2::scale_color_continuous(low = "#f5efef", high = "red")

  if (plot_by_cluster & "Cluster" %in% colnames(result_df)) {
    g <- g + ggplot2::facet_grid(result_df$Cluster ~ .,
      scales = "free_y", space = "free", drop = TRUE
    )
  } else if (plot_by_cluster) {
    message("For plotting by cluster, there must a column named `Cluster` in the input data frame!")
  }

  return(g)
}


#' Plot Term-Gene Graph
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
#' @param layout The type of layout to create (see \code{\link[ggraph]{ggraph}} for details. Default = "auto")
#' @param use_description Boolean argument to indicate whether term descriptions
#'  (in the "Term_Description" column) should be used. (default = \code{FALSE})
#' @param node_size Argument to indicate whether to use number of significant genes ("num_genes")
#'  or the -log10(lowest p value) ("p_val") for adjusting the node sizes (default = "num_genes")
#'
#' @return a  \code{\link[ggraph]{ggraph}} object containing the term-gene graph.
#'  Each node corresponds to an enriched term (beige), an up-regulated gene (green)
#'  or a down-regulated gene (red). An edge between a term and a gene indicates
#'  that the given term involves the gene. Size of a term node is proportional
#'  to either the number of genes (if \code{node_size = "num_genes"}) or
#'  the -log10(lowest p value) (if \code{node_size = "p_val"}).
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
#' p <- term_gene_graph(RA_output)
#' p <- term_gene_graph(RA_output, num_terms = 5)
#' p <- term_gene_graph(RA_output, node_size = "p_val")
term_gene_graph <- function(result_df, num_terms = 10,
                            layout = "auto", use_description = FALSE,
                            node_size = "num_genes") {

  ############ Argument Checks
  ### Check num_terms is NULL or numeric
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
  ### Check necessary columnns
  necessary_cols <- c("Up_regulated", "Down_regulated", "lowest_p", ID_column)

  if (!all(necessary_cols %in% colnames(result_df))) {
    stop(paste(c("All of", paste(necessary_cols, collapse = ", "),
                 "must be present in `results_df`!"), collapse = " "))
  }

  ############ Initial steps
  ### set num_terms to NULL if number of enriched terms is smaller than num_terms
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
    genes <- c(up_genes, down_genes)

    for (gene in genes) {
      graph_df <- rbind(graph_df,
                         data.frame(Term = result_df[i, ID_column],
                                    Gene = gene))
    }
  }


  up_genes <- lapply(result_df$Up_regulated,
                     function(x) unlist(strsplit(x, ", ")))
  up_genes <- unlist(up_genes)

  ############ Create graph object and plot
  ### create igraph object
  g <- igraph::graph_from_data_frame(graph_df, directed = FALSE)
  cond_term <- names(igraph::V(g)) %in% result_df[, ID_column]
  cond_up_gene <- names(igraph::V(g)) %in% up_genes
  igraph::V(g)$type <- ifelse(cond_term, "term",
                              ifelse(cond_up_gene, "up", "down"))
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
  igraph::V(g)$color <- ifelse(igraph::V(g)$type == "term", "#E5D7BF",
                               ifelse(igraph::V(g)$type == "up", "green",
                                      "red"))

  ### Create graph
  p <- ggraph::ggraph(g, layout = layout)
  p <- p + ggraph::geom_edge_link(alpha = .8, colour = "darkgrey")
  p <- p + ggraph::geom_node_point(ggplot2::aes_(color = ~ I(color), size = ~size))
  p <- p + ggplot2::scale_size(range = c(5, 10),
                               breaks = round(seq(round(min(igraph::V(g)$size)),
                                                  round(max(igraph::V(g)$size)),
                                                  length.out = 4)),
                               name = size_label)
  p <- p + ggplot2::theme_void()
  p <- p + ggraph::geom_node_text(ggplot2::aes_(label = ~name), nudge_y = .2)
  p <- p + ggplot2::scale_colour_manual(values = unique(igraph::V(g)$color),
                                        name = NULL,
                                        labels = c("enriched term",
                                                   "up-regulated gene",
                                                   "down-regulated gene"))
  if (is.null(num_terms)) {
    p <- p + ggplot2::ggtitle("Term-Gene Graph")
  } else {
    p <- p + ggplot2::ggtitle("Term-Gene Graph",
                              subtitle = paste(c("Top", num_terms, "terms"),
                                               collapse = " "))
  }

  p <- p + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                          plot.subtitle = ggplot2::element_text(hjust = 0.5))

  return(p)
}
