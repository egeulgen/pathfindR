#' Convert an upper-tail p-value to a z-score
#'
#' Returns the value \code{z} such that \eqn{P(Z > z) = p} for a standard
#' normal variable, i.e. \code{qnorm(1 - p)}. Significant (small) p-values map
#' to large positive z-scores; \code{p = 0.5} maps to \code{0}.
#'
#' @param p Numeric vector of p-values, assumed already clamped away from 0/1.
#'
#' @return Numeric vector of z-scores.
.upper_tail_zscore <- function(p) {
  stats::qnorm(1 - p)
}

#' Parse the experiment input into a clean gene / p-value data frame
#'
#' Accepts a data frame. If columns named \code{gene} and \code{pvalue}
#' (case-insensitive) exist they are used, otherwise the first two columns are
#' taken as gene and p-value respectively. Gene names are upper-cased.
#'
#' @param experiment A data frame of gene / p-value pairs.
#'
#' @return A data frame with character column \code{gene} and numeric column
#'   \code{pvalue}.
.parse_experiment <- function(experiment) {
  if (!is.data.frame(experiment)) {
    stop("`experiment` must be a data frame of gene / p-value pairs.")
  }
  cn <- tolower(names(experiment))
  gene_col <- match("gene", cn)
  pval_col <- match("pvalue", cn)
  if (is.na(gene_col)) gene_col <- 1L
  if (is.na(pval_col)) pval_col <- 2L
  if (ncol(experiment) < 2L) {
    stop("`experiment` must have at least two columns (gene, p-value).")
  }
  data.frame(
    gene = toupper(as.character(experiment[[gene_col]])),
    pvalue = as.numeric(experiment[[pval_col]]),
    stringsAsFactors = FALSE
  )
}

#' Build the undirected interaction network from a SIF file
#'
#' Reads a Simple Interaction Format (SIF) file and converts it into an
#' undirected \code{\link[igraph]{igraph}} graph object. The SIF file is expected to
#' contain at least three columns where the first and third columns represent
#' interacting genes or nodes. The second column (interaction type) is ignored.
#'
#' @param sif_path Character string specifying the path to the SIF file.
#'   The file should be tab-delimited and contain at least 3 columns.
#'
#' @return A list with elements: \code{g} (an \code{igraph} graph),
#'   \code{nodes} (character vector of node names in graph order),
#'   \code{nbr} (named list of neighbour-name vectors) and \code{name2id}
#'   (named integer vector mapping node name to vertex id).
.build_network <- function(sif_path) {
  sif <- utils::read.delim(file = sif_path, quote = "", header = FALSE)

  if (ncol(sif) < 3) {
    stop("SIF file must contain at least 3 columns.")
  }

  edges <- sif[, c(1, 3)]
  colnames(edges) <- c("source", "target")

  g <- igraph::graph_from_data_frame(edges, directed = FALSE)
  g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)

  node_names <- igraph::V(g)$name
  if (is.null(node_names)) node_names <- character(0)

  nbr <- lapply(igraph::adjacent_vertices(g, igraph::V(g)), igraph::as_ids)
  names(nbr) <- node_names

  result <- list(
    g       = g,
    nodes   = node_names,
    nbr     = nbr,
    name2id = stats::setNames(seq_along(node_names), node_names)
  )
  return(result)
}
