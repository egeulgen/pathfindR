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
#'
#' @export
build_network <- function(sif_path) {
  sif <- utils::read.delim(file = sif_path, quote = "", header = FALSE)

  if (ncol(sif) < 3) {
    stop("SIF file must contain at least 3 columns.")
  }

  src <- toupper(as.character(sif[[1]]))
  tgt <- toupper(as.character(sif[[3]]))

  # Track first-seen insertion order for MC calibration.
  all_mentioned <- c(rbind(src, tgt))
  insertion_order <- all_mentioned[!duplicated(all_mentioned)]

  # Remove self-loops and duplicate edges before handing off to igraph.
  # Doing this here (rather than via igraph::simplify) preserves the vertex
  # insertion order that igraph::graph_from_data_frame produces, which must
  # match the SIF insertion order above.
  keep <- src != tgt
  src <- src[keep]
  tgt <- tgt[keep]
  keep2 <- !duplicated(paste0(pmin(src, tgt), "\x01", pmax(src, tgt)))
  src <- src[keep2]
  tgt <- tgt[keep2]

  g <- igraph::graph_from_data_frame(
    data.frame(source = src, target = tgt, stringsAsFactors = FALSE),
    directed = FALSE
  )

  node_names <- igraph::V(g)$name
  if (is.null(node_names)) node_names <- character(0)

  nbr <- lapply(igraph::adjacent_vertices(g, igraph::V(g)), igraph::as_ids)
  names(nbr) <- node_names

  list(
    g               = g,
    nodes           = node_names,
    nbr             = nbr,
    name2id         = stats::setNames(seq_along(node_names), node_names),
    insertion_order = insertion_order
  )
}
