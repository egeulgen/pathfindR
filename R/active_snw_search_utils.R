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

#' Build the undirected interaction network from an adjacency list
#'
#' The adjacency list is a named list mapping each node to a character vector
#' of its neighbours, e.g. \code{list(A = c("B", "C"), C = "X")}. Edges are
#' treated as undirected, self-interactions are dropped, duplicate edges are
#' collapsed, and all node names are upper-cased. Nodes that only ever appear
#' in a self-interaction are not added.
#'
#' @param adjacency A named list describing the interactions.
#'
#' @return A list with elements: \code{g} (an \code{igraph} graph),
#'   \code{nodes} (character vector of node names in graph order),
#'   \code{nbr} (named list of neighbour-name vectors) and \code{name2id}
#'   (named integer vector mapping node name to vertex id).
.build_network <- function(adjacency) {
  if (!is.list(adjacency)) {
    stop("`adjacency` must be a named list mapping each node to its neighbours.")
  }
  if (length(adjacency) > 0L && is.null(names(adjacency))) {
    stop("`adjacency` must be a named list mapping each node to its neighbours.")
  }

  src <- toupper(names(adjacency))
  tgt <- lapply(adjacency, function(x) toupper(as.character(x)))
  lens <- lengths(tgt)

  from <- rep(src, lens)
  to <- unlist(tgt, use.names = FALSE)

  keep <- if (length(from) == 0L) logical(0) else from != to
  from <- from[keep]
  to <- to[keep]

  all_nodes <- unique(c(from, to))

  if (length(all_nodes) == 0L) {
    g <- igraph::make_empty_graph(n = 0, directed = FALSE)
  } else {
    edges_df <- data.frame(from = from, to = to, stringsAsFactors = FALSE)
    vertices_df <- data.frame(name = all_nodes, stringsAsFactors = FALSE)
    g <- igraph::graph_from_data_frame(edges_df,
      directed = FALSE,
      vertices = vertices_df
    )
    g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
  }

  node_names <- igraph::V(g)$name
  if (is.null(node_names)) node_names <- character(0)

  nbr <- lapply(igraph::adjacent_vertices(g, igraph::V(g)), igraph::as_ids)
  names(nbr) <- node_names

  list(
    g       = g,
    nodes   = node_names,
    nbr     = nbr,
    name2id = stats::setNames(seq_along(node_names), node_names)
  )
}
