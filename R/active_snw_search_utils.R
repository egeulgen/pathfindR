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
#' undirected \code{\link[igraph]{igraph}} graph object. Following the Java
#' reference's SIFReader, the column count is taken from the first line: a 2-column
#' file uses columns 1 and 2 as the interacting nodes, a 3-column file uses columns
#' 1 and 3 (the middle interaction-type column is ignored). Node names are
#' upper-cased and self-interactions are discarded.
#'
#' To reproduce the Java implementation's greedy search bit-for-bit, this function
#' also reconstructs two order-sensitive structures that the Java code derives from
#' its HashMap/HashSet traversal:
#'   * \code{nodes}: the node order of Java's \code{networkNodeList}
#'     (\code{adjacency.keySet()} iteration order), via \code{java_node_order()}.
#'   * \code{nbr_idx}: per-node neighbour lists in Java's \code{HashSet} iteration
#'     order, via \code{java_neighbour_order()}.
#' These orders drive both the Monte-Carlo calibration (which shuffles z-scores in
#' node order) and the greedy expansion/removal, so matching them is what makes the
#' R/C++ output align with the Java reference. The \code{igraph} object is retained
#' for the SA / GA algorithms, whose component scoring is order-independent.
#'
#' @param sif_path Character string specifying the path to the SIF file.
#'   The file should be whitespace/tab-delimited with 2 or 3 columns.
#'
#' @return A list with elements: \code{g} (an \code{igraph} graph),
#'   \code{nodes} (node names in Java \code{networkNodeList} order),
#'   \code{nbr} (named list of neighbour-name vectors in Java HashSet order),
#'   \code{nbr_idx} (list of 1-based neighbour-id vectors aligned to \code{nodes},
#'   in Java HashSet order, ready for \code{run_greedy_search()}) and
#'   \code{name2id} (named integer vector mapping node name to its index in
#'   \code{nodes}) and \code{csr_offsets} / \code{csr_nbrs} (a compressed
#'   sparse-row, 0-based adjacency used by the SA / GA component scorer).
#'
#' @export
build_network <- function(sif_path) {
  # --- Read raw lines and split exactly like Java's SIFReader -----------------
  # Java: line.split("[ \\t]"); column count fixed from the first non-empty line;
  # lines with a different column count are discarded.
  raw_lines <- readLines(sif_path, warn = FALSE)
  raw_lines <- raw_lines[nzchar(raw_lines)]
  if (length(raw_lines) == 0L) {
    stop("SIF file is empty.")
  }

  split_line <- function(ln) strsplit(ln, "[ \t]")[[1]]

  first <- split_line(raw_lines[1])
  column_number <- length(first)
  if (!column_number %in% c(2L, 3L)) {
    stop("SIF file must have 2 or 3 columns.")
  }

  parts <- strsplit(raw_lines, "[ \t]")
  ok <- vapply(parts, length, integer(1)) == column_number
  parts <- parts[ok]

  mat <- do.call(rbind, parts)
  src <- toupper(mat[, 1])
  tgt <- toupper(mat[, if (column_number == 3L) 3L else 2L])

  # Discard self-interactions (Java does the same).
  keep <- src != tgt
  src <- src[keep]
  tgt <- tgt[keep]

  if (length(src) == 0L) {
    stop("SIF file contains no non-self interactions.")
  }

  # --- Java node order (networkNodeList = adjacency.keySet() iteration) --------
  # java_node_order() takes the raw, upper-cased, in-file-order endpoint columns
  # (NOT de-duplicated) and returns node names in Java HashMap key order.
  nodes <- as.character(java_node_order(src, tgt))
  N <- length(nodes)
  name2id <- stats::setNames(seq_along(nodes), nodes)

  # --- Per-node neighbour insertion order -------------------------------------
  # For edge (n1,n2) on line L, Java appends n2 to n1's neighbour set and n1 to
  # n2's, both "at time L". A node's neighbour insertion order is therefore the
  # other endpoint across the lines in which it appears, in line order, de-duped
  # keeping the first occurrence. We build that with a stable order on (node,
  # line) and a duplicate drop.
  ne <- length(src)
  from <- c(src, tgt) # the node whose neighbour list we extend
  to <- c(tgt, src) # the neighbour being added
  line <- c(seq_len(ne), seq_len(ne)) # originating line index

  o <- order(from, line, method = "radix")
  from <- from[o]
  to <- to[o]
  dup <- duplicated(paste0(from, "\x01", to))
  from <- from[!dup]
  to <- to[!dup]

  # Per-node neighbour names in insertion order (named list aligned to `nodes`).
  nbr_ins <- split(to, factor(from, levels = nodes))

  # --- Reorder each neighbour list into Java HashSet iteration order ----------
  # and convert to 1-based ids into `nodes` for the greedy search.
  nbr <- vector("list", N)
  nbr_idx <- vector("list", N)
  names(nbr) <- nodes
  for (i in seq_len(N)) {
    ins_names <- nbr_ins[[i]]
    if (length(ins_names) == 0L) {
      nbr[[i]] <- character(0)
      nbr_idx[[i]] <- integer(0)
      next
    }
    ins_ids <- name2id[ins_names]
    ord_ids <- java_neighbour_order(as.integer(ins_ids), ins_names) # Java HashSet order
    nbr_idx[[i]] <- as.integer(ord_ids)
    nbr[[i]] <- nodes[ord_ids]
  }

  keep2 <- !duplicated(paste0(pmin(src, tgt), "\x01", pmax(src, tgt)))
  g <- igraph::graph_from_data_frame(
    data.frame(source = src[keep2], target = tgt[keep2], stringsAsFactors = FALSE),
    directed = FALSE,
    vertices = data.frame(name = nodes, stringsAsFactors = FALSE)
  )

  # --- CSR adjacency for the SA / GA component scorer -------------------------
  # Flatten the (Java-ordered) per-node neighbour id lists into compressed
  # sparse row form with 0-based ids, so component_scores_sorted() can scan the
  # graph without any per-call allocation. Built once here and reused across all
  # SA iterations / GA individual evaluations.
  csr_offsets <- integer(N + 1L)
  csr_offsets[1L] <- 0L
  for (i in seq_len(N)) {
    csr_offsets[i + 1L] <- csr_offsets[i] + length(nbr_idx[[i]])
  }
  csr_nbrs <- integer(csr_offsets[N + 1L])
  if (length(csr_nbrs) > 0L) {
    # nbr_idx holds 1-based ids; the C++ scorer expects 0-based.
    csr_nbrs <- as.integer(unlist(nbr_idx, use.names = FALSE)) - 1L
  }

  list(
    g           = g,
    nodes       = nodes,
    nbr         = nbr,
    nbr_idx     = nbr_idx,
    name2id     = name2id,
    csr_offsets = csr_offsets,
    csr_nbrs    = csr_nbrs
  )
}
