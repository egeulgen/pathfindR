#' Build the score context
#'
#' Computes per-node z-scores from the experiment p-values and the
#' Monte-Carlo mean / standard deviation of the raw subnetwork score at every
#' possible subnetwork size, which are used to calibrate (normalise) scores.
#'
#' P-values are clamped to a safe range, the smallest p-value is kept when a
#' gene appears more than once, and network nodes without a p-value are given
#' \code{params$p_for_nonsignificant}.
#'
#' @param network A network as returned by \code{.build_network()}.
#' @param experiment A data frame of gene / p-value pairs.
#' @param params A list of run parameters.
#'
#' @return A list with elements \code{z} (named z-score vector),
#'   \code{means} and \code{stds} (numeric vectors indexed by subnetwork size),
#'   and \code{nodes}.
.build_score_context <- function(network, experiment, params) {
  set.seed(params$seed)
  nodes <- network$nodes
  N <- length(nodes)

  MIN_SIG <- 1e-13
  MAX_SIG <- 1 - MIN_SIG

  exp_df <- .parse_experiment(experiment)

  pmap <- stats::setNames(rep(NA_real_, N), nodes)
  in_net <- exp_df$gene %in% nodes
  genes <- exp_df$gene[in_net]
  pvals <- exp_df$pvalue[in_net]
  pvals <- pmin(pmax(pvals, MIN_SIG), MAX_SIG)

  if (length(genes) > 0L) {
    agg <- tapply(pvals, genes, min)
    pmap[names(agg)] <- as.numeric(agg)
  }
  pmap[is.na(pmap)] <- params$p_for_nonsignificant

  z <- .upper_tail_zscore(pmap)
  names(z) <- nodes

  # --- Monte-Carlo calibration (vectorised) ---------------------------------
  # Instead of a pure-R loop over `trials`, we shuffle the full z-vector in a
  # matrix (rows = trials) and compute prefix sums in one vectorised pass.
  means <- numeric(N)
  stds <- numeric(N)

  if (N > 0L) {
    set.seed(params$seed)
    trials <- 2000L
    zvec <- as.numeric(z)
    sqrtk <- sqrt(seq_len(N))

    # Build a (trials x N) matrix where each row is a random permutation
    idx_mat <- matrix(0L, nrow = trials, ncol = N)
    for (t in seq_len(trials)) {
      idx_mat[t, ] <- sample.int(N)
    }
    # z-values for all trials at once
    z_mat <- matrix(zvec[idx_mat], nrow = trials, ncol = N)
    # Cumulative sums across columns (within each trial/row)
    cs_mat <- t(apply(z_mat, 1L, cumsum))
    # Divide each column k by sqrt(k)
    score_mat <- sweep(cs_mat, 2L, sqrtk, "/")
    # Single-node subnetworks score 0
    score_mat[, 1L] <- 0

    means <- colMeans(score_mat)
    stds <- sqrt(colMeans(score_mat^2) - means^2 + 1e-7)
  }

  list(
    z     = z,
    means = means,
    stds  = stds,
    nodes = nodes
  )
}

#' Score a subnetwork from its size and z-score sum
#'
#' Single-node subnetworks always score 0. Otherwise the raw score is
#' \code{zsum / sqrt(n)}, optionally calibrated against the Monte-Carlo
#' distribution and optionally penalised for size.
#'
#' @param sc A score context from \code{.build_score_context()}.
#' @param n Number of nodes in the subnetwork.
#' @param zsum Sum of the z-scores of the subnetwork nodes.
#' @param normalize Logical; whether to calibrate the score.
#'
#' @return A numeric score.
.score_subnetwork <- function(sc, n, zsum, normalize) {
  if (n == 1L) {
    return(0)
  }
  score <- zsum / sqrt(n)
  if (isTRUE(normalize)) {
    score <- (score - sc$means[n]) / sc$stds[n]
  }
  score
}

#' Create a subnetwork object from a vector of node names
#'
#' @param sc A score context.
#' @param nodes Character vector of node names.
#'
#' @return A list with elements \code{nodes}, \code{zsum} and \code{score}.
.make_subnetwork <- function(sc, nodes) {
  if (length(nodes) == 0L) {
    return(list(nodes = character(0), zsum = 0, score = 0))
  }
  zsum <- sum(sc$z[nodes])
  list(
    nodes = nodes,
    zsum  = zsum,
    score = .score_subnetwork(sc, length(nodes), zsum, TRUE)
  )
}

#' Find the connected components among a set of "on" nodes
#'
#' Returns one group of node names per connected component of the subgraph
#' induced by \code{on_names} (isolated nodes form singleton components).
#'
#' @param network A network from \code{.build_network()}.
#' @param on_names Character vector of node names that are switched on.
#'
#' @return A list of character vectors, one per component.
.find_components_named <- function(network, on_names) {
  if (length(on_names) == 0L) {
    return(list())
  }
  sub <- igraph::induced_subgraph(network$g, on_names)
  comp <- igraph::components(sub)
  unname(split(names(comp$membership), comp$membership))
}

#' Find the scored subnetworks among a set of "on" nodes
#'
#' @param network A network from \code{.build_network()}.
#' @param sc A score context.
#' @param on_names Character vector of node names that are switched on.
#'
#' @return A list of subnetwork objects.
.find_subnetworks <- function(network, sc, on_names) {
  groups <- .find_components_named(network, on_names)
  lapply(groups, function(ns) .make_subnetwork(sc, ns))
}

#' Sort a list of subnetwork objects by score, descending
#'
#' @param subs A list of subnetwork objects.
#'
#' @return The list reordered so the highest scoring subnetwork is first.
.sort_subnetworks_desc <- function(subs) {
  if (length(subs) <= 1L) {
    return(subs)
  }
  scores <- vapply(subs, function(s) s$score, numeric(1))
  subs[order(scores, decreasing = TRUE)]
}

# -----------------------------------------------------------------------------
# Mutable component object (used by the greedy search). It is an environment so
# that add / remove operations are reflected in place across recursive calls.
# -----------------------------------------------------------------------------

#' Create a mutable component seeded with a single node
#'
#' @param sc A score context.
#' @param start_node A single node name to start from.
#'
#' @return An environment representing the component.
.component_new <- function(sc, start_node) {
  comp <- new.env(parent = emptyenv())
  comp$sc <- sc
  comp$nodes <- start_node
  comp$members <- new.env(parent = emptyenv())
  assign(start_node, TRUE, envir = comp$members)
  comp$zsum <- sc$z[[start_node]]
  comp$score <- .score_subnetwork(sc, 1L, comp$zsum, TRUE)
  comp
}

#' Whether a component contains a node
#' @noRd
.component_contains <- function(comp, node) {
  exists(node, envir = comp$members, inherits = FALSE)
}

#' Add a node to a component, updating its score in place
#' @noRd
.component_add <- function(comp, node) {
  comp$nodes <- c(comp$nodes, node)
  assign(node, TRUE, envir = comp$members)
  comp$zsum <- comp$zsum + comp$sc$z[[node]]
  comp$score <- .score_subnetwork(comp$sc, length(comp$nodes), comp$zsum, TRUE)
  invisible(NULL)
}

#' Remove a node from a component, updating its score in place
#' @noRd
.component_remove <- function(comp, node) {
  if (.component_contains(comp, node)) {
    comp$nodes <- comp$nodes[comp$nodes != node]
    rm(list = node, envir = comp$members)
    comp$zsum <- comp$zsum - comp$sc$z[[node]]
    comp$score <- .score_subnetwork(comp$sc, length(comp$nodes), comp$zsum, TRUE)
  }
  invisible(NULL)
}
