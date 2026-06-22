#' Build the score context
#'
#' Computes per-node z-scores from the experiment p-values and the Monte-Carlo
#' mean / standard deviation of the raw subnetwork score at every possible
#' subnetwork size, used to calibrate (normalise) scores.
#'
#' The Monte-Carlo step is an exact replica of the Java reference's
#' \code{ScoreCalculations.calculateMeanAndStdForMonteCarlo}: it shuffles the
#' z-score vector 2000 times with a \code{java.util.Random}-equivalent generator
#' and Fisher-Yates shuffle (see \code{get_java_mc_calibration()} in
#' \file{score_calculation_utils.cpp}). Because \code{Collections.shuffle} starts
#' from \code{networkNodeList} order, the z-score vector MUST be in Java node
#' order; \code{network$nodes} (from \code{build_network()}) provides exactly that,
#' so \code{z} is built over it directly. With the matching seed this reproduces
#' the Java means/stds to floating-point precision, which is what aligns the
#' calibrated scores and rankings with the Java output.
#'
#' P-values are clamped to a safe range, the smallest p-value is kept when a
#' gene appears more than once, and network nodes without a p-value are given
#' \code{params$p_for_nonsignificant}.
#'
#' @param network A network as returned by \code{build_network()}.
#' @param experiment A data frame of gene / p-value pairs.
#' @param params A list of run parameters. \code{params$seed} must match the Java
#'   \code{-seedForRandom} value (Java default 1234) for identical calibration.
#'
#' @return A list with elements \code{z} (named z-score vector),
#'   \code{means} and \code{stds} (numeric vectors indexed by subnetwork size),
#'   and \code{nodes}.
#'
#' @export
build_score_context <- function(network, experiment, params) {
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

  # Use the Java approximation for Z-scores
  # this is to match version <=2.X.Y behavior
  z_vec <- get_java_zscores(as.numeric(pmap))
  z <- stats::setNames(z_vec, nodes)

  # Exact Java Monte-Carlo replica. z_vec is in Java node order (network$nodes)
  mc_data <- get_java_mc_calibration(z_vec, trials = 2000, seed = params$seed)

  list(
    z     = z,
    means = mc_data$means,
    stds  = mc_data$stds,
    nodes = nodes
  )
}

#' Score a subnetwork from its size and z-score sum
#'
#' Single-node subnetworks always score 0. Otherwise the raw score is
#' \code{zsum / sqrt(n)}, optionally calibrated against the Monte-Carlo
#' distribution and optionally penalized for size.
#'
#' @param sc A score context from \code{build_score_context()}.
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
#' @param network A network from \code{build_network()}.
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
#' @param network A network from \code{build_network()}.
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
