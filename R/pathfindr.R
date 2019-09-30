.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "##############################################################################
                        Welcome to pathfindR!

Please cite the article below if you use pathfindR in published reseach:

Ulgen E, Ozisik O, Sezerman OU. 2019. pathfindR: An R Package for Comprehensive
Identification of Enriched Pathways in Omics Data Through Active Subnetworks.
Front. Genet. doi:10.3389/fgene.2019.00858

##############################################################################"
  )
}

#' pathfindR: A package for Pathway Enrichment Analysis Utilizing Active
#' Subnetworks
#'
#' The pathfindR package provides two important main functions: \code{run_pathfindR} and
#' \code{choose_clusters}.
#'
#' @section run_pathfindR: This function is the wrapper function for the pathfindR
#'   workflow. It takes in a data frame consisting of Gene Symbol,
#'   log-fold-change (optional) and adjusted-p values. After input testing, any gene
#'   symbols that are not in the PIN are converted to alias symbols if the alias
#'   is in the PIN. Next, active subnetwork search is performed. Enrichment
#'   analysis is performed using the genes in each of the active
#'   subnetworks. Enriched terms with adjusted-p values lower than
#'   \code{enrichment_threshold} are discarded. The lowest adjusted-p value
#'   (over all subnetworks) for each term is kept. This process of active
#'   subnetwork search and enrichment is repeated  for a selected number of
#'   \code{iterations}, which is executed in parallel. Over all iterations, the
#'   lowest and the highest adjusted-p values, as well as number of occurrences
#'   are reported for each enriched term.
#'
#' @section cluster_pathways: This function first calculates the pairwise
#'   kappa statistics between the terms in the \code{result_df} data frame. By default,
#'   hierarchical clustering is performed and the optimal number of clusters is chosen.
#'   Optionally, a fuzzy partitioning algorithm can also be used. The function returns
#'   a data frame with cluster assignments.
#'
#' @seealso See \code{\link{run_pathfindR}} and \code{\link{cluster_pathways}}
#'   for more details.
#'
#' @docType package
#' @name pathfindR
NULL
