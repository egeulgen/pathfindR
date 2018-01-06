.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
"##############################################################################
                        Welcome to pathfindr
##############################################################################")
}

#' pathfindr: A package for Pathway Enrichment Analysis Utilizing Active
#' Subnetworks
#'
#' The pathfindr package provides two important functions: \code{pathfindr} and
#' \code{choose_clusters}.
#'
#' @section pathfindr: This function is the wrapper function for the pathfindr
#'   workflow. It takes in a data frame consisting of Gene Symbol,
#'   log-fold-change and adjusted-p values. After input testing, any gene
#'   symbols that are not in the PIN are converted to alias symbols if the alias
#'   is in the PIN. Next, active subnetwork search is performed. Pathway
#'   enrichment analysis is performed using the genes in each of the active
#'   subnetworks. Pathways with adjusted-p values lower than
#'   \code{enrichment_threshold} are discarded. The lowest adjusted-p value
#'   (over all subnetworks) for each pathway is kept. This process of active
#'   subnetwork search and enrichment is repeated  for a selected number of
#'   \code{iterations}, which is done in parallel. Over all iterations, the
#'   lowest and the highest adjusted-p values, as well as number of occurences
#'   are reported for each enriched pathway.
#'
#' @section choose_clusters: This function first calculates the pairwise
#'   distances between the pathways in the \code{result_df} data frame. Via a
#'   shiny HTML document, the hierarchical clustering dendrogram is visualized.
#'   In this HTML document, the user can select the value at which to cut the
#'   tree and the resulting representative pathways (chosen by smallest lowest p
#'   value) are presented as a table and pathways with cluster assignments are
#'   saved as a csv file to the current directory.
#'
#' @seealso See \code{\link{run_pathfindr}} and \code{\link{choose_clusters}}
#'   for more details.
#'
#' @docType package
#' @name pathfindr
NULL
