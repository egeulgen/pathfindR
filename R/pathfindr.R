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

#' pathfindR: A package for Enrichment Analysis Utilizing Active Subnetworks
#'
#' pathfindR is a tool for active-subnetwork-oriented gene set enrichment analysis.
#' The main aim of the package is to identify active subnetworks in a
#' protein-protein interaction network using a user-provided list of genes
#' and associated p values then performing enrichment analyses on the identified
#' subnetworks, discovering enriched terms (i.e. pathways, gene ontology, TF target
#' gene sets etc.) that possibly underlie the phenotype of interest.
#'
#' pathfindR also offers functionalities to cluster the enriched terms and
#' identify representative terms in each cluster, to score the enriched terms
#' per sample and to visualize analysis results.
#'
#' @seealso See \code{\link{run_pathfindR}} for details on the pathfindR
#' active-subnetwork-oriented enrichment analysis
#' See \code{\link{cluster_enriched_terms}} for details on methods of enriched
#' terms clustering to define clusters of biologically-related terms
#' See \code{\link{score_terms}} for details on agglomerated score calculation
#' for enriched terms to investigate how a gene set is altered in a given sample
#' (or in cases vs. controls)
#' See \code{\link{term_gene_graph}} for details on visualizing terms and
#' term-related genes as a graph to determine the degree of overlap between the
#' enriched terms by identifying shared and/or distinct significant genes
#'
#' @docType package
#' @name pathfindR
NULL
