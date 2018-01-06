#' Get Current KEGG Pathway Genes
#'
#' @param kegg_update boolean value indicating whether or not to update KEGG
#'   pathway genes (default = FALSE)
#'
#' @return list containing genes for each pathway. The names of the list are
#'   KEGG IDs.
#' @export
#' @seealso \code{\link{enrichment}} for pathway enrichment analysis. See
#'   \code{\link{run_pathfindr}} for the wrapper function of the pathfindr
#'   workflow
#' @examples
#' \dontrun{
#' current_KEGG()
#' }
current_KEGG <- function(kegg_update = FALSE) {
  if (!file.exists("KEGG_pws.Rdata") | kegg_update) {
    # created named list, eg:  path:map00010: "Glycolysis / Gluconeogenesis"
    pathways_list <- KEGGREST::keggList("pathway", "hsa")

    # make them into KEGG-style human pathway identifiers
    pathway_codes <- sub("path:", "", names(pathways_list))

    # parse pathway genes
    genes_by_pathway <- sapply(pathway_codes, function(pwid){
      pw <- KEGGREST::keggGet(pwid)
      pw <- pw[[1]]$GENE[c(FALSE, TRUE)] ## get gene symbols
      pw <- sub(";.+", "", pw)
      pw <- pw[grepl("^[a-zA-Z0-9_-]*$", pw)] ## remove mistaken lines
      pw
    })
    save(genes_by_pathway, file = "KEGG_pws.RData")
  }
  else
    load("KEGG_pws.RData")
  return(genes_by_pathway)
}
