current_KEGG <- function() {

  # created named list, eg:  path:map00010: "Glycolysis / Gluconeogenesis"
  pathways_list <- KEGGREST::keggList("pathway", "hsa")

  # make them into KEGG-style human pathway identifiers
  pathway_codes <- sub("path:", "", names(pathways_list))

  # parse pathway genes
  genes_by_pathway <- sapply(pathway_codes, function(pwid){
    pw <- KEGGREST::keggGet(pwid)
    pw <- pw[[1]]$GENE[c(F, T)]
    pw <- sub(";.+", "", pw)
    pw <- pw[grepl("^[a-zA-Z0-9_-]*$", pw)] ## removing mistaken lines
    pw
  })
  return(genes_by_pathway)
}
