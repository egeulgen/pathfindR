#' Annotate Involved Genes In Pathways and Visualize Pathways
#'
#' @param pw_table Data frame of enrichment results. Columns are: "ID",
#'   "Pathway", "occurence", "lowest_p", "highest_p".
#' @param gene_data Single column data frame containing change values (e.g.
#'   log(fold change) values) for significant genes.
#'
#' @return Data frame of enrichment results with genes involved in each pathway
#'   presented. Columns are: "ID", "Pathway", "occurence", "lowest_p",
#'   "highest_p","Involved_genes". The function also creates visualizations of
#'   the pathways with the package "pathview".
#' @export
#' @seealso \code{\link[pathview]{pathview}} for pathway based data integration
#'   and visualization.
#' @examples
#' pathmap(pathway_table, change_data)
pathmap <- function(pw_table, gene_data) {
  suppressPackageStartupMessages(library(pathview)) ## cannot find work-around

  ## Load required data for pathview
  data("gene.idtype.bods", package = "pathview")
  data("gene.idtype.list", package = "pathview")
  data("cpd.simtypes", package = "pathview")

  pw_table$Involved_genes <- ""
  pw_table$Pathway <- gsub("\\/", "-", pw_table$Pathway) ## fix kegg names such as "Glycolysis / Gluconeogenesis"

  dir.create("pathway_maps")
  setwd("pathway_maps")

  for (i in 1:nrow(pw_table)) {
    tmp <- pathview(gene.data = gene_data,
                              gene.idtype = "SYMBOL",
                              pathway.id = pw_table$ID[i],
                              species = "hsa",
                              out.suffix = pw_table$Pathway[i],
                              keys.align = "y", kegg.native = T,
                              key.pos = "topright", same.layer = F)

    tmp <- tmp$plot.data.gene$all.mapped
    tmp <- tmp[tmp != ""]
    tmp <- unlist(strsplit(tmp, split = ","))
    tmp <- unique(tmp)
    if (any(tmp != ""))
      tmp <- AnnotationDbi::select(org.Hs.eg.db:::org.Hs.eg.db, tmp,
                                   "SYMBOL", "ENTREZID")[, 2]

    pw_table$Involved_genes[i] <- paste(tmp, collapse = ", ")
  }
  setwd("..")
  return(pw_table)
}
