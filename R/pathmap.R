pathmap <- function(pw_table, gene_data) {
  ## Load required data for pathview
  data("gene.idtype.bods", package = "pathview")
  data("gene.idtype.list", package = "pathview")
  data("cpd.simtypes", package = "pathview")

  pw_table$Involved_genes <- ""
  pw_table$Pathway <- gsub("\\/", "-", pw_table$Pathway) ## fix kegg names such as "Glycolysis / Gluconeogenesis"

  dir.create("pathway_maps")
  setwd("pathway_maps")

  for (i in 1:nrow(pw_table)) {
    tmp <- pathview::pathview(gene.data = gene_data,
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
  ## remove required data for pathview
  rm(gene.idtype.bods, gene.idtype.list, cpd.simtypes)
  setwd("..")
  return(pw_table)
}
