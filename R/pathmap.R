pathmap <- function(pw_table, gene_data) {
  suppressPackageStartupMessages(library(org.Hs.eg.db))
  suppressPackageStartupMessages(library(pathview))
  pw_table$Involved_genes <- ""
  pw_table$Pathway <- gsub("\\/", "-", pw_table$Pathway)

  dir.create("pathway_maps")
  setwd("pathway_maps")

  for (i in 1:nrow(pw_table)) {

    tmp <- pathview(gene.data = gene_data, gene.idtype = "SYMBOL",
                    pathway.id = pw_table$ID[i], species = "hsa",
                    out.suffix = pw_table$Pathway[i],
                    keys.align = "y", kegg.native = T,
                    key.pos = "topright", same.layer = F)

    tmp <- tmp$plot.data.gene$all.mapped
    tmp <- tmp[tmp != ""]
    tmp <- unlist(strsplit(tmp, split = ","))
    tmp <- unique(tmp)
    if (any(tmp != ""))
      tmp <- select(org.Hs.eg.db, tmp, "SYMBOL", "ENTREZID")[, 2]

    pw_table$Involved_genes[i] <- paste(tmp, collapse = ", ")
  }
  setwd("..")
  detach("package:pathview", unload = TRUE)
  return(pw_table)
}
