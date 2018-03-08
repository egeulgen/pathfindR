#' Annotate Involved Genes In Pathways and Visualize Pathways
#'
#' @param pw_table Data frame of enrichment results. Columns are: "ID",
#'   "Pathway", "occurrence", "lowest_p", "highest_p".
#' @param gene_data Single column data frame containing change values (e.g.
#'   log(fold change) values) for significant genes.
#'
#' @return Data frame of enrichment results with genes involved in each pathway
#'   presented. Columns are: "ID", "Pathway", "occurrence", "lowest_p",
#'   "highest_p","Involved_genes". The function also creates visualizations of
#'   the pathways with the package \code{pathview} and saves them in the folder
#'   "pathway_maps" in the folder "pathfindr_Results" under the current working
#'   directory.
#'
#' @export
#' @import pathview
#' @seealso \code{\link[pathview]{pathview}} for pathway-based data integration
#'   and visualization. See \code{\link{run_pathfindR}} for the wrapper function
#'   of the pathfindR workflow
#' @examples
#' \dontrun{
#' pathmap(pathway_table, change_data)
#' }
pathmap <- function(pw_table, gene_data) {

  pw_table$Up_regulated <- ""
  pw_table$Down_regulated <- ""

  ## fix KEGG names such as "Glycolysis / Gluconeogenesis"
  pw_table$Pathway <- gsub("\\/", "-", pw_table$Pathway)

  upreg <- rownames(gene_data)[gene_data > 0] ## need log2 ratio
  downreg <- rownames(gene_data)[gene_data < 0]

  dir.create("pathway_maps")
  setwd("pathway_maps")

  for (i in 1:nrow(pw_table)) {
    tmp <- pathview::pathview(gene.data = gene_data,
                              gene.idtype = "SYMBOL",
                              pathway.id = pw_table$ID[i],
                              species = "hsa",
                              out.suffix = pw_table$Pathway[i],
                              keys.align = "y", kegg.native = TRUE,
                              key.pos = "topright", same.layer = FALSE)

   if (length(tmp) != 1) {
     tmp <- tmp$plot.data.gene$all.mapped
     tmp <- tmp[tmp != ""]
     tmp <- unlist(strsplit(tmp, split = ","))
     tmp <- unique(tmp)
     if (any(tmp != ""))
       tmp <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, tmp,
                                    "SYMBOL", "ENTREZID")[, 2]

     pw_table$Up_regulated[i] <- paste(tmp[tmp %in% upreg], collapse = ", ")
     pw_table$Down_regulated[i] <- paste(tmp[tmp %in% downreg], collapse = ", ")
   }
  }
  setwd("..")
  return(pw_table)
}
