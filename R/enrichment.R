enrichment <- function(genes_by_pathway, genes_of_interest,
                       pathways_list, adj_method = "bonferroni") {
  hyperg <- Category:::.doHyperGInternal
  hyperg_test <- function(pw_genes, chosen_genes, all_genes, over=TRUE) {
    pw_genes_selected <- length(intersect(chosen_genes, pw_genes))
    pw_genes_in_pool <- length(pw_genes)
    tol_genes_n_pool <- length(all_genes)
    non_pw_genes_in_pool <- tol_genes_n_pool - pw_genes_in_pool
    num_selected_genes <- length(chosen_genes)
    hyperg(pw_genes_in_pool, non_pw_genes_in_pool,
           num_selected_genes, pw_genes_selected, over)
  }

  suppressPackageStartupMessages(library(org.Hs.eg.db))
  all_genes <- AnnotationDbi::keys(org.Hs.eg.db, "SYMBOL")

  enrichment_res <- t(sapply(genes_by_pathway, hyperg_test,
                             genes_of_interest, all_genes))
  enrichment_res <- as.data.frame(enrichment_res)

  enrichment_res$p <- unlist(enrichment_res$p)
  enrichment_res$odds <- unlist(enrichment_res$odds)
  enrichment_res$expected <- unlist(enrichment_res$expected)

  enrichment_res <- enrichment_res[order(enrichment_res$p), ]
  enrichment_res$adj_p <- p.adjust(enrichment_res$p, method = adj_method)

  ## pathway IDs
  enrichment_res$ID <- rownames(enrichment_res)

  ## Pathway desriptions
  idx <- match(paste0("path:", enrichment_res$ID), names(pathways_list))
  enrichment_res$Pathway <- pathways_list[idx]
  enrichment_res$Pathway <- sub(" - Homo sapiens \\(human\\)", "",
                                enrichment_res$Pathway)

  ## reorder columns
  remaining <- setdiff(colnames(enrichment_res), c("ID", "Pathway"))
  enrichment_res <- enrichment_res[, c("ID", "Pathway", remaining)]

  return(enrichment_res)
}
