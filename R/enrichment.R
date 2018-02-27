#' Perform Enrichment Analysis
#'
#' @param genes_by_pathway List that contains genes for each pathway. Names of
#'   this list are KEGG IDs.
#' @param genes_of_interest The set of gene symbols to be used for enrichment
#'   analysis. In the scope of this package, these are genes that were
#'   identified for an active subnetwork.
#' @param pathways_list List that contains pathway descriptions for KEGG pathway
#'   IDs. Names of this list are KEGG IDs.
#' @param adj_method correction method to be used for adjusting p-values.
#' @param enrichment_threshold adjusted-p value threshold used when filtering
#'   pathway enrichment results
#' @param pin_path path to the Protein-Protein Interaction Network (PIN) file used in
#'   the analysis
#'
#' @return A data frame that contains enrichment results.
#' @export
#' @seealso \code{\link[stats]{p.adjust}} for adjustment of p values. See
#'   \code{\link{run_pathfindR}} for the wrapper function of the pathfindR
#'   workflow.
#' @examples
#' pin_path <- return_pin_path("KEGG")
#' enrichment(genes_by_pathway, c("PER1", "PER2", "CRY1", "CREB1"), pathways_list,
#'            "bonferroni", 0.05, pin_path)
enrichment <- function(genes_by_pathway, genes_of_interest,
                       pathways_list, adj_method = "bonferroni",
                       enrichment_threshold, pin_path) {
  hyperg_test <- function(pw_genes, chosen_genes, all_genes) {
    pw_genes_selected <- length(intersect(chosen_genes, pw_genes))
    pw_genes_in_pool <- length(pw_genes)
    tot_genes_in_pool <- length(all_genes)
    non_pw_genes_in_pool <- tot_genes_in_pool - pw_genes_in_pool
    num_selected_genes <- length(chosen_genes)

    stats::phyper(pw_genes_selected - 1, pw_genes_in_pool,
                  non_pw_genes_in_pool, num_selected_genes,
                  lower.tail = FALSE)
  }

  pin <- utils::read.delim(file = pin_path, header = FALSE,
                           stringsAsFactors = FALSE)
  all_genes <- unique(c(pin$V1, pin$V2))

  enrichment_res <- sapply(genes_by_pathway, hyperg_test,
                             genes_of_interest, all_genes)
  enrichment_res <- as.data.frame(enrichment_res)
  colnames(enrichment_res) <- "p_value"

  idx <- order(enrichment_res$p_value)
  enrichment_res <- enrichment_res[idx, , drop = FALSE]
  enrichment_res$adj_p <- stats::p.adjust(enrichment_res$p, method = adj_method)

  if (all(enrichment_res$adj_p > enrichment_threshold))
    return(NULL)
  else {
    cond <- enrichment_res$adj_p <= enrichment_threshold
    enrichment_res <- enrichment_res[cond, ]

    ## pathway IDs
    enrichment_res$ID <- rownames(enrichment_res)

    ## Pathway descriptions
    idx <- match(paste0("path:", enrichment_res$ID), names(pathways_list))
    enrichment_res$Pathway <- pathways_list[idx]
    enrichment_res$Pathway <- sub(" - Homo sapiens \\(human\\)", "",
                                  enrichment_res$Pathway)

    ## reorder columns
    remaining <- setdiff(colnames(enrichment_res), c("ID", "Pathway"))
    enrichment_res <- enrichment_res[, c("ID", "Pathway", remaining)]

    return(enrichment_res)
  }
}
