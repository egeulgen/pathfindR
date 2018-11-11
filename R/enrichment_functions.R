#' Hypergeometic Distribution-based Hypothesis Testing
#'
#' @param pw_genes vector of genes in the selected pathway
#' @param chosen_genes vector containing the set of input genes
#' @param all_genes vector of "all" genes (background)
#'
#' @return the p-value as determined using the hypergeometric distribution.
#'
#' @details To determine whether the `chosen genes` are enriched
#' (compared to a background pool of genes) in the `pw_genes`, the
#' hypergeometric distribution is assumed and the appropriate p value
#' (the value under the right tail) is calculated and returned.
#'
#'
#' @export
#'
#' @examples
#' hyperg_test(letters[1:5], letters[2:5], letters)
#' hyperg_test(letters[1:5], letters[2:10], letters)
#' hyperg_test(letters[1:5], letters[2:13], letters)
hyperg_test <- function(pw_genes, chosen_genes, all_genes) {
  pw_genes_selected <- length(intersect(chosen_genes, pw_genes))
  pw_genes_in_pool <- length(pw_genes)
  tot_genes_in_pool <- length(all_genes)
  non_pw_genes_in_pool <- tot_genes_in_pool - pw_genes_in_pool
  num_selected_genes <- length(chosen_genes)

  p_val <- stats::phyper(pw_genes_selected - 1, pw_genes_in_pool,
                         non_pw_genes_in_pool, num_selected_genes,
                         lower.tail = FALSE)
  return(p_val)
}

#' Perform Enrichment Analysis for a Single Gene Set
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
#' @param DEG_vec vector of differentially-expressed gene symbols
#'
#' @return A data frame that contains enrichment results.
#' @export
#' @seealso \code{\link[stats]{p.adjust}} for adjustment of p values. See
#'   \code{\link{run_pathfindR}} for the wrapper function of the pathfindR
#'   workflow. \code{\link{hyperg_test}} for the details on hypergeometric
#'   distribution-based hypothesis testing.
#' @examples
#' pin_path <- return_pin_path("KEGG")
#' enrichment(kegg_genes, c("PER1", "PER2", "CRY1", "CREB1"), kegg_pathways,
#'            "bonferroni", 0.05, pin_path, c("PER1"))
enrichment <- function(genes_by_pathway, genes_of_interest,
                       pathways_list, adj_method = "bonferroni",
                       enrichment_threshold, pin_path, DEG_vec) {

  pin <- utils::read.delim(file = pin_path, header = FALSE,
                           stringsAsFactors = FALSE)
  all_genes <- unique(c(pin$V1, pin$V2))

  ## Hypergeometric test for p value
  enrichment_res <- sapply(genes_by_pathway, pathfindR::hyperg_test,
                             genes_of_interest, all_genes)
  enrichment_res <- as.data.frame(enrichment_res)
  colnames(enrichment_res) <- "p_value"

  ## Fold enrinchment
  fe_calc <- function(x) {
    A <- sum(genes_of_interest %in% x) / length(genes_of_interest)
    B <- sum(all_genes %in% x) / length(all_genes)
    return(A / B)
  }
  enrichment_res$Fold_Enrichment <- sapply(genes_by_pathway, fe_calc)


  idx <- order(enrichment_res$p_value)
  enrichment_res <- enrichment_res[idx,, drop = FALSE]
  enrichment_res$adj_p <- stats::p.adjust(enrichment_res$p, method = adj_method)

  if (all(enrichment_res$adj_p > enrichment_threshold))
    return(NULL)
  else {
    cond <- enrichment_res$adj_p <= enrichment_threshold
    enrichment_res <- enrichment_res[cond, ]

    ## pathway IDs
    enrichment_res$ID <- rownames(enrichment_res)

    ## Pathway descriptions
    idx <- match(enrichment_res$ID, names(pathways_list))
    enrichment_res$Pathway <- pathways_list[idx]

    ## non-DEG Active Subnetwork Genes
    tmp <- genes_of_interest[!genes_of_interest %in% DEG_vec]

    for (i in 1:nrow(enrichment_res)) {
      tmp2 <- tmp[tmp %in% genes_by_pathway[[enrichment_res$ID[i]]]]
      enrichment_res$non_DEG_Active_Snw_Genes[i] <- paste(tmp2, collapse = ", ")
    }

    ## reorder columns
    to_order <- c("ID", "Pathway", "Fold_Enrichment",
                  "p_value", "adj_p", "non_DEG_Active_Snw_Genes")
    enrichment_res <- enrichment_res[, to_order]

    return(enrichment_res)
  }
}
