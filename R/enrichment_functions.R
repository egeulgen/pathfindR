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

#' Perform Enrichment Analyses on the Input Subnetworks
#'
#' @param snws a list of subnetwork genes (i.e., vectors of genes for each subnetwork)
#' @param gene_sets the gene sets used for enrichment analysis. Available gene sets
#'  are KEGG, Reactome, BioCarta, GO-All, GO-BP, GO-CC, GO-MF or Custom. If "Custom", the arguments
#'  custom_genes and custom pathways must be specified. (Default = "KEGG")
#' @param custom_genes a list containing the genes involved in each custom pathway. Each element
#' is a vector of gene symbols located in the given pathway. Names correspond to
#' the ID of the pathway.
#' @param custom_pathways A list containing the descriptions for each custom pathway. Names of the
#' list correspond to the ID of the pathway.
#' @param pin_path path to the Protein Interaction Network (PIN) file used in
#'   the analysis
#' @param input_genes vector of input gene symbols (that were used during active subnetwork genes)
#' @param adj_method correction method to be used for adjusting p-values of
#'  pathway enrichment results (Default: 'bonferroni', see ?p.adjust)
#' @param enrichment_threshold threshold used when filtering individual iterations' pathway
#'  enrichment results
#' @param list_active_snw_genes boolean value indicating whether or not to report
#' the non-DEG active subnetwork genes for the active subnetwork which was enriched for
#' the given pathway with the lowest p value (default = FALSE)
#'
#'
#' @return a dataframe of combined enrichment results. Columns are: \describe{
#'   \item{ID}{KEGG ID of the enriched pathway}
#'   \item{Pathway}{Description of the enriched pathway}
#'   \item{Fold_Enrichment}{Fold enrichment value for the enriched pathway}
#'   \item{p_value}{p value of enrichment}
#'   \item{adj_p}{adjusted p value of enrichment}
#'   \item{non_DEG_Active_Snw_Genes (OPTIONAL)}{the non-DEG active subnetwork genes, comma-separated}
#' }
#'
#' @export
#'
#' @seealso \code{\link{enrichment}} for the enrichment analysis for a single gene set
#'
#' @examples
#' \dontrun{
#' enr_res <- enrichment_analyses(snws, input_genes = my_input$GENE,
#' gene_sets = "KEGG", pin_path = "path/to/PIN")
#' }
#'
enrichment_analyses <- function(snws, input_genes,
                                gene_sets = "KEGG", custom_genes = NULL, custom_pathways = NULL,
                                pin_path,
                                adj_method = "bonferroni", enrichment_threshold = 5e-2,
                                list_active_snw_genes = FALSE) {
  ############ Load Gene Set Data
  gene_sets_df <- data.frame("Gene Set Name" = c("KEGG", "Reactome", "BioCarta",
                                                 "GO-All", "GO-BP", "GO-CC", "GO-MF"),
                             "genes_by_pathway" = c("kegg_genes", "reactome_genes", "biocarta_genes",
                                                    "go_all_genes", "go_bp_genes", "go_cc_genes", "go_mf_genes"),
                             "pathways_list" = c("kegg_pathways", "reactome_pathways", "biocarta_pathways",
                                                 "go_all_pathways", "go_bp_pathways", "go_cc_pathways", "go_mf_pathways"))

  if (gene_sets %in% gene_sets_df$Gene.Set.Name) {
    idx <- which(gene_sets_df$Gene.Set.Name == gene_sets)
    genes_name <- gene_sets_df$genes_by_pathway[idx]
    pws_name <-  gene_sets_df$pathways_list[idx]

    genes_by_pathway <- base::eval(parse(text = paste0("pathfindR::", genes_name)))
    pathways_list <- base::eval(parse(text = paste0("pathfindR::", pws_name)))

  } else if (gene_sets == "Custom") {
    genes_by_pathway <- custom_genes
    pathways_list <- custom_pathways
  }

  ############ Enrichment per subnetwork
  enrichment_res <- lapply(snws, function(x)
    pathfindR::enrichment(genes_by_pathway, x, pathways_list,
                          adj_method, enrichment_threshold,
                          pin_path, DEG_vec = input_genes))

  ############ Combine Enrichments Results for All Subnetworks
  enrichment_res <- Reduce(rbind, enrichment_res)

  ############ Process if non-empty
  if (!is.null(enrichment_res)) {
    ## delete non_DEG_Active_Snw_Genes if list_active_snw_genes == FALSE
    if (!list_active_snw_genes)
      enrichment_res$non_DEG_Active_Snw_Genes <- NULL

    ## keep lowest p for each pathway
    idx <- order(enrichment_res$adj_p)
    enrichment_res <- enrichment_res[idx, ]
    enrichment_res <- enrichment_res[!duplicated(enrichment_res$ID), ]
  }
  return(enrichment_res)
}


#' Summarize Enrichment Results
#'
#' @param enrichment_res a dataframe of combined enrichment results. Columns are: \describe{
#'   \item{ID}{KEGG ID of the enriched pathway}
#'   \item{Pathway}{Description of the enriched pathway}
#'   \item{Fold_Enrichment}{Fold enrichment value for the enriched pathway}
#'   \item{p_value}{p value of enrichment}
#'   \item{adj_p}{adjusted p value of enrichment}
#'   \item{non_DEG_Active_Snw_Genes (OPTIONAL)}{the non-DEG active subnetwork genes, comma-separated}
#' }
#' @param list_active_snw_genes boolean value indicating whether or not to report
#' the non-DEG active subnetwork genes for the active subnetwork which was enriched for
#' the given pathway with the lowest p value (default = FALSE)
#'
#' @return a dataframe of summarized enrichment results (over multiple iterations). Columns are: \describe{
#'   \item{ID}{KEGG ID of the enriched pathway}
#'   \item{Pathway}{Description of the enriched pathway}
#'   \item{Fold_Enrichment}{Fold enrichment value for the enriched pathway}
#'   \item{occurrence}{the number of iterations that the given pathway was found to enriched over all iterations}
#'   \item{lowest_p}{the lowest adjusted-p value of the given pathway over all iterations}
#'   \item{highest_p}{the highest adjusted-p value of the given pathway over all iterations}
#'   \item{non_DEG_Active_Snw_Genes (OPTIONAL)}{the non-DEG active subnetwork genes, comma-separated}
#' }
#' @export
#'
#' @examples
#' \dontrun{
#' summarize_enrichment_results(enrichment_res)
#' }
summarize_enrichment_results <- function(enrichment_res, list_active_snw_genes = FALSE) {
  ## Annotate lowest p, highest p and occurrence
  final_res <- enrichment_res
  lowest_p <- tapply(enrichment_res$adj_p, enrichment_res$ID, min)
  highest_p <- tapply(enrichment_res$adj_p, enrichment_res$ID, max)
  occurrence <- tapply(enrichment_res$adj_p, enrichment_res$ID, length)

  matched_idx <- match(final_res$ID, names(lowest_p))
  final_res$lowest_p <- as.numeric(lowest_p[matched_idx])

  matched_idx <- match(final_res$ID, names(highest_p))
  final_res$highest_p <- as.numeric(highest_p[matched_idx])

  matched_idx <- match(final_res$ID, names(occurrence))
  final_res$occurrence <- as.numeric(occurrence[matched_idx])

  ## reorder columns
  keep <- c("ID", "Pathway", "Fold_Enrichment", "occurrence", "lowest_p", "highest_p")
  if (list_active_snw_genes)
    keep <- c(keep, "non_DEG_Active_Snw_Genes")
  final_res <- final_res[, keep]

  ## keep data with lowest p-value over all iterations
  final_res <- final_res[order(final_res$lowest_p), ]
  final_res <- final_res[!duplicated(final_res$ID), ]
  rownames(final_res) <- NULL

  return(final_res)
}
