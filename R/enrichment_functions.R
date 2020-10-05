#' Hypergeometric Distribution-based Hypothesis Testing
#'
#' @param term_genes vector of genes in the selected term gene set
#' @param chosen_genes vector containing the set of input genes
#' @param background_genes vector of background genes (i.e. universal set of
#' genes in the experiment)
#'
#' @return the p-value as determined using the hypergeometric distribution.
#'
#' @details To determine whether the \code{chosen_genes} are enriched
#' (compared to a background pool of genes) in the \code{term_genes}, the
#' hypergeometric distribution is assumed and the appropriate p value
#' (the value under the right tail) is calculated and returned.
#'
#' @export
#'
#' @examples
#' hyperg_test(letters[1:5], letters[2:5], letters)
#' hyperg_test(letters[1:5], letters[2:10], letters)
#' hyperg_test(letters[1:5], letters[2:13], letters)
hyperg_test <- function(term_genes, chosen_genes, background_genes) {

  #### Argument checks
  if (!is.atomic(term_genes)) {
    stop("`term_genes` should be a vector")
  }
  if (!is.atomic(chosen_genes)) {
    stop("`chosen_genes` should be a vector")
  }
  if(!is.atomic(background_genes)) {
    stop("`background_genes` should be a vector")
  }

  if (length(term_genes) > length(background_genes)) {
    stop("`term_genes` cannot be larger than `background_genes`!")
  }
  if (length(chosen_genes) > length(background_genes)) {
    stop("`chosen_genes` cannot be larger than `background_genes`!")
  }

  #### Calculate p value
  term_genes_selected <- sum(chosen_genes %in% term_genes)
  term_genes_in_pool <- sum(term_genes %in% background_genes)
  tot_genes_in_pool <- length(background_genes)
  non_term_genes_in_pool <- tot_genes_in_pool - term_genes_in_pool
  num_selected_genes <- length(chosen_genes)

  p_val <- stats::phyper(term_genes_selected - 1, term_genes_in_pool,
                         non_term_genes_in_pool, num_selected_genes,
                         lower.tail = FALSE)
  return(p_val)
}

#' Perform Enrichment Analysis for a Single Gene Set
#'
#' @param input_genes The set of gene symbols to be used for enrichment
#'   analysis. In the scope of this package, these are genes that were
#'   identified for an active subnetwork
#' @param genes_by_term List that contains genes for each gene set. Names of
#'   this list are gene set IDs (default = kegg_genes)
#' @param term_descriptions Vector that contains term descriptions for the
#'   gene sets. Names of this vector are gene set IDs (default = kegg_descriptions)
#' @param adj_method correction method to be used for adjusting p-values.
#'   (default = "bonferroni")
#' @param enrichment_threshold adjusted-p value threshold used when filtering
#'   enrichment results (default = 0.05)
#' @param sig_genes_vec vector of significant gene symbols. In the scope of this
#'   package, these are the input genes that were used for active subnetwork search
#' @param background_genes vector of background genes. In the scope of this package,
#'   the background genes are taken as all genes in the PIN
#'   (see \code{\link{enrichment_analyses}})
#'
#' @return A data frame that contains enrichment results
#' @export
#' @seealso \code{\link[stats]{p.adjust}} for adjustment of p values. See
#'   \code{\link{run_pathfindR}} for the wrapper function of the pathfindR
#'   workflow. \code{\link{hyperg_test}} for the details on hypergeometric
#'   distribution-based hypothesis testing.
#' @examples
#' enrichment(input_genes = c("PER1", "PER2", "CRY1", "CREB1"),
#'            sig_genes_vec = "PER1",
#'            background_genes = unlist(pathfindR.data::kegg_genes))
enrichment <- function(input_genes,
                       genes_by_term = pathfindR.data::kegg_genes,
                       term_descriptions = pathfindR.data::kegg_descriptions,
                       adj_method = "bonferroni",
                       enrichment_threshold = 5e-2,
                       sig_genes_vec,
                       background_genes) {

  #### Argument checks
  ## input genes
  if (!is.atomic(input_genes)) {
    stop("`input_genes` should be a vector of gene symbols")
  }

  ## gene sets data
  if (!is.list(genes_by_term)) {
    stop("`genes_by_term` should be a list of term gene sets")
  }
  if (is.null(names(genes_by_term))) {
    stop("`genes_by_term` should be a named list (names are gene set IDs)")
  }

  if (!is.atomic(term_descriptions)) {
    stop("`term_descriptions` should be a vector of term gene descriptions")
  }
  if (is.null(names(term_descriptions))) {
    stop("`term_descriptions` should be a named vector (names are gene set IDs)")
  }

  if (length(genes_by_term) != length(term_descriptions)) {
    stop("The lengths of `genes_by_term` and `term_descriptions` should be the same")
  }
  if (any(names(genes_by_term) != names(term_descriptions))) {
    stop("The names of `genes_by_term` and `term_descriptions` should all be the same")
  }

  ## enrichment threshold
  if (!is.numeric(enrichment_threshold)) {
    stop("`enrichment_threshold` should be a numeric value between 0 and 1")
  }
  if (enrichment_threshold < 0 | enrichment_threshold > 1) {
    stop("`enrichment_threshold` should be between 0 and 1")
  }

  ## signif. genes and background (universal set) genes
  if (!is.atomic(sig_genes_vec)) {
    stop("`sig_genes_vec` should be a vector")
  }
  if (!is.atomic(background_genes)) {
    stop("`background_genes` should be a vector")
  }

  #### Obtain p values
  enrichment_res <- vapply(genes_by_term, pathfindR::hyperg_test, 0.1,
                           input_genes, background_genes)
  enrichment_res <- as.data.frame(enrichment_res)
  colnames(enrichment_res) <- "p_value"

  # Adjust p values
  idx <- order(enrichment_res$p_value)
  enrichment_res <- enrichment_res[idx, , drop = FALSE]
  enrichment_res$adj_p <- stats::p.adjust(enrichment_res$p, method = adj_method)


  #### Filter by adj-p
  cond <- enrichment_res$adj_p <= enrichment_threshold
  # Empty case (if all adj-p > threshold)
  if (sum(cond) == 0)
    return(NULL)
  enrichment_res <- enrichment_res[cond, ]

  #### Add other columns
  # Term IDs
  enrichment_res$ID <- rownames(enrichment_res)

  ## Term descriptions
  idx <- match(enrichment_res$ID, names(term_descriptions))
  enrichment_res$Term_Description <- term_descriptions[idx]

  # Fold enrinchment
  gset_for_fe <- genes_by_term[rownames(enrichment_res)]
  A <- vapply(gset_for_fe, function(gset) length(intersect(sig_genes_vec, gset)), 1L) / length(sig_genes_vec)
  B <- vapply(gset_for_fe, function(gset) length(intersect(background_genes, gset)), 1L) / length(background_genes)
  enrichment_res$Fold_Enrichment <- A / B

  # Non-significant Subnetwork Genes
  non_sig_snw_genes <- base::setdiff(input_genes, sig_genes_vec)
  for (i in base::seq_len(nrow(enrichment_res))) {
    tmp <- intersect(non_sig_snw_genes, genes_by_term[[enrichment_res$ID[i]]])
    enrichment_res$non_Signif_Snw_Genes[i] <- paste(tmp, collapse = ", ")
  }

  ## reorder columns
  to_order <- c("ID", "Term_Description", "Fold_Enrichment",
                "p_value", "adj_p", "non_Signif_Snw_Genes")
  enrichment_res <- enrichment_res[, to_order]

  return(enrichment_res)
}

#' Perform Enrichment Analyses on the Input Subnetworks
#'
#' @param snws a list of subnetwork genes (i.e., vectors of genes for each subnetwork)
#' @inheritParams enrichment
#' @inheritParams return_pin_path
#' @param list_active_snw_genes boolean value indicating whether or not to report
#' the non-significant active subnetwork genes for the active subnetwork which was enriched for
#' the given term with the lowest p value (default = \code{FALSE})
#'
#' @return a dataframe of combined enrichment results. Columns are: \describe{
#'   \item{ID}{ID of the enriched term}
#'   \item{Term_Description}{Description of the enriched term}
#'   \item{Fold_Enrichment}{Fold enrichment value for the enriched term}
#'   \item{p_value}{p value of enrichment}
#'   \item{adj_p}{adjusted p value of enrichment}
#'   \item{support}{the support (proportion of active subnetworks leading to enrichment over all subnetworks) for the gene set}
#'   \item{non_Signif_Snw_Genes (OPTIONAL)}{the non-significant active subnetwork genes, comma-separated}
#' }
#'
#' @export
#'
#' @seealso \code{\link{enrichment}} for the enrichment analysis for a single gene set
#'
#' @examples
#' enr_res <- enrichment_analyses(snws = example_active_snws[1:2],
#'                                sig_genes_vec = RA_input$Gene.symbol[1:25],
#'                                pin_name_path = "KEGG")
enrichment_analyses <- function(snws,
                                sig_genes_vec,
                                pin_name_path = "Biogrid",
                                genes_by_term = pathfindR.data::kegg_genes,
                                term_descriptions = pathfindR.data::kegg_descriptions,
                                adj_method = "bonferroni",
                                enrichment_threshold = 0.05,
                                list_active_snw_genes = FALSE) {

  ### Argument check
  if (!is.logical(list_active_snw_genes)) {
    stop("`list_active_snw_genes` should be either TRUE or FALSE")
  }

  ### Load PIN Data
  pin_path <- return_pin_path(pin_name_path)
  pin <- utils::read.delim(file = pin_path, header = FALSE,
                           stringsAsFactors = FALSE)

  background_genes <- unique(c(pin[, 1], pin[, 3]))

  # turn all to upper case for best match
  genes_by_term <- lapply(genes_by_term, base::toupper)
  sig_genes_vec <- base::toupper(sig_genes_vec)
  background_genes <- base::toupper(background_genes)

  ############ Enrichment per subnetwork
  enrichment_res <- lapply(snws, function(x)
    pathfindR::enrichment(input_genes = base::toupper(x),
                          genes_by_term = genes_by_term,
                          term_descriptions = term_descriptions,
                          adj_method = adj_method,
                          enrichment_threshold = enrichment_threshold,
                          sig_genes_vec = sig_genes_vec,
                          background_genes = background_genes))

  ### indices for snw.s
  if (length(enrichment_res) != 0) {
    for (i in 1:length(enrichment_res)) {
      if (!is.null(enrichment_res[[i]])) {
        enrichment_res[[i]]$snw_idx <- i
      }
    }
  }

  ############ Combine Enrichments Results for All Subnetworks
  enrichment_res <- Reduce(rbind, enrichment_res)

  ############ Process if non-empty
  if (!is.null(enrichment_res)) {

    ## calculate support values
    support <- tapply(enrichment_res$snw_idx, enrichment_res$ID, length)
    support <- support / length(snws)
    enrichment_res$support <- support[match(enrichment_res$ID, names(support))]
    enrichment_res$snw_idx <- NULL

    ## delete non_Signif_Snw_Genes if list_active_snw_genes == FALSE
    if (!list_active_snw_genes) {
      enrichment_res$non_Signif_Snw_Genes <- NULL
    }

    ## keep lowest p for each term
    idx <- order(enrichment_res$adj_p)
    enrichment_res <- enrichment_res[idx, ]
    enrichment_res <- enrichment_res[!duplicated(enrichment_res$ID), ]
  }
  return(enrichment_res)
}


#' Summarize Enrichment Results
#'
#' @param enrichment_res a dataframe of combined enrichment results. Columns are: \describe{
#'   \item{ID}{ID of the enriched term}
#'   \item{Term_Description}{Description of the enriched term}
#'   \item{Fold_Enrichment}{Fold enrichment value for the enriched term}
#'   \item{p_value}{p value of enrichment}
#'   \item{adj_p}{adjusted p value of enrichment}
#'   \item{non_Signif_Snw_Genes (OPTIONAL)}{the non-significant active subnetwork genes, comma-separated}
#' }
#' @inheritParams enrichment_analyses
#'
#' @return a dataframe of summarized enrichment results (over multiple iterations). Columns are: \describe{
#'   \item{ID}{ID of the enriched term}
#'   \item{Term_Description}{Description of the enriched term}
#'   \item{Fold_Enrichment}{Fold enrichment value for the enriched term}
#'   \item{occurrence}{the number of iterations that the given term was found to enriched over all iterations}
#'   \item{support}{the median support (proportion of active subnetworks leading to enrichment within an iteration) over all iterations}
#'   \item{lowest_p}{the lowest adjusted-p value of the given term over all iterations}
#'   \item{highest_p}{the highest adjusted-p value of the given term over all iterations}
#'   \item{non_Signif_Snw_Genes (OPTIONAL)}{the non-significant active subnetwork genes, comma-separated}
#' }
#' @export
#'
#' @examples
#' \dontrun{
#' summarize_enrichment_results(enrichment_res)
#' }
summarize_enrichment_results <- function(enrichment_res,
                                         list_active_snw_genes = FALSE) {
  ## Argument checks
  if(!is.logical(list_active_snw_genes)) {
    stop("`list_active_snw_genes` should be either TRUE or FALSE")
  }

  nec_cols <- c("ID", "Term_Description", "Fold_Enrichment", "p_value", "adj_p", "support")
  if (list_active_snw_genes) {
    nec_cols <- c(nec_cols, "non_Signif_Snw_Genes")
  }

  if (!is.data.frame(enrichment_res)) {
    stop("`enrichment_res` should be a data frame")
  }

  if (ncol(enrichment_res) != length(nec_cols)) {
    stop("`enrichment_res` should have exactly ", length(nec_cols), " columns")
  }

  if (!all(nec_cols %in% colnames(enrichment_res))) {
    stop("`enrichment_res` should have column names ",
         paste(dQuote(nec_cols), collapse = ", "))
  }

  ## Annotate lowest p, highest p, occurrence and median support
  final_res <- enrichment_res
  lowest_p <- tapply(enrichment_res$adj_p, enrichment_res$ID, min)
  highest_p <- tapply(enrichment_res$adj_p, enrichment_res$ID, max)
  occurrence <- tapply(enrichment_res$adj_p, enrichment_res$ID, length)
  support <- tapply(enrichment_res$support, enrichment_res$ID, stats::median)

  matched_idx <- match(final_res$ID, names(lowest_p))
  final_res$lowest_p <- as.numeric(lowest_p[matched_idx])

  matched_idx <- match(final_res$ID, names(highest_p))
  final_res$highest_p <- as.numeric(highest_p[matched_idx])

  matched_idx <- match(final_res$ID, names(occurrence))
  final_res$occurrence <- as.numeric(occurrence[matched_idx])

  matched_idx <- match(final_res$ID, names(support))
  final_res$support <- as.numeric(support[matched_idx])

  ## reorder columns
  keep <- c("ID", "Term_Description", "Fold_Enrichment",
            "occurrence", "support", "lowest_p", "highest_p")
  if (list_active_snw_genes) {
    keep <- c(keep, "non_Signif_Snw_Genes")
  }
  final_res <- final_res[, keep]

  ## keep data with lowest p-value over all iterations
  final_res <- final_res[order(final_res$lowest_p), ]
  final_res <- final_res[!duplicated(final_res$ID), ]
  rownames(final_res) <- NULL

  return(final_res)
}
