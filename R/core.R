#' Wrapper Function for pathfindR - Active-Subnetwork-Oriented Enrichment Workflow
#'
#' \code{run_pathfindR} is the wrapper function for the pathfindR workflow
#'
#' This function takes in a data frame consisting of Gene Symbol, log-fold-change
#' and adjusted-p values. After input testing, any gene symbols that are not in
#' the PIN are converted to alias symbols if the alias is in the PIN. Next,
#' active subnetwork search is performed. Enrichment analysis is
#' performed using the genes in each of the active subnetworks. Terms with
#' adjusted-p values lower than \code{enrichment_threshold} are discarded. The
#' lowest adjusted-p value (over all subnetworks) for each term is kept. This
#' process of active subnetwork search and enrichment is repeated  for a selected
#' number of \code{iterations}, which is done in parallel. Over all iterations,
#' the lowest and the highest adjusted-p values, as well as number of occurrences
#' are reported for each enriched term.
#'
#' @inheritParams input_processing
#' @inheritParams fetch_gene_set
#' @inheritParams enrichment_analyses
#' @param plot_enrichment_chart boolean value. If TRUE, a bubble chart displaying
#'  the enrichment results is plotted. (default = TRUE)
#' @param output_dir the directory to be created where the output and intermediate
#'  files are saved (default = \code{NULL}, a temporary directory is used)
#' @param ... additional arguments for \code{\link{active_snw_enrichment_wrapper}}
#'
#' @return Data frame of pathfindR enrichment results. Columns are: \describe{
#'   \item{ID}{ID of the enriched term}
#'   \item{Term_Description}{Description of the enriched term}
#'   \item{Fold_Enrichment}{Fold enrichment value for the enriched term (Calculated using ONLY the input genes)}
#'   \item{occurrence}{the number of iterations that the given term was found to enriched over all iterations}
#'   \item{support}{the median support (proportion of active subnetworks leading to enrichment within an iteration) over all iterations}
#'   \item{lowest_p}{the lowest adjusted-p value of the given term over all iterations}
#'   \item{highest_p}{the highest adjusted-p value of the given term over all iterations}
#'   \item{non_Signif_Snw_Genes (OPTIONAL)}{the non-significant active subnetwork genes, comma-separated}
#'   \item{Up_regulated}{the up-regulated genes (as determined by `change value` > 0, if the `change column` was provided) in the input involved in the given term's gene set, comma-separated. If change column not provided, all affected are listed here.}
#'   \item{Down_regulated}{the down-regulated genes (as determined by `change value` < 0, if the `change column` was provided) in the input involved in the given term's gene set, comma-separated}
#' }
#'  The function also creates an HTML report with the pathfindR enrichment
#'  results linked to the visualizations of the enriched terms in addition to
#'  the table of converted gene symbols. This report can be found in
#'  '\code{output_dir}/results.html' under the current working directory.
#'
#'  By default, a bubble chart of top 10 enrichment results are plotted. The x-axis
#'  corresponds to fold enrichment values while the y-axis indicates the enriched
#'  terms. Sizes of the bubbles indicate the number of significant genes in the given terms.
#'  Color indicates the -log10(lowest-p) value; the more red it is, the more
#'  significant the enriched term is. See \code{\link{enrichment_chart}}.
#'
#' @import knitr
#' @import rmarkdown
#' @import parallel
#' @import doParallel
#' @import foreach
#' @import graphics
#'
#' @export
#'
#' @section Warning: Especially depending on the protein interaction network,
#'  the algorithm and the number of iterations you choose, 'active subnetwork
#'  search + enrichment' component of \code{run_pathfindR} may take a long time to finish.
#'
#' @seealso
#' \code{\link{input_testing}} for input testing, \code{\link{input_processing}} for input processing,
#' \code{\link{active_snw_search}} for active subnetwork search and subnetwork filtering,
#' \code{\link{enrichment_analyses}} for enrichment analysis (using the active subnetworks),
#' \code{\link{summarize_enrichment_results}} for summarizing the active-subnetwork-oriented enrichment results,
#' \code{\link{annotate_term_genes}} for annotation of affected genes in the given gene sets,
#' \code{\link{visualize_terms}} for visualization of enriched terms,
#' \code{\link{enrichment_chart}} for a visual summary of the pathfindR enrichment results,
#' \code{\link[foreach]{foreach}} for details on parallel execution of looping constructs,
#' \code{\link{cluster_enriched_terms}} for clustering the resulting enriched terms and partitioning into clusters.
#'
#' @examples
#' \dontrun{
#' run_pathfindR(example_pathfindR_input)
#' }
run_pathfindR <- function(input, gene_sets = "KEGG", min_gset_size = 10, max_gset_size = 300,
    custom_genes = NULL, custom_descriptions = NULL, pin_name_path = "Biogrid", p_val_threshold = 0.05,
    enrichment_threshold = 0.05, convert2alias = TRUE, plot_enrichment_chart = TRUE,
    output_dir = NULL, list_active_snw_genes = FALSE, ...) {
    ############ Argument checks
    if (!is.logical(plot_enrichment_chart)) {
        stop("`plot_enrichment_chart` should be either TRUE or FALSE")
    }
    if (!is.logical(list_active_snw_genes)) {
        stop("`list_active_snw_genes` should be either TRUE or FALSE")
    }

    gset_list <- fetch_gene_set(gene_sets = gene_sets, min_gset_size = min_gset_size,
        max_gset_size = max_gset_size, custom_genes = custom_genes, custom_descriptions = custom_descriptions)

    ## absolute path to PIN
    pin_path <- return_pin_path(pin_name_path)

    ## create output dir
    output_dir_org <- output_dir
    output_dir <- configure_output_dir(output_dir)
    # on exit, set working directory back to original working directory
    org_dir <- getwd()
    on.exit(setwd(org_dir))
    # create and change working directory into the output directory
    dir.create(output_dir, recursive = TRUE)
    output_dir <- normalizePath(output_dir)
    setwd(output_dir)

    input_testing(input, p_val_threshold)

    input_processed <- input_processing(input, p_val_threshold, pin_path, convert2alias)

    combined_res <- active_snw_enrichment_wrapper(input_processed, pin_path, gset_list,
        enrichment_threshold, list_active_snw_genes, ...)
    setwd(output_dir)

    ## In case no enrichment was found
    if (is.null(combined_res)) {
        warning("Did not find any enriched terms!", call. = FALSE)
        return(data.frame())
    }

    final_res <- summarize_enrichment_results(combined_res, list_active_snw_genes)


    final_res <- annotate_term_genes(result_df = final_res, input_processed = input_processed,
        genes_by_term = gset_list$genes_by_term)

    if (!is.null(output_dir_org)) {
        create_HTML_report(input = input, input_processed = input_processed, final_res = final_res,
            dir_for_report = output_dir)
    }

    if (plot_enrichment_chart) {
        graphics::plot(enrichment_chart(result_df = final_res))
    }

    message(paste0("Found ", nrow(final_res), " enriched terms\n\n"))
    message("You may run:\n")
    message("- cluster_enriched_terms() for clustering enriched terms\n")
    message("- visualize_terms() for visualizing enriched term diagrams\n\n")

    return(final_res)
}
