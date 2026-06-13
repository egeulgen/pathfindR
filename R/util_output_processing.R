#' Annotate the Affected Genes in the Provided Enriched Terms
#'
#' Function to annotate the involved affected (input) genes in each term.
#'
#' @param result_df data frame of enrichment results.
#'  The only must-have column is 'ID'.
#' @param input_processed input data processed via \code{\link{input_processing}}
#' @param genes_by_term List that contains genes for each gene set. Names of
#'   this list are gene set IDs (default = kegg_genes)
#'
#' @return The original data frame with two additional columns:  \describe{
#'   \item{Up_regulated}{the up-regulated genes in the input involved in the given term's gene set, comma-separated}
#'   \item{Down_regulated}{the down-regulated genes in the input involved in the given term's gene set, comma-separated}
#' }
#' @export
#'
#' @examples
#' example_gene_data <- example_pathfindR_input
#' colnames(example_gene_data) <- c("GENE", "CHANGE", "P_VALUE")
#'
#' annotated_result <- annotate_term_genes(
#'   result_df = example_pathfindR_output,
#'   input_processed = example_gene_data
#' )
annotate_term_genes <- function(result_df, input_processed, genes_by_term = pathfindR.data::kegg_genes) {
  message("## Annotating involved genes and visualizing enriched terms")
  ### Argument checks
  if (!is.data.frame(result_df)) {
    stop("`result_df` should be a data frame")
  }
  if (!"ID" %in% colnames(result_df)) {
    stop("`result_df` should contain an \"ID\" column")
  }

  if (!is.data.frame(input_processed)) {
    stop("`input_processed` should be a data frame")
  }
  if (!all(c("GENE", "CHANGE") %in% colnames(input_processed))) {
    stop("`input_processed` should contain the columns \"GENE\" and \"CHANGE\"")
  }

  if (!is.list(genes_by_term)) {
    stop("`genes_by_term` should be a list of term gene sets")
  }
  if (is.null(names(genes_by_term))) {
    stop("`genes_by_term` should be a named list (names are gene set IDs)")
  }

  ### Annotate up/down-regulated term-related genes Up/Down-regulated genes
  upreg <- base::toupper(input_processed$GENE[input_processed$CHANGE >= 0])
  downreg <- base::toupper(input_processed$GENE[input_processed$CHANGE < 0])

  ## Annotation
  annotated_df <- result_df
  annotated_df$Down_regulated <- annotated_df$Up_regulated <- NA
  for (i in base::seq_len(nrow(annotated_df))) {
    idx <- which(names(genes_by_term) == annotated_df$ID[i])
    temp <- genes_by_term[[idx]]
    annotated_df$Up_regulated[i] <- paste(temp[base::toupper(temp) %in% upreg],
      collapse = ", "
    )
    annotated_df$Down_regulated[i] <- paste(temp[base::toupper(temp) %in% downreg],
      collapse = ", "
    )
  }

  return(annotated_df)
}

#' Create HTML Report of pathfindR Results
#'
#' @inheritParams run_pathfindR
#' @param input_processed processed input data frame
#' @param final_res final pathfindR result data frame
#' @param dir_for_report directory to render the report in
create_HTML_report <- function(input, input_processed, final_res, dir_for_report) {
  message("## Creating HTML report")
  rmarkdown::render(
    input = system.file("rmd", "results.Rmd", package = "pathfindR"),
    output_dir = dir_for_report
  )
  rmarkdown::render(
    input = system.file("rmd", "enriched_terms.Rmd", package = "pathfindR"),
    params = list(df = final_res), output_dir = dir_for_report
  )
  rmarkdown::render(
    input = system.file("rmd", "conversion_table.Rmd", package = "pathfindR"),
    params = list(df = input_processed, original_df = input), output_dir = dir_for_report
  )
}
