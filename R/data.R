#' Example Input for the pathfindr Enrichment Workflow
#'
#' A dataset containing differentially-expressed genes along with the associated
#' log2-fold-changes and adjusted p-values for the GEO data set GSE50772. This dataset
#' aimed to characterize gene expression profiles in the peripheral blood mononuclear
#' cells of SLE patients versus healthy blood donors.
#'
#' @format A data frame with 572 rows and 3 variables:
#' \describe{
#'   \item{Gene.symbol}{HGNC gene symbols of the differentially-expressed genes}
#'   \item{logFC}{log2-fold-change values}
#'   \item{adj.P.Val}{adjusted p values, via the Benjamini & Hochberg(1995) method}
#' }
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50772}
"example_input"
