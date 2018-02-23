#' Example Input for the pathfindR Enrichment Workflow - Rheumatoid Arthritis
#'
#' A dataset containing differentially-expressed genes along with the associated
#' log2-fold-change values and adjusted p-values for the GEO data set GSE15573.
#' This dataset aimed to characterize gene expression profiles in the peripheral
#' blood mononuclear cells of 18 rheumatoid arthritis (RA) patients versus 15
#' healthy subjects. Differentially-expressed genes with adj.P.Val <= 0.05 are
#' presented in this dataset.
#'
#' @format A data frame with 571 rows and 3 variables:
#' \describe{
#'   \item{Gene.symbol}{HGNC gene symbols of the differentially-expressed genes}
#'   \item{logFC}{log2-fold-change values}
#'   \item{adj.P.Val}{adjusted p values, via the Benjamini & Hochberg (1995) method}
#' }
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE15573}
#' @seealso \code{\link{RA_output}} for example output.
"RA_input"

#' Example Output for the pathfindR Enrichment Workflow - Rheumatoid Arthritis
#'
#' A dataset containing the results of pathfindR's active-subnetwork-oriented
#' pathway enrichment workflow performed on the rheumatoid arthritis
#' differential-expression dataset \code{RA_input}. Active subnetwork search
#' was performed with Greedy Algorithm using the Biogrid PIN.
#'
#' @format A data frame with 51 rows and 7 columns:
#' \describe{
#'   \item{ID}{KEGG ID of the enriched pathway}
#'   \item{Pathway}{Description of the enriched pathway}
#'   \item{occurrence}{the number of iterations that the given pathway was found to enriched over all iterations}
#'   \item{lowest_p}{the lowest adjusted-p value of the given pathway over all iterations}
#'   \item{highest_p}{the highest adjusted-p value of the given pathway over all iterations}
#'   \item{Up_regulated}{the up-regulated genes in the input involved in the given pathway, comma-separated}
#'   \item{Down_regulated}{the down-regulated genes in the input involved in the given pathway, comma-separated}
#' }
#' @seealso \code{\link{RA_input}} for example input.
"RA_output"
