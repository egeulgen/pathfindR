#' Example Input for the pathfindR Enrichment Workflow - Rheumatoid Arthritis
#'
#' A dataset containing the differentially-expressed genes along with the
#' associated log2-fold-change values and adjusted p-values for the GEO dataset
#' GSE15573. The microarray dataset aimed to characterize gene expression profiles in the
#' peripheral blood mononuclear cells of 18 rheumatoid arthritis (RA) patients
#' versus 15 healthy subjects. Differentially-expressed genes with adj.P.Val <
#' 0.05 are presented in this dataset.
#'
#' @format A data frame with 572 rows and 3 variables: \describe{
#'   \item{Gene.symbol}{HGNC gene symbols of the differentially-expressed genes}
#'   \item{logFC}{log2-fold-change values}
#'   \item{adj.P.Val}{adjusted p values, via the Benjamini & Hochberg (1995) method}
#' }
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE15573}
#' @seealso \code{\link{RA_output}} for example output of the enrichment workflow.
#' \code{\link{RA_clustered}} for example output of the clustering workflow.
"RA_input"

#' Example Output for the pathfindR Enrichment Workflow - Rheumatoid Arthritis
#'
#' A dataset containing the results of pathfindR's active-subnetwork-oriented
#' pathway enrichment workflow performed on the rheumatoid arthritis
#' differential-expression dataset \code{RA_input}. Active subnetwork search
#' was performed with Greedy Algorithm using the Biogrid PIN.
#'
#' @format A data frame with 33 rows and 8 columns:
#' \describe{
#'   \item{ID}{KEGG ID of the enriched pathway}
#'   \item{Pathway}{Description of the enriched pathway}
#'   \item{Fold_Enrichment}{Fold enrichment value for the enriched pathway}
#'   \item{occurrence}{the number of iterations that the given pathway was found to enriched over all iterations}
#'   \item{lowest_p}{the lowest adjusted-p value of the given pathway over all iterations}
#'   \item{highest_p}{the highest adjusted-p value of the given pathway over all iterations}
#'   \item{Up_regulated}{the up-regulated genes in the input involved in the given pathway, comma-separated}
#'   \item{Down_regulated}{the down-regulated genes in the input involved in the given pathway, comma-separated}
#' }
#' @seealso \code{\link{RA_input}} for example input of the enrichment workflow.
#' \code{\link{RA_clustered}} for example output of the clustering workflow.
"RA_output"

#' Example Output for the pathfindR Clustering Workflow - Rheumatoid Arthritis
#'
#' A dataset containing the results of pathfindR's pathway clustering and
#' partitioning  workflow performed on the rheumatoid arthritis
#' enrichment results \code{RA_output}. The number of clusters were detected
#' automatically as 11 and the agglomeration method was "average".
#'
#' @format A data frame with 33 rows and 8 columns:
#' \describe{
#'   \item{ID}{KEGG ID of the enriched pathway}
#'   \item{Pathway}{Description of the enriched pathway}
#'   \item{Fold_Enrichment}{Fold enrichment value for the enriched pathway}
#'   \item{occurrence}{the number of iterations that the given pathway was found to enriched over all iterations}
#'   \item{lowest_p}{the lowest adjusted-p value of the given pathway over all iterations}
#'   \item{highest_p}{the highest adjusted-p value of the given pathway over all iterations}
#'   \item{Up_regulated}{the up-regulated genes in the input involved in the given pathway, comma-separated}
#'   \item{Down_regulated}{the down-regulated genes in the input involved in the given pathway, comma-separated}
#'   \item{Cluster}{the cluster to which the pathway is assigned}
#'   \item{Status}{whether the pathway is the "Representative" pathway in its cluster or only a "Member"}
#' }
#' @seealso \code{\link{RA_input}} for example input of the enrichment workflow.
#' \code{\link{RA_output}} for example output of the enrichment workflow.
"RA_clustered"

#' Example Input for pathfindR - pathway z-scores
#'
#' A matrix containing the log2-normalized expression values of the differentially-expressed genes
#' for 18 rheumatoid arthritis (RA) patients and 15 healthy subjects. Expression values of
#' differentially-expressed genes with adj.P.Val <= 0.05 are presented in this dataset.
#'
#' @format A matrix with 572 rows and 33 columns.
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE15573}
"RA_exp_mat"



#' KEGG Gene Sets
#'
#' A list containing the genes involved in each human KEGG pathway. Each element
#' is a vector of gene symbols located in the given pathway. Names correspond to
#' the KEGG ID of the pathway. Pathways that did not contain any genes were
#' discarded. This data was retrieved on May 13, 2018.
#'
#' @format list containing 319 vectors of gene symbols. Each vector corresponds
#'   to a pathway.
"kegg_genes"

#' KEGG Pathway Descriptions
#'
#' A list containing the descriptions for each human KEGG pathway. Names of the
#' list correspond to the KEGG ID of the pathway. Pathways that did not contain
#' any genes were discarded. This data was retrieved on May 13, 2018.
#'
#' @format list containing 319 character values, the descriptions for the given
#'   pathways.
"kegg_pathways"

#' Reactome Gene Sets
#'
#' A list containing the genes involved in each human Reactome pathway. Each
#' element is a vector of gene symbols located in the given pathway. Names
#' correspond to the Reactome ID of the pathway. This data was retrieved on May 13,
#' 2018.
#'
#' @format list containing 2022 vectors of gene symbols. Each vector corresponds
#'   to a pathway.
"reactome_genes"

#' Reactome Pathway Descriptions
#'
#' A list containing the descriptions for each human Reactome pathway. Names of the
#' list correspond to the Reactome ID of the pathway. This data was retrieved on May
#' 13, 2018.
#'
#' @format list containing 2022 character values, the descriptions for the given
#'   pathways.
"reactome_pathways"

#' BioCarta Gene Sets
#'
#' A list containing the genes involved in each human BioCarta pathway. Each
#' element is a vector of gene symbols located in the given pathway. This data
#' was retrieved on May 13, 2018.
#'
#' @format list containing 217 vectors of gene symbols. Each vector corresponds
#'   to a gene set.
"biocarta_genes"

#' BioCarta Pathway Descriptions
#'
#' A list containing the descriptions for each human Reactome pathway. This data
#' was retrieved on May 13, 2018.
#'
#' @format list containing 217 character values, the descriptions for the given
#'   pathways.
"biocarta_pathways"

#' Gene Ontology - Biological Process Ontology Gene Sets
#'
#' A list containing the genes involved in each GO Biological Process. Each
#' element is a vector of gene symbols located in the given gene set. This data
#' was retrieved on May 13, 2018.
#'
#' @format list containing 3941 vectors of gene symbols. Each vector corresponds
#'   to a gene set.
"go_bp_genes"

#' Gene Ontology - Biological Process Ontology Descriptions
#'
#' A list containing the descriptions for each human GO Biological Process. This
#' data was retrieved on May 13, 2018.
#'
#' @format list containing 3941 character values, the descriptions for the given
#'   pathways.
"go_bp_pathways"

#' Gene Ontology - Cellular Component Ontology Gene Sets
#'
#' A list containing the genes involved in each GO Cellular Component. Each
#' element is a vector of gene symbols located in the given gene set. This data
#' was retrieved on May 13, 2018.
#'
#' @format list containing 470 vectors of gene symbols. Each vector corresponds
#'   to a gene set.
"go_cc_genes"

#' Gene Ontology - Cellular Component Ontology Descriptions
#'
#' A list containing the descriptions for each human GO Cellular Component. This
#' data was retrieved on May 13, 2018.
#'
#' @format list containing 470 character values, the descriptions for the given
#'   pathways.
"go_cc_pathways"

#' Gene Ontology - Molecular Function Ontology Gene Sets
#'
#' A list containing the genes involved in each GO Molecular Function. Each
#' element is a vector of gene symbols located in the given gene set. This data
#' was retrieved on May 13, 2018.
#'
#' @format list containing 713 vectors of gene symbols. Each vector corresponds
#'   to a gene set.
"go_mf_genes"

#' Gene Ontology - Molecular Function Ontology Descriptions
#'
#' A list containing the descriptions for each human GO Molecular Function. This
#' data was retrieved on May 13, 2018.
#'
#' @format list containing 713 character values, the descriptions for the given
#'   pathways.
"go_mf_pathways"

#' Custom Gene Set Enrichment Results
#'
#' A data frame consisting of pathfindR enrichment results on the example TF target data.
#'
#' @format data frame containing 2 rows and 8 columns. Each row is a gene set (the TF target gene sets).
"custom_result"
