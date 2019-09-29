#' Example Input for the pathfindR Enrichment Workflow - Rheumatoid Arthritis
#'
#' A dataset containing the differentially-expressed genes along with the
#' associated log2-fold-change values and adjusted p-values for the GEO dataset
#' GSE15573. The microarray dataset aimed to characterize gene expression profiles in the
#' peripheral blood mononuclear cells of 18 rheumatoid arthritis (RA) patients
#' versus 15 healthy subjects. Differentially-expressed genes with adj.P.Val <
#' 0.05 are presented in this dataset.
#' \emph{Generated on Sep 28, 2019.}
#'
#' @format A data frame with 572 rows and 3 variables: \describe{
#'   \item{Gene.symbol}{HGNC gene symbols of the differentially-expressed genes}
#'   \item{logFC}{log2-fold-change values}
#'   \item{adj.P.Val}{adjusted p values, via the Benjamini & Hochberg (1995) method}
#' }
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE15573}
#'
#' @seealso \code{\link{RA_output}} for example output of the enrichment workflow.
#' \code{\link{RA_clustered}} for example output of the clustering workflow.
#' \code{\link{run_pathfindR}} for details on the patfindR enrichment analysis.
"RA_input"

#' Example Output for the pathfindR Enrichment Workflow - Rheumatoid Arthritis
#'
#' A dataset containing the results of pathfindR's active-subnetwork-oriented
#' pathway enrichment workflow performed on the rheumatoid arthritis
#' differential-expression result dataset \code{RA_input}. Analysis via
#' \code{run_pathfindR} was performed using the default settings.
#' \emph{Generated on Sep 28, 2019.}
#'
#' @format A data frame with 89 rows and 8 columns:
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
#' \code{\link{run_pathfindR}} for details on the patfindR enrichment analysis.
"RA_output"

#' Example Output for the pathfindR Clustering Workflow - Rheumatoid Arthritis
#'
#' A dataset containing the results of pathfindR's pathway clustering and
#' partitioning  workflow performed on the rheumatoid arthritis
#' enrichment results \code{RA_output}. The clustering and partitioning
#' function \code{cluster_pathways} was used with the default settings
#' (i.e. hierarchical clustering was performed and the agglomeration method
#' was "average"). The optimal number of clusters (yielding the highest average
#' silhouette width) was determined to be 22. Finally, the pathways with the
#' lowest p values in each cluster were assigned as representative pathways
#' for that cluster.
#' \emph{Generated on Sep 28, 2019.}
#'
#' @format A data frame with 89 rows and 10 columns:
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
#' \code{\link{cluster_pathways}} for details on the patfindR clustering approaches.
"RA_clustered"

#' Example Input for pathfindR - pathway z-scores
#'
#' A matrix containing the log2-normalized expression values of the differentially-expressed genes
#' for 18 rheumatoid arthritis (RA) patients and 15 healthy subjects. Expression values of
#' differentially-expressed genes with adj.P.Val <= 0.05 are presented in this dataset.
#' \emph{Generated on Sep 28, 2019.}
#'
#' @format A matrix with 572 rows and 33 columns.
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE15573}
"RA_exp_mat"



#' KEGG Gene Sets
#'
#' A list containing the genes involved in each human KEGG pathway. Each element
#' is a vector of gene symbols located in the given pathway. Names correspond to
#' the KEGG ID of the pathway. Pathways that did not contain any genes were
#' discarded.
#' \emph{Generated on Sep 30, 2019.}
#'
#' @format list containing 326 vectors of gene symbols. Each vector corresponds
#'   to a pathway.
"kegg_genes"

#' KEGG Pathway Descriptions
#'
#' A list containing the descriptions for each human KEGG pathway. Names of the
#' list correspond to the KEGG ID of the pathway. Pathways that did not contain
#' any genes were discarded.
#' \emph{Generated on Sep 30, 2019.}
#'
#' @format vector containing 326 character values, the descriptions for the given
#'   pathways.
"kegg_descriptions"

#' Reactome Gene Sets
#'
#' A list containing the genes involved in each human Reactome pathway. Each
#' element is a vector of gene symbols located in the given pathway. Names
#' correspond to the Reactome ID of the pathway.
#' \emph{Generated on Sep 30, 2019.}
#'
#' @format list containing 2263 vectors of gene symbols. Each vector corresponds
#'   to a pathway.
"reactome_genes"

#' Reactome Pathway Descriptions
#'
#' A list containing the descriptions for each human Reactome pathway. Names of the
#' list correspond to the Reactome ID of the pathway.
#' \emph{Generated on Sep 30, 2019.}
#'
#' @format list containing 2263 character values, the descriptions for the given
#'   pathways.
"reactome_descriptions"

#' BioCarta Gene Sets
#'
#' A list containing the genes involved in each human BioCarta pathway. Each
#' element is a vector of gene symbols located in the given pathway.
#' \emph{Generated on Sep 30, 2019.}
#'
#' @format list containing 314 vectors of gene symbols. Each vector corresponds
#'   to a gene set.
"biocarta_genes"

#' BioCarta Pathway Descriptions
#'
#' A list containing the descriptions for each human BioCarta pathway.
#' \emph{Generated on Sep 30, 2019.}
#'
#' @format list containing 314 character values, the descriptions for the given
#'   pathways.
"biocarta_descriptions"

#' Gene Ontology - All Ontology Gene Sets
#'
#' A list containing the genes involved in each GO ontology term. Each
#' element is a vector of gene symbols located in the given gene set.
#' \emph{Generated on Sep 30, 2019.}
#'
#' @format list containing 14586 vectors of gene symbols. Each vector corresponds
#'   to a gene set.
"go_all_genes"

#' Gene Ontology - All Ontology Descriptions
#'
#' A list containing the descriptions for each human GO ontology term.
#' \emph{Generated on Sep 30, 2019.}
#'
#' @format list containing 14586 character values, the descriptions for the given
#' go gene set.
"go_all_descriptions"



#' Custom Gene Set Enrichment Results
#'
#' A data frame consisting of pathfindR enrichment results on the example TF target data.
#' \emph{Generated on Sep 28, 2019.}
#' @format data frame containing 2 rows and 8 columns. Each row is a gene set (the TF target gene sets).
"custom_result"

#' Example Active Subnetworks
#'
#' A list of vectors containing genes for each active subnetwork that passed
#' the filtering step.
#' \emph{Generated on Sep 28, 2019.}
#'
#' @format list containing 167 vectors. Each vector is the set of genes for the
#' given active subnetwork.
"example_active_snws"
