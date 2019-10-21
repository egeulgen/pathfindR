#' Example Input for pathfindR - Enriched Term Scoring
#'
#' A matrix containing the log2-normalized expression values of the differentially-expressed genes
#' for 18 rheumatoid arthritis (RA) patients and 15 healthy subjects. Expression values of 572
#' differentially-expressed genes with adj.P.Val <= 0.05 are presented in this dataset.
#' \emph{Generated on Sep 28, 2019.}
#'
#' @format A matrix with 572 rows and 33 columns.
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE15573}
"RA_exp_mat"

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
#' enrichment workflow performed on the rheumatoid arthritis
#' differential-expression result dataset \code{RA_input}. Analysis via
#' \code{run_pathfindR} was performed using the default settings.
#' \emph{Generated on Oct 22, 2019.}
#'
#' @format A data frame with 101 rows and 8 columns:
#' \describe{
#'   \item{ID}{ID of the enriched term}
#'   \item{Term_Description}{Description of the enriched term}
#'   \item{Fold_Enrichment}{Fold enrichment value for the enriched term}
#'   \item{occurrence}{the number of iterations that the given term was found to enriched over all iterations}
#'   \item{lowest_p}{the lowest adjusted-p value of the given term over all iterations}
#'   \item{highest_p}{the highest adjusted-p value of the given term over all iterations}
#'   \item{Up_regulated}{the up-regulated genes in the input involved in the given term, comma-separated}
#'   \item{Down_regulated}{the down-regulated genes in the input involved in the given term, comma-separated}
#' }
#' @seealso \code{\link{RA_input}} for example input of the enrichment workflow.
#' \code{\link{RA_clustered}} for example output of the clustering workflow.
#' \code{\link{run_pathfindR}} for details on the patfindR enrichment workflow.
"RA_output"

#' Example Output for the pathfindR Clustering Workflow - Rheumatoid Arthritis
#'
#' A dataset containing the results of pathfindR's clustering and
#' partitioning  workflow performed on the rheumatoid arthritis
#' enrichment results \code{RA_output}. The clustering and partitioning
#' function \code{cluster_enriched_terms} was used with the default settings
#' (i.e. hierarchical clustering was performed and the agglomeration method
#' was "average"). The optimal number of clusters (yielding the highest average
#' silhouette width) was determined to be 16. Finally, the enriched terms with the
#' lowest p values in each cluster were assigned as representative terms
#' for that cluster.
#' \emph{Generated on Oct 22, 2019.}
#'
#' @format A data frame with 101 rows and 10 columns:
#' \describe{
#'   \item{ID}{ID of the enriched term}
#'   \item{Term_Description}{Description of the enriched term}
#'   \item{Fold_Enrichment}{Fold enrichment value for the enriched term}
#'   \item{occurrence}{the number of iterations that the given term was found to enriched over all iterations}
#'   \item{lowest_p}{the lowest adjusted-p value of the given term over all iterations}
#'   \item{highest_p}{the highest adjusted-p value of the given term over all iterations}
#'   \item{Up_regulated}{the up-regulated genes in the input involved in the given term, comma-separated}
#'   \item{Down_regulated}{the down-regulated genes in the input involved in the given term, comma-separated}
#'   \item{Cluster}{the cluster to which the enriched term is assigned}
#'   \item{Status}{whether the enriched term is the "Representative" term in its cluster or only a "Member"}
#' }
#' @seealso \code{\link{RA_input}} for example input of the enrichment workflow.
#' \code{\link{RA_output}} for example output of the enrichment workflow.
#' \code{\link{cluster_enriched_terms}} for details on the patfindR clustering approaches.
"RA_clustered"

#' KEGG Pathways - Gene Sets
#'
#' A list containing the genes involved in each human KEGG pathway. Each element
#' is a vector of gene symbols located in the given pathway. Names correspond to
#' the KEGG ID of the pathway. Pathways that did not contain any genes were
#' discarded.
#' \emph{Generated on Oct 12, 2019.}
#'
#' @format list containing 327 vectors of gene symbols. Each vector corresponds
#'   to a pathway.
"kegg_genes"

#' KEGG Pathways - Descriptions
#'
#' A list containing the descriptions for each human KEGG pathway. Names of the
#' list correspond to the KEGG ID of the pathway. Pathways that did not contain
#' any genes were discarded.
#' \emph{Generated on Oct 12, 2019.}
#'
#' @format vector containing 327 character values, the descriptions for the given
#'   pathways.
"kegg_descriptions"

#' Reactome Pathways - Gene Sets
#'
#' A list containing the genes involved in each human Reactome pathway. Each
#' element is a vector of gene symbols located in the given pathway. Names
#' correspond to the Reactome ID of the pathway.
#' \emph{Generated on Oct 12, 2019.}
#'
#' @format list containing 2263 vectors of gene symbols. Each vector corresponds
#'   to a pathway.
"reactome_genes"

#' Reactome Pathways - Descriptions
#'
#' A list containing the descriptions for each human Reactome pathway. Names of the
#' list correspond to the Reactome ID of the pathway.
#' \emph{Generated on Oct 12, 2019.}
#'
#' @format list containing 2263 character values, the descriptions for the given
#'   pathways.
"reactome_descriptions"

#' BioCarta Pathways - Gene Sets
#'
#' A list containing the genes involved in each human BioCarta pathway. Each
#' element is a vector of gene symbols located in the given pathway.
#' \emph{Generated on Oct 12, 2019.}
#'
#' @format list containing 314 vectors of gene symbols. Each vector corresponds
#'   to a gene set.
"biocarta_genes"

#' BioCarta Pathways - Descriptions
#'
#' A list containing the descriptions for each human BioCarta pathway.
#' \emph{Generated on Oct 12, 2019.}
#'
#' @format list containing 314 character values, the descriptions for the given
#'   pathways.
"biocarta_descriptions"

#' Gene Ontology - All Gene Ontology Gene Sets
#'
#' A list containing the genes involved in each GO ontology term. Each
#' element is a vector of gene symbols located in the given gene set.
#' \emph{Generated on Sep 30, 2019.}
#'
#' @format list containing 14586 vectors of gene symbols. Each vector corresponds
#'   to a gene set.
"go_all_genes"

#' Gene Ontology - All Gene Ontology Descriptions
#'
#' A data frame containing descriptions of Gene Ontology terms (for all categories)
#' \emph{Generated on Sep 30, 2019.}
#'
#' @format data frame containing 14586 rows and 3 columns. Columns are \describe{
#' \item{GO_ID}{ID of the GO term}
#' \item{GO_term}{Description the GO term}
#' \item{Category}{Category of the GO term (i.e., "Component", "Function" or "Process")}
#' }
"GO_all_terms_df"



#' Custom Gene Set Enrichment Results
#'
#' A data frame consisting of pathfindR enrichment analysis results on the
#' example TF target genes data (target gene sets of CREB and MYC).
#' \emph{Generated on Oct 22, 2019.}
#' @format data frame containing 2 rows and 8 columns. Each row is a gene set (the TF target gene sets).
"custom_result"

#' Example Active Subnetworks
#'
#' A list of vectors containing genes for each active subnetwork that passed
#' the filtering step.
#' \emph{Generated on Oct 22, 2019.}
#'
#' @format list containing 175 vectors. Each vector is the set of genes for the
#' given active subnetwork.
"example_active_snws"

#' BioGRID PIN Adjacency List
#'
#' An adjacency list of vectors containing interactors B for each interactor A
#' in the BioGRID protein-protein interaction network (The designations
#' "interactor A" and "interactor B" are arbitrary).
#' \emph{Generated on Oct 12, 2019.}
#'
#' @format list containing 15231 vectors. Each vector is the set of gene symbols
#' of interactors B for each interactor A.
"biogrid_adj_list"

#' GeneMania PIN Adjacency List
#'
#' An adjacency list of vectors containing interactors B for each interactor A
#' in the GeneMania protein-protein interaction network (The designations
#' "interactor A" and "interactor B" are arbitrary).
#' \emph{Generated on Oct 12, 2019.}
#'
#' @format list containing 12345 vectors. Each vector is the set of gene symbols
#' of interactors B for each interactor A.
"gene_mania_adj_list"

#' KEGG PIN Adjacency List
#'
#' An adjacency list of vectors containing interactors B for each interactor A
#' in the KEGG protein-protein interaction network (The designations
#' "interactor A" and "interactor B" are arbitrary).
#' \emph{Generated on Oct 12, 2019.}
#'
#' @format list containing 4507 vectors. Each vector is the set of gene symbols
#' of interactors B for each interactor A.
"kegg_adj_list"

#' IntAct PIN Adjacency List
#'
#' An adjacency list of vectors containing interactors B for each interactor A
#' in the IntAct protein-protein interaction network (The designations
#' "interactor A" and "interactor B" are arbitrary).
#' \emph{Generated on Oct 19, 2019.}
#'
#' @format list containing 15057 vectors. Each vector is the set of gene symbols
#' of interactors B for each interactor A.
"intact_adj_list"
