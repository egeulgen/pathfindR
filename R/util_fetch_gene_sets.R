#' Fetch Gene Set Objects
#'
#' Function for obtaining the gene sets per term and the term descriptions to
#' be used for enrichment analysis.
#'
#' @param gene_sets Name of the gene sets to be used for enrichment analysis.
#'  Available gene sets are 'KEGG', 'Reactome', 'BioCarta', 'GO-All',
#'  'GO-BP', 'GO-CC', 'GO-MF', 'cell_markers', 'mmu_KEGG' or 'Custom'.
#'  If 'Custom', the arguments \code{custom_genes} and \code{custom_descriptions}
#'  must be specified. (Default = 'KEGG')
#' @param min_gset_size minimum number of genes a term must contain (default = 10)
#' @param max_gset_size maximum number of genes a term must contain (default = 300)
#' @param custom_genes a list containing the genes involved in each custom
#'  term. Each element is a vector of gene symbols located in the given custom
#'  term. Names should correspond to the IDs of the custom terms.
#' @param custom_descriptions A vector containing the descriptions for each
#'  custom  term. Names of the vector should correspond to the IDs of the custom
#'  terms.
#'
#' @return a list containing 2 elements \describe{
#'   \item{genes_by_term}{list of vectors of genes contained in each term}
#'   \item{term_descriptions}{vector of descriptions per each term}
#' }
#'
#' @export
#'
#' @examples
#' KEGG_gset <- fetch_gene_sets()
#' GO_MF_gset <- fetch_gene_sets("GO-MF", min_gset_size = 20, max_gset_size = 100)
fetch_gene_sets <- function(gene_sets = "KEGG", min_gset_size = 10, max_gset_size = 300,
                            custom_genes = NULL, custom_descriptions = NULL) {
  ### Argument checks
  all_gs_opts <- c(
    "KEGG", "Reactome", "BioCarta", "GO-All", "GO-BP", "GO-CC",
    "GO-MF", "cell_markers", "mmu_KEGG", "Custom"
  )
  if (!gene_sets %in% all_gs_opts) {
    stop("`gene_sets` should be one of ", paste(dQuote(all_gs_opts), collapse = ", "))
  }

  if (!is.numeric(min_gset_size)) {
    stop("`min_gset_size` should be numeric")
  }
  if (!is.numeric(max_gset_size)) {
    stop("`max_gset_size` should be numeric")
  }


  ### Custom Gene Sets
  if (gene_sets == "Custom") {
    if (is.null(custom_genes) | is.null(custom_descriptions)) {
      stop("`custom_genes` and `custom_descriptions` must be provided if `gene_sets = \"Custom\"`")
    }

    if (!is.list(custom_genes)) {
      stop("`custom_genes` should be a list of term gene sets")
    }
    if (is.null(names(custom_genes))) {
      stop("`custom_genes` should be a named list (names are gene set IDs)")
    }

    if (!is.atomic(custom_descriptions)) {
      stop("`custom_descriptions` should be a vector of term gene descriptions")
    }
    if (is.null(names(custom_descriptions))) {
      stop("`custom_descriptions` should be a named vector (names are gene set IDs)")
    }

    # filter by size
    gset_lens <- vapply(custom_genes, length, 1)
    keep <- which(gset_lens >= min_gset_size & gset_lens <= max_gset_size)
    custom_genes <- custom_genes[keep]
    custom_descriptions <- custom_descriptions[names(custom_genes)]

    return(list(genes_by_term = custom_genes, term_descriptions = custom_descriptions))
  }

  ### Built-in Gene Sets GO gene sets
  if (grepl("^GO", gene_sets)) {
    genes_by_term <- pathfindR.data::go_all_genes

    GO_df <- pathfindR.data:::GO_all_terms_df
    term_descriptions <- GO_df$GO_term
    names(term_descriptions) <- GO_df$GO_ID

    if (gene_sets == "GO-BP") {
      tmp <- GO_df$GO_ID[GO_df$Category == "Process"]
      genes_by_term <- genes_by_term[tmp]
      term_descriptions <- term_descriptions[tmp]
    } else if (gene_sets == "GO-CC") {
      tmp <- GO_df$GO_ID[GO_df$Category == "Component"]
      genes_by_term <- genes_by_term[tmp]
      term_descriptions <- term_descriptions[tmp]
    } else if (gene_sets == "GO-MF") {
      tmp <- GO_df$GO_ID[GO_df$Category == "Function"]
      genes_by_term <- genes_by_term[tmp]
      term_descriptions <- term_descriptions[tmp]
    }

    ## non-GO (KEGG, Reactome, BioCarta, mmu_KEGG)
  } else {
    if (gene_sets == "KEGG") {
      genes_by_term <- pathfindR.data::kegg_genes
      term_descriptions <- pathfindR.data::kegg_descriptions
    } else if (gene_sets == "Reactome") {
      genes_by_term <- pathfindR.data::reactome_genes
      term_descriptions <- pathfindR.data::reactome_descriptions
    } else if (gene_sets == "BioCarta") {
      genes_by_term <- pathfindR.data::biocarta_genes
      term_descriptions <- pathfindR.data::biocarta_descriptions
    } else if (gene_sets == "mmu_KEGG") {
      genes_by_term <- pathfindR.data::mmu_kegg_genes
      term_descriptions <- pathfindR.data::mmu_kegg_descriptions
    } else {
      genes_by_term <- pathfindR.data::cell_markers_gsets
      term_descriptions <- pathfindR.data::cell_markers_descriptions
    }
  }

  # filter by size
  term_lens <- vapply(genes_by_term, length, 1)
  keep <- which(term_lens >= min_gset_size & term_lens <= max_gset_size)
  genes_by_term <- genes_by_term[keep]
  term_descriptions <- term_descriptions[names(genes_by_term)]

  return(list(genes_by_term = genes_by_term, term_descriptions = term_descriptions))
}
