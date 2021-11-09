#' Process Data frame of Protein-protein Interactions
#'
#' @param pin_df data frame of protein-protein interactions with 2 columns:
#' "Interactor_A" and "Interactor_B"
#'
#' @return processed PIN data frame (removes self-interactions and
#' duplicated interactions)
process_pin <- function(pin_df) {
  # remove self-interactions
  pin_df <- pin_df[pin_df$Interactor_A != pin_df$Interactor_B, ]

  # remove duplicated inteactions (including symmetric ones)
  pin_df <- unique(t(apply(pin_df, 1, sort)))

  pin_df <- as.data.frame(pin_df)
  colnames(pin_df) <- c("Interactor_A", "Interactor_B")
  return(pin_df)
}

#' Retrieve the Requested Release of Organism-specific BioGRID PIN
#'
#' @param org organism name. BioGRID naming requires underscores for spaces so
#' "Homo sapiens" becomes "Homo_sapiens", "Mus musculus" becomes "Mus_musculus"
#' etc. See \url{https://wiki.thebiogrid.org/doku.php/statistics} for a full
#' list of available organisms (default = "Homo_sapiens")
#' @param path2pin the path of the file to save the PIN data. By default, the
#' PIN data is saved in a temporary file
#' @param release the requested BioGRID release (default = "4.4.200")
#'
#' @return the path of the file in which the PIN data was saved. If
#' \code{path2pin} was not supplied by the user, the PIN data is saved in a
#' temporary file
get_biogrid_pin <- function(org = "Homo_sapiens", path2pin, release = "4.4.200") {
  # check organism name
  all_org_names <- c("Anopheles_gambiae_PEST", "Apis_mellifera",
                     "Arabidopsis_thaliana_Columbia", "Bacillus_subtilis_168",
                     "Bos_taurus", "Caenorhabditis_elegans",
                     "Candida_albicans_SC5314", "Canis_familiaris",
                     "Cavia_porcellus", "Chlamydomonas_reinhardtii",
                     "Chlorocebus_sabaeus", "Cricetulus_griseus",
                     "Danio_rerio", "Dictyostelium_discoideum_AX4",
                     "Drosophila_melanogaster", "Emericella_nidulans_FGSC_A4",
                     "Equus_caballus", "Escherichia_coli_K12_MC4100_BW2952",
                     "Escherichia_coli_K12_MG1655", "Escherichia_coli_K12_W3110",
                     "Escherichia_coli_K12", "Gallus_gallus", "Glycine_max",
                     "Hepatitus_C_Virus", "Homo_sapiens", "Human_Herpesvirus_1",
                     "Human_Herpesvirus_2", "Human_Herpesvirus_3",
                     "Human_Herpesvirus_4", "Human_Herpesvirus_5",
                     "Human_Herpesvirus_6A", "Human_Herpesvirus_6B",
                     "Human_Herpesvirus_7", "Human_Herpesvirus_8",
                     "Human_Immunodeficiency_Virus_1",
                     "Human_Immunodeficiency_Virus_2", "Human_papillomavirus_10",
                     "Human_papillomavirus_16", "Human_papillomavirus_6b",
                     "Leishmania_major_Friedlin", "Macaca_mulatta",
                     "Meleagris_gallopavo", "Mus_musculus",
                     "Mycobacterium_tuberculosis_H37Rv",
                     "Neurospora_crassa_OR74A", "Nicotiana_tomentosiformis",
                     "Oryctolagus_cuniculus", "Oryza_sativa_Japonica",
                     "Ovis_aries", "Pan_troglodytes", "Pediculus_humanus",
                     "Plasmodium_falciparum_3D7", "Rattus_norvegicus",
                     "Ricinus_communis", "Saccharomyces_cerevisiae_S288c",
                     "Schizosaccharomyces_pombe_972h",
                     "Selaginella_moellendorffii",
                     "Simian_Immunodeficiency_Virus", "Simian_Virus_40",
                     "Solanum_lycopersicum", "Solanum_tuberosum",
                     "Streptococcus_pneumoniae_ATCCBAA255",
                     "Strongylocentrotus_purpuratus", "Sus_scrofa",
                     "Tobacco_Mosaic_Virus", "Ustilago_maydis_521",
                     "Vaccinia_Virus", "Vitis_vinifera", "Xenopus_laevis",
                     "Zea_mays")
  if (!org %in% all_org_names)
    stop(paste(org, "is not a valid Biogrid organism.",
               "Available organisms are listed on: https://wiki.thebiogrid.org/doku.php/statistics"))

  # release directory for download
  rel_dir <- paste0("BIOGRID-", release)

  # choose tab2 vs. tab3
  tab_v <- ifelse(utils::compareVersion(release, "3.5.183") == -1, ".tab2", ".tab3")

  # download tab2 format organism files
  tmp <- tempfile()
  fname <- paste0("BIOGRID-ORGANISM-", release, tab_v)
  biogrid_url <- paste0("https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/", rel_dir, "/", fname, ".zip")
  utils::download.file(biogrid_url, tmp, method = getOption("download.file.method"), quiet = TRUE)

  # parse organism names
  all_org_files <- utils::unzip(tmp, list = TRUE)
  all_org_files$Organism <- sub("\\.tab\\d\\.txt", "", all_org_files$Name)
  all_org_files$Organism <- sub("BIOGRID-ORGANISM-", "", all_org_files$Organism)
  all_org_files$Organism <- sub("-.*\\d+$", "", all_org_files$Organism)

  org_file <- all_org_files$Name[all_org_files$Organism == org]

  # process and save organism PIN file
  biogrid_df <- utils::read.delim(unz(tmp, org_file),
                                  check.names = FALSE,
                                  colClasses = "character",
                                  stringsAsFactors = FALSE)
  biogrid_pin <- data.frame(Interactor_A = biogrid_df[, "Official Symbol Interactor A"],
                            Interactor_B = biogrid_df[, "Official Symbol Interactor B"],
                            stringsAsFactors = FALSE)
  biogrid_pin <- process_pin(biogrid_pin)

  final_pin <- data.frame(intA = biogrid_pin$Interactor_A,
                          pp = "pp",
                          intB = biogrid_pin$Interactor_B,
                          stringsAsFactors = FALSE)

  if (missing(path2pin))
    path2pin <- tempfile()
  utils::write.table(final_pin,
                     path2pin,
                     sep = "\t",
                     row.names = FALSE, col.names = FALSE, quote = FALSE)
  return(path2pin)
}

#' Retrieve Organism-specific PIN data
#'
#' @param source As of this version, this function is implemented to get data
#' from "BioGRID" only. This argument (and this wrapper function) was implemented
#' for future utility
#' @inheritParams get_biogrid_pin
#' @param ... additional arguments for \code{\link{get_biogrid_pin}}
#'
#' @return the path of the file in which the PIN data was saved. If
#' \code{path2pin} was not supplied by the user, the PIN data is saved in a
#' temporary file
#' @export
#'
#' @examples
#' \dontrun{
#' pin_path <- get_pin_file()
#' }
get_pin_file <- function(source = "BioGRID", org = "Homo_sapiens", path2pin, ...) {
  ## TODO
  if (source != "BioGRID")
    stop("As of this version, this function is implemented to get data from BioGRID only")

  path2pin <- get_biogrid_pin(org = org, path2pin = path2pin, ...)
  return(path2pin)
}

#' Retrieve Gene Sets from GMT-format File
#'
#' @param path2gmt path to the gmt file
#'
#' @return list containing 2 elements: \itemize{
#' \item{gene_sets}{A list containing the genes involved in each gene set}
#' \item{descriptions}{A named vector containing the descriptions for each gene set}
#' }
gset_list_from_gmt <- function(path2gmt) {
  gmt_lines <- readLines(path2gmt)

  ## Genes list
  genes_list <- lapply(gmt_lines, function(x) {
    x <- unlist(strsplit(x, "\t"))
    x <- unique(x[3:length(x)])
    return(x)
  })

  names(genes_list) <- vapply(gmt_lines, function(x) {
    x <- unlist(strsplit(x, "\t"))
    return(x[2])
  }, "a")

  ## Decriptions vector
  descriptions_vec <- vapply(gmt_lines, function(x) {
    x <- unlist(strsplit(x, "\t"))
    return(x[1])
  }, "a")

  names(descriptions_vec) <- vapply(gmt_lines, function(x) {
    x <- unlist(strsplit(x, "\t"))
    return(x[2])
  }, "a")

  # remove empty gene sets (if any)
  genes_list <- genes_list[vapply(genes_list, length, 1) != 0]
  descriptions_vec <- descriptions_vec[names(genes_list)]

  return(list(gene_sets = genes_list, descriptions = descriptions_vec))
}

#' Retrieve Organism-specific KEGG Pathway Gene Sets
#'
#' @param org_code KEGG organism code for the selected organism. For a full list
#' of all available organisms, see \url{https://www.genome.jp/kegg/catalog/org_list.html}
#'
#' @return list containing 2 elements: \itemize{
#' \item{gene_sets}{A list containing the genes involved in each KEGG pathway}
#' \item{descriptions}{A named vector containing the descriptions for each KEGG pathway}
#' }
get_kegg_gsets <- function(org_code = "hsa") {
  # created named list, eg:  path:map00010: "Glycolysis / Gluconeogenesis"
  pathways_list <- KEGGREST::keggList("pathway", org_code)

  # make them into KEGG-style pathway identifiers
  pathway_codes <- sub("path:", "", names(pathways_list))

  # parse pathway genes
  genes_by_pathway <- lapply(pathway_codes, function(pwid) {
    pw <- KEGGREST::keggGet(pwid)

    ## get gene symbols
    all_entries <- pw[[1]]$GENE
    if(is.null(all_entries))
      return(NULL)
    tmp <- c(TRUE, FALSE)
    if (grepl(";", all_entries[2])) {
      tmp <- c(FALSE, TRUE)
    }
    pw <- all_entries[tmp]

    pw <- sub(";.+", "", pw) ## discard any description
    pw <- pw[grep("^[A-Za-z0-9_-]+(\\@)?$", pw)] ## remove mistaken lines
    pw <- unique(pw) ## keep unique genes
    return(pw)
  })

  names(genes_by_pathway) <- pathway_codes

  # remove empty gene sets (metabolic pathways)
  kegg_genes <- genes_by_pathway[vapply(genes_by_pathway, length, 1) != 0]

  kegg_descriptions <- pathways_list
  names(kegg_descriptions) <- sub("path:", "", names(kegg_descriptions))
  kegg_descriptions <- sub(" & .*$", "", sub("-([^-]*)$", "&\\1", kegg_descriptions))
  kegg_descriptions <- kegg_descriptions[names(kegg_descriptions) %in% names(kegg_genes)]

  result <- list(gene_sets = kegg_genes, descriptions = kegg_descriptions)
  return(result)
}

#' Retrieve Reactome Pathway Gene Sets
#'
#' @return Gets the latest Reactome pathways gene sets in gmt format. Parses the
#' gmt file and returns a list containing 2 elements: \itemize{
#' \item{gene_sets}{A list containing the genes involved in each Reactome pathway}
#' \item{descriptions}{A named vector containing the descriptions for each Reactome pathway}
#' }
#'
get_reactome_gsets <- function() {
  tmp <- tempfile()
  reactome_url <- "https://reactome.org/download/current/ReactomePathways.gmt.zip"
  utils::download.file(reactome_url, tmp, method = getOption("download.file.method"))

  reactome_gmt <- unz(tmp, "ReactomePathways.gmt")
  result <- gset_list_from_gmt(reactome_gmt)
  close(reactome_gmt)

  # fix illegal char(s)
  result$descriptions <- gsub("[^ -~]", "", result$descriptions)
  return(result)
}

#' Retrieve Organism-specific MSigDB Gene Sets
#'
#' @param species species name, such as Homo sapiens, Mus musculus, etc.
#' See \code{\link[msigdbr]{msigdbr_show_species}} for all the species available in
#' the msigdbr package
#' @param collection collection. i.e., H, C1, C2, C3, C4, C5, C6, C7.
#' @param subcollection sub-collection, such as CGP, BP, etc. (default = NULL,
#' i.e. list all gene sets in collection)
#'
#' @return Retrieves the MSigDB gene sets and returns a list containing 2 elements: \itemize{
#' \item{gene_sets}{A list containing the genes involved in each of the selected MSigDB gene sets}
#' \item{descriptions}{A named vector containing the descriptions for each selected MSigDB gene set}
#' }
#'
#' @details this function utilizes the function \code{\link[msigdbr]{msigdbr}}
#' from the \code{msigdbr} package to retrieve the 'Molecular Signatures Database'
#' (MSigDB) gene sets (Subramanian et al. 2005 <doi:10.1073/pnas.0506580102>,
#' Liberzon et al. 2015 <doi:10.1016/j.cels.2015.12.004>).
#' Available collections are: H: hallmark gene sets, C1: positional gene sets,
#' C2: curated gene sets, C3: motif gene sets, C4: computational gene sets,
#' C5: GO gene sets, C6: oncogenic signatures and C7: immunologic signatures
get_mgsigdb_gsets <- function(species = "Homo sapiens", collection, subcollection = NULL) {
  # arg check
  all_collections <- c("H", "C1", "C2", "C3", "C4", "C5", "C6", "C7")
  if (!collection %in% all_collections)
    stop("`collection` should be one of ", paste(dQuote(all_collections), collapse = ", "))

  # retrieve msigdbr df
  msig_df <- msigdbr::msigdbr(species = species,
                              category = collection,
                              subcategory = subcollection)

  ### create gene sets list
  all_gs_ids <- unique(msig_df$gs_id)
  msig_gsets_list <- list()
  for (id in all_gs_ids) {
    sub_df <- msig_df[msig_df$gs_id == id, ]
    msig_gsets_list[[id]] <- sub_df$gene_symbol
  }
  ### create gene sets descriptions
  msig_gsets_descriptions <- msig_df[, c("gs_name", "gs_id")]
  msig_gsets_descriptions <- unique(msig_gsets_descriptions)
  tmp <- msig_gsets_descriptions$gs_id
  msig_gsets_descriptions <- msig_gsets_descriptions$gs_name
  names(msig_gsets_descriptions) <- tmp

  result <- list(gene_sets = msig_gsets_list, descriptions = msig_gsets_descriptions)
  return(result)
}

#' Retrieve Organism-specific Gene Sets List
#'
#' @param source As of this version, either "KEGG", "Reactome" or "MSigDB" (default = "KEGG")
#' @param org_code (Used for "KEGG" only) KEGG organism code for the selected organism. For a full list
#' of all available organisms, see \url{https://www.genome.jp/kegg/catalog/org_list.html}
#' @param species (Used for "MSigDB" only) species name, such as Homo sapiens, Mus musculus, etc.
#' See \code{\link[msigdbr]{msigdbr_show_species}} for all the species available in
#' the msigdbr package (default = "Homo sapiens")
#' @param collection (Used for "MSigDB" only) collection. i.e., H, C1, C2, C3, C4, C5, C6, C7.
#' @param subcollection (Used for "MSigDB" only)  sub-collection, such as CGP, MIR, BP, etc. (default = NULL,
#' i.e. list all gene sets in collection)
#'
#' @return A list containing 2 elements: \itemize{
#' \item{gene_sets}{A list containing the genes involved in each gene set}
#' \item{descriptions}{A named vector containing the descriptions for each gene set}
#' }. For "KEGG" and "MSigDB", it is possible to choose a specific organism. For a full list
#' of all available KEGG organisms, see \url{https://www.genome.jp/kegg/catalog/org_list.html}.
#' See \code{\link[msigdbr]{msigdbr_show_species}} for all the species available in
#' the msigdbr package used for obtaining "MSigDB" gene sets.
#' For Reactome, there is only one collection of pathway gene sets.
#' @export
#'
get_gene_sets_list <- function(source = "KEGG", org_code = "hsa",
                               species = "Homo sapiens",
                               collection, subcollection = NULL) {
  if (source == "KEGG") {
    return(get_kegg_gsets(org_code))
  } else if (source == "Reactome") {
    message("For Reactome, there is only one collection of pathway gene sets.")
    return(get_reactome_gsets())
  } else if (source == "MSigDB") {
    return(get_mgsigdb_gsets(species = species,
                             collection = collection,
                             subcollection = subcollection))
  } else {
    stop("As of this version, this function is implemented to get data from KEGG, Reactome and MSigDB only")
  }
}
