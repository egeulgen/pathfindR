.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "##############################################################################
                        Welcome to pathfindR!

Please cite the article below if you use pathfindR in published reseach:

Ulgen E, Ozisik O, Sezerman OU. 2019. pathfindR: An R Package for Comprehensive
Identification of Enriched Pathways in Omics Data Through Active Subnetworks.
Front. Genet. doi:10.3389/fgene.2019.00858

##############################################################################"
  )

  ### Check Java version
  check_java_version()
}

#' Obtain Java Version
#'
#' @return character vector containing the output of "java -version"
#'
#' @details this function was adapted from the CRAN package \code{sparklyr}
fetch_java_version <- function() {
  java_home <- Sys.getenv("JAVA_HOME", unset = NA)
  if (!is.na(java_home)) {
    java <- file.path(java_home, "bin", "java")
    if (identical(.Platform$OS.type, "windows")) {
      java <- paste0(java, ".exe")
    }
    if (!file.exists(java)) {
      java <- ""
    }
  } else {
    java <- Sys.which("java")
  }

  if (java == "") {
    stop(
      "Java version not detected. Please download and install Java from ",
      dQuote("https://www.java.com/en/")
    )
  }

  version <- system2(java, "-version", stderr = TRUE, stdout = TRUE)
  if (length(version) < 1) {
    stop(
      "Java version not detected. Please download and install Java from ",
      dQuote("https://www.java.com/en/")
    )
  }

  return(version)
}

#' Check Java Version
#'
#' @param version character vector containing the output of "java -version". If
#' NULL, result of \code{\link{fetch_java_version}} is used (default = NULL)
#'
#' @return only parses and checks whether the java version is >= 1.8
#'
#' @details this function was adapted from the CRAN package \code{sparklyr}
check_java_version <- function(version = NULL) {
  if (is.null(version)) {
    version <- fetch_java_version()
  }

  # find line with version info
  versionLine <- version[grepl("version", version)]
  if (length(versionLine) != 1) {
    stop("Java version detected but couldn't parse version from ", paste(version, collapse = " - "))
  }

  # transform to usable R version string
  vers_string <- strsplit(versionLine, "\\s+", perl = TRUE)[[1]]
  vers_string <- vers_string[grepl("[0-9]+\\.[0-9]+\\.[0-9]+", vers_string, perl = TRUE)]
  if (length(vers_string) != 1) {
    vers_string <- strsplit(versionLine, "\\s+", perl = TRUE)[[1]]
    vers_string <- vers_string[grepl("[0-9]+", vers_string, perl = TRUE)]
    vers_string <- vers_string[!grepl("-", vers_string)]

    if (length(vers_string) != 1) {
      stop("Java version detected but couldn't parse version from: ", versionLine)
    }
  }

  parsedVersion <- gsub("^\"|\"$", "", vers_string)
  parsedVersion <- gsub("_", ".", parsedVersion)
  parsedVersion <- gsub("[^0-9.]+", "", parsedVersion)

  # ensure Java 1.8 (8) or higher
  if (utils::compareVersion(parsedVersion, "1.8") < 0) {
    stop(
      "Java version", parsedVersion, " detected but Java >=8 is required.
  Please download and install Java from ",
      dQuote("https://www.java.com/en/")
    )
  }
}
