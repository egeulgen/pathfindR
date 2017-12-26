if (!"org.Hs.eg.db" %in% installed.packages()) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("org.Hs.eg.db")
}

if (!"KEGGREST" %in% installed.packages()) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("KEGGREST")
}

if (!"Category" %in% installed.packages()) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("Category")
}

if (!"pathview" %in% installed.packages()) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("pathview")
}



if (!"foreach" %in% installed.packages())
  install.packages("foreach")

if (!"doSNOW" %in% installed.packages())
  install.packages("doSNOW")

if (!"rmarkdown" %in% installed.packages())
  install.packages("rmarkdown")

if (!"knitr" %in% installed.packages())
  install.packages("knitr")
