parsejActive <- function(jactive_output, signif_genes) {
  ends <- grep("SPOTPvaluesig", jactive_output$V1)
  starts <- c(1, ends + 1)
  starts <- starts[-length(starts)]
  ends <- ends - 1

  subnetworks <- list()
  scores <- c()
  for (i in 1:length(starts)){
    subnetworks[[i]] <- jactive_output$V1[starts[i]:ends[i]][-1]
    scores <- c(scores, as.numeric(jactive_output$V1[starts[i]:ends[i]][1]))
  }

  # exclude subnetworks with negative or zero score
  subnetworks <- subnetworks[scores > 0]
  # select subnetworks with at least 2 significant genes
  subnetworks <- subnetworks[sapply(subnetworks,
                                    function(x) sum(x %in% signif_genes)) >= 2]

  return(subnetworks)
}
