# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
# R scripts generated  Sat Dec 30 15:44:22 EST 2017

################################################################
#   Differential expression analysis with limma
library(Biobase)
library(GEOquery)
library(limma)

# load series and platform data from GEO

gset <- getGEO("GSE61635", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- phenoData(gset)$title
sml <- ifelse(grepl("SLE", gsms), "SLE", "CONTROL")

# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# set up the data and proceed with analysis
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(SLE-CONTROL, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)

tT <- topTable(fit2, adjust="fdr", sort.by="P", number=Inf)
tT <- tT[abs(tT$logFC) >= 1, ]
tT <- tT[tT$adj.P.Val <= 0.05, ]
missing <- tT$ID[!grepl("^[a-zA-Z0-9_-]*$", tT$Gene.symbol)]
library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
syms <- getBM(attributes = c('affy_hg_u133a_2', 'hgnc_symbol'),
                  filters = 'affy_hg_u133a_2',
                  values = missing,
                  mart = ensembl)
syms <- syms[syms$hgnc_symbol != "",]
syms <- syms[!duplicated(syms$affy_hg_u133a_2),]

tT$Gene.symbol[match(syms$affy_hg_u133a_2, tT$ID)] <- syms$hgnc_symbol
tT$Gene.symbol[!grepl("^[a-zA-Z0-9_-]*$", tT$Gene.symbol)] <- sapply(tT$Gene.symbol[!grepl("^[a-zA-Z0-9_-]*$", tT$Gene.symbol)],
                                                                     function(x) unlist(strsplit(x, "\\/\\/\\/"))[1])


final <- subset(tT, select=c("Gene.symbol", "logFC", "adj.P.Val"))
final <- final[order(final$adj.P.Val), ]

final <- final[final$Gene.symbol != "", ]
final <- final[!duplicated(final$Gene.symbol), ]

write.csv(final, "SLE.csv", row.names=F)
example_input <- final
save(example_input, file = "../data/example_input.RData")
