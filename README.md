# pathfindr : An R Package for Pathway Enrichment Analysis Utilizing Active Subnetworks
pathfindr is a tool for pathway enrichment analysis via active subnetworks. The package also offers the option to cluster the resulting the enriched pathway and choose representative pathways.

## Overview of the Enrichment Workflow

![pathfindr Enrichment Workflow](./vignettes/pathfindr.png?raw=true "pathfindr Enrichment Workflow")
This workflow takes in a data frame consisting of Gene Symbol, log-fold-change and adjusted-p values. After input testing, any gene symbols that are not in the PIN are converted to alias symbols if the alias is in the PIN. Next, active subnetwork search is performed. Pathway enrichment analyses are performed using the genes in each of the active subnetworks. Pathways with adjusted-p values larger than enrichment_threshold are discarded. The lowest adjusted-p value (over all subnetworks) for each pathway is kept. This process of active subnetwork search and enrichment analyses is repeated for a selected number of iterations, which is done in parallel. Over all iterations, the lowest and the highest adjusted-p values, as well as number of occurences are reported for each enriched pathway.

## Overview of the Clustering Workflow

![Pathway Clustering Workflow](./vignettes/pw_clustering.png?raw=true "Pathway Clustering Workflow")
This workflow first calculates the pairwise distances between the pathways in the resulting data frame. Via a shiny app, presented as an HTML document, the hierarchical clustering dendrogram is visualized. In this HTML document, the user can select the agglomeration method and the distance value at which to cut the tree. The dendrogram with the cut-off value marked with a red line is dynamically visualized and the resulting cluster assignments of the pathways along with annotation of representative pathways (chosen by smallest lowest p value) are presented as a table. This table can be saved as a csv file via pressing the button `Get Pathways w\ Cluster Info`.

## Dependencies
For the active subnetwork search component to work, the user must have [JAVA](https://www.java.com/en/download/manual.jsp) installed and path/to/java must be in the PATH environment variable.

### Resources
The protein interaction networks (PINs) were gathered from various resources:
- [Biogrid](https://downloads.thebiogrid.org/BioGRID)
- [BioPlex 2.0](http://bioplex.hms.harvard.edu/downloadInteractions.php)
- [GeneMania](http://genemania.org/data/): only interactions with weights >= 0.0006 were kept
- [STRING](https://string-db.org/cgi/download.pl?UserId=eCoJ8Fv0OZg6&sessionId=H23WzOsHidKz&species_text=Homo+sapiens): only interactions with combined_score >= 800 were kept.
- [IntAct](https://www.ebi.ac.uk/intact/)
- [HitPredict](http://hintdb.hgc.jp/htp/download.html)
- KEGG PIN - created via an in-house script
