# pathfindr : An R Package for Pathway Enrichment Analysis Utilizing Active Subnetworks
pathfindr is a tool for pathway enrichment analysis via active subnetworks. The package also offers the option to cluster the resulting the enriched pathway and choose representative pathways.

![pathfindr Enrichment Workflow](./pathfindr.jpg?raw=true "pathfindr Enrichment Workflow")
This workflow takes in a data frame consisting of Gene Symbol, log-fold-change and adjusted-p values. After input testing, any gene symbols that are not in the PIN are converted to alias symbols if the alias is in the PIN. Next, active subnetwork search is performed. Pathway enrichment analysis is performed using the genes in each of the active subnetworks. Pathways with adjusted-p values lower than enrichment_threshold are discarded. The lowest adjusted-p value (over all subnetworks) for each pathway is kept. This process of active subnetwork search and enrichment is repeated for a selected number of iterations, which is done in parallel. Over all iterations, the lowest and the highest adjusted-p values, as well as number of occurences are reported for each enriched pathway.

![Pathway Clustering Workflow](./pw_clustering.jpg?raw=true "Pathway Clustering Workflow")
This function first calculates the pairwise distances between the pathways in the result_df data frame. Via a shiny HTML document, the hierarchical clustering dendrogram is visualized. In this HTML document, the user can select the value at which to cut the tree and the resulting representative pathways (chosen by smallest lowest p value) are presented as a table and pathways with cluster assignments are saved as a csv file to the current directory.

## Overview of the Workflow

## Dependencies
For the active subnetwork search component to work, the user must have [JAVA](https://www.java.com/en/download/manual.jsp) installed and path/to/java must be in the path.

### Resources
The protein interaction networks (PINs) were gathered from various resources:
- [Biogrid](https://downloads.thebiogrid.org/BioGRID)
- [BioPlex 2.0](http://bioplex.hms.harvard.edu/downloadInteractions.php)
- [GeneMania](http://genemania.org/data/)
- [STRING](https://string-db.org/cgi/download.pl?UserId=eCoJ8Fv0OZg6&sessionId=H23WzOsHidKz&species_text=Homo+sapiens)