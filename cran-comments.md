## Test environments
* local OS X 10.14.1, R 3.5.1
* Ubuntu 14.04.5 LTS (on travis-ci), R 3.5.1
* win-builder (devel and release)
* R-hub

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:
* checking installed package size ... NOTE
  installed size is 13.7Mb
  sub-directories of 1Mb or more:
    data      2.6Mb
    doc       1.3Mb
    extdata   9.5Mb

  The data in the directory "extdata" are protein-protein interaction networks, which are needed to run the active subnetwork search and other functions, so we could not eliminate those oversized files. The data in the "data" directory are example input and output, and individual gene sets. Since the gene sets are necessary for gene set enrichment analyses and are used by functions, we could not eliminate these files.

## Downstream dependencies
There are currently no downstream dependencies for this package.
