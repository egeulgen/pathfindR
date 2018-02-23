## Test environments
* local OS X install, R 3.4.3
* Ubuntu 14.04.5 (on travis-ci), R 3.4.2
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* checking installed package size ... NOTE
  installed size is 10.0Mb
  sub-directories of 1Mb or more:
    extdata   9.4Mb

 This is our first submission. The data in this directory are protein-protein interaction networks, needed to run the active subnetwork search and other functions so we could not eliminate those over-sized files. 

## Downstream dependencies
There are currently no downstream dependencies for this package.