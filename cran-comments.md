## Test environments
* local OS X 10.13.3, R 3.5.0
* Ubuntu 14.04.5 LTS (on travis-ci), R 3.5.0

* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:
* checking installed package size ... NOTE
  installed size is 11.5Mb
  sub-directories of 1Mb or more:
    data      1.4Mb
    extdata   9.5Mb

  The data in the directory "extdata" are protein-protein interaction networks, which are needed to run the active subnetwork search and other functions, so we could not eliminate those oversized files.

## Downstream dependencies
There are currently no downstream dependencies for this package.