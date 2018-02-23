## Test environments
* local OS X install, R 3.4.3
* Ubuntu 14.04.5 LTS (on travis-ci), R 3.4.2
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs. 

There were 2 NOTEs:
* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Ege Ulgen <egeulgen@gmail.com>'

  New submission

  Non-FOSS package license (file LICENSE)

  Possibly mis-spelled words in DESCRIPTION:
    Subnetworks (3:53)

* checking installed package size ... NOTE
  installed size is 10.0Mb
  sub-directories of 1Mb or more:
    extdata   9.4Mb

 This is our first submission. The data in this directory are protein-protein interaction networks, needed to run the active subnetwork search and other functions so we could not eliminate those over-sized files. 

## Downstream dependencies
There are currently no downstream dependencies for this package.