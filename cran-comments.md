## Test environments
* local OS X 10.13.3, R 3.4.3
* Ubuntu 14.04.5 LTS (on travis-ci), R 3.4.2
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs. 

There were 2 NOTEs:
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Ege Ulgen <egeulgen@gmail.com>'

New submission

Possibly mis-spelled words in DESCRIPTION:
  OU (14:33)
  Ozisik (14:14)
  Sezerman (14:24)
  Subnetworks (3:53, 15:42)
  Ulgen (14:5)
  bioRxiv (15:55)
  pathfindR (8:31, 12:5, 14:43)
  subnetworks (9:22, 9:56, 11:63)

* checking installed package size ... NOTE
  installed size is 10.2Mb
  sub-directories of 1Mb or more:
    extdata   9.5Mb

 This is our first submission. The data in the directory "extdata" are protein-protein interaction networks, which are needed to run the active subnetwork search and other functions, so we could not eliminate those oversized files.

## Downstream dependencies
There are currently no downstream dependencies for this package.