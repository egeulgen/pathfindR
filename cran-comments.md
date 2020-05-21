## Test environments
* local OS X 10.15.4, R 4.0.0
* Ubuntu 16.04.6 LTS (on Travis-CI), R 4.0.0
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:
* checking installed package size ... NOTE
  installed size is  7.3Mb
  sub-directories of 1Mb or more:
    data   1.3Mb
    doc    2.8Mb
    R      2.7Mb
  
  As requested, we reduced the check time and the size of the package. However, 
  we could not eliminate or furtherr reduce the size of some data. The remaining 
  large data are protein-protein interaction networks, needed to run the active 
  subnetwork search (the core functionality) so we could not eliminate those 
  over-sized files. As suggested before, if we later need to add more built-in 
  data, we will be creating a separate dataset package.

## Downstream dependencies
  There are currently no downstream dependencies for this package.
