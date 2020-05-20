## Test environments
* local OS X 10.15.4, R 4.0.0
* Ubuntu 16.04.6 LTS (on Travis-CI), R 4.0.0
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:
* checking installed package size ... NOTE
  installed size is  8.3Mb
  sub-directories of 1Mb or more:
    data   1.3Mb
    doc    3.8Mb
    R      2.7Mb

  As suggested before, if we later need to add more built-in data, we will be creating a separate dataset package .

## Downstream dependencies
  There are currently no downstream dependencies for this package.
