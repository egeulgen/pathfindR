## Test environments
* local OS X 10.15.1, R 3.6.1
* Ubuntu 16.04.6 LTS (on Travis-CI), R 3.6.1
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:
* checking installed package size ... NOTE
    installed size is  6.8Mb
    sub-directories of 1Mb or more:
      R      3.1Mb
      data   1.2Mb
      doc    2.1Mb

  By turning the protein-protein interaction networks into internal R data (now residing under R), we reduced the size from approx. 14Mb to 7 Mb. If we later want to add more built-in data, we will be creating a separate dataset package (as suggested before).

## Downstream dependencies
There are currently no downstream dependencies for this package.
