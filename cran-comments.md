## Test environments
* local OS X 10.15.1, R 3.6.1
* Ubuntu 16.04.6 LTS (on Travis-CI), R 3.6.1
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs. 

* There was only 1 NOTE in local and travis-ci checks:
* checking installed package size ... NOTE
  installed size is  6.8Mb
  sub-directories of 1Mb or more:
    R      3.1Mb
    data   1.2Mb
    doc    2.1Mb

  If we later need to add more built-in data, we will be creating a separate dataset package (as suggested before).

* Upon the warning of Kurt Hornik, with the current r-devel, pathfindR threw an error because it was using `class(.) == *` to check the class of an object, this was fixed by checking the classes more properly.

## Downstream dependencies
  There are currently no downstream dependencies for this package.
