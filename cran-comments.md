## Test environments
* local OS X 10.15.1, R 3.6.1
* Ubuntu 16.04.6 LTS (on Travis-CI), R 3.6.1
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs. 

* There was only 1 NOTE:
  checking installed package size ... NOTE
  installed size is  6.8Mb
  sub-directories of 1Mb or more:
    R      3.1Mb
    data   1.2Mb
    doc    2.1Mb

  If we later need to add more built-in data, we will be creating a separate dataset package (as suggested before).

* Upon the warning of Prof. Brian D. Ripley, the Java version in SystemRequirements was corrected to "Java (>= 8.0)" and the Java version is checked via rJava (as suggested on https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Writing-portable-packages). 

## Downstream dependencies
  There are currently no downstream dependencies for this package.
