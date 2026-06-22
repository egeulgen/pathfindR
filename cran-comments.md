## Test environments
* local macOS 15.7.4, R 4.6.0
* macOS-latest (GitHub Actions), R 4.6.0
* windows-latest (GitHub Actions), R 4.6.0
* ubuntu-latest (GitHub Actions), R 4.6.0
* ubuntu-latest (GitHub Actions), R oldrel
* ubuntu-latest (GitHub Actions), R devel
* win-builder (release and devel)

## R CMD check results
0 errors | 0 warnings | 0 notes

This is a major release. It re-implements the active subnetwork search in R and
C++ (via Rcpp) and removes the Java dependency, so the package no longer
declares 'SystemRequirements: Java'. This release also renames a few exported
functions and includes minor fixes; see NEWS.md for details.

## Downstream dependencies
There are currently no downstream dependencies for this package.
