## Test environments
* local OS X 13.3, R 4.3.0
* macOS-latest (on GitHub-Actions), R 4.3.0
* windows-latest (on GitHub-Actions), R 4.3.0
* ubuntu-latest (on GitHub-Actions), 4.2.3, R 4.3.0, devel
* win-builder (devel and release)
* R-hub (via check_for_cran())

## R CMD check results
  There were no ERRORs, WARNINGs or NOTEs.
  
  This is a patch version update for 'pathfindR' containing a bug fix to address
  problems on Debian systems, where a function attempted to write to user’s home 
  filespace. The relevant function was not essential and is now removed. Hence, 
  the issue should have been resolved.
  
## Downstream dependencies
  There are currently no downstream dependencies for this package.
