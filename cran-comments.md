## Test environments
* local macOS 15.7.4, R 4.6.0
* macOS-latest (GitHub Actions), R 4.6.0
* windows-latest (GitHub Actions), R 4.6.0
* ubuntu-latest (GitHub Actions), R 4.6.0
* ubuntu-latest (GitHub Actions), R oldrel
* ubuntu-latest (GitHub Actions), R devel
* win-builder (release and devel)
* r-hub (clang-ubsan)

## R CMD check results
0 errors | 0 warnings | 0 notes

This is a patch release. This release addresses the UBSAN issues reported on 
the CRAN additional checks (M1-SAN) at
https://www.stats.ox.ac.uk/pub/bdr/M1-SAN/pathfindR

Specifically, a signed integer overflow in the bundled C++ active subnetwork
search code (`src/java_compat.h`, in `java_string_hashcode`) has been fixed.
The hash accumulation now uses defined unsigned-integer arithmetic, which
removes the undefined behaviour while producing identical results. The same
class of issue was preemptively addressed in `JavaRandom::nextInt` and
`java_cap_for`.

## Downstream dependencies
There are currently no downstream dependencies for this package.
