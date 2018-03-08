## Test environments
* local OS X 10.13.3, R 3.4.3
* Ubuntu 14.04.5 LTS (on travis-ci), R 3.4.2
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs. 

There were 3 NOTEs:
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Ege Ulgen <egeulgen@gmail.com>'

New submission

License components with restrictions and base license permitting such:
  MIT + file LICENSE
File 'LICENSE':
  MIT License
  
  Copyright (c) 2018 Ege Ulgen, Ozan Ozisik
  
  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:
  
  The above copyright notice and this permission notice shall be included in all
  copies or substantial portions of the Software.
  
  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.

Possibly mis-spelled words in DESCRIPTION:
  OU (12:70)
  Ozisik (12:51)
  Sezerman (12:61)
  Subnetworks (3:53, 14:12)
  Ulgen (12:42)
  bioRxiv (14:25)
  pathfindR (7:14, 13:5)

Found the following (possibly) invalid URLs:
  URL: https://doi.org/10.1101/272450
    From: inst/doc/pathfindr_vignette.html
    Status: 404
    Message: Not Found

Found the following (possibly) invalid DOIs:
  DOI: 10.1101/272450
    From: DESCRIPTION
    Status: Not Found
    Message: 404

* checking installed package size ... NOTE
  installed size is 10.2Mb
  sub-directories of 1Mb or more:
    extdata   9.5Mb

* checking DESCRIPTION meta-information ... NOTE
Malformed Description field: should contain one or more complete sentences.

 This is our first submission. The DOI was provided by bioRxiv and the link will likely take a while to be up (see https://www.biorxiv.org/content/early/2018/03/07/272450). The data in the directory "extdata" are protein-protein interaction networks, needed to run the active subnetwork search and other functions so we could not eliminate those over-sized files. The final NOTE states "should contain one or more complete sentences" but there are 5 complete sentences in the Description field.

## Downstream dependencies
There are currently no downstream dependencies for this package.