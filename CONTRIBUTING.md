# Contributing to pathfindR development

The goal of this guide is to help you in contributing to pathfindR. The guide is divided into two main pieces:

1. Filing a bug report or feature request in an issue.
1. Suggesting a change via a pull request.

Please note that pathfindR is released with a [Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this project, 
you agree to abide by its terms.

## Issues

When filing an issue, the most important thing is to include a minimal 
reproducible example so that we can quickly verify the problem, and then figure 
out how to fix it. There are three things you need to include to make your 
example reproducible: required packages, data, code.

1.  **Packages** should be loaded at the top of the script, so it's easy to
    see which ones the example needs.
  
1.  The easiest way to include **data** is to use `dput()` to generate the R code 
    to recreate it. For example, to recreate the `mtcars` dataset in R,
    I'd perform the following steps:
  
       1. Run `dput(mtcars)` in R
       2. Copy the output
       3. In my reproducible script, type `mtcars <- ` then paste.
       
    But even better is if you can create a `data.frame()` with just a handful
    of rows and columns that still illustrates the problem.
    
    For more complex **data**, you can use `saveRDS()` to save the object and attach it with the issue.
  
1.  Spend a little bit of time ensuring that your **code** is easy for others to
    read:
  
    * make sure you've used spaces and your variable names are concise, but
      informative
  
    * use comments to indicate where your problem lies
  
    * do your best to remove everything that is not related to the problem.  
     The shorter your code is, the easier it is to understand.

You can check you have actually made a reproducible example by starting up a 
fresh R session and pasting your script in.

## Pull requests

To contribute a change to pathfindR, you follow these steps:

1. Create a branch in git and make your changes.
1. Push branch to github and issue pull request (PR).
1. Discuss the pull request.
1. Iterate until either we accept the PR or decide that it's not
   a good fit for pathfindR.

If you're not familiar with git or github, please start by reading <http://r-pkgs.had.co.nz/git.html>

# Attribution
This Contributing guide was adapted from [ggplot2](https://github.com/tidyverse/ggplot2)
