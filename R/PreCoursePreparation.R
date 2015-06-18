## # Course Preparation - Phylogenetics and comparative methods in R
## This document will guide you through the steps of preparing your R environment so that
## you can have a minimum number of headaches  working your way through the course. Please
## try to follow all steps before arriving for the course. New to R? A very useful site to 
## familiarize yourself with the "basics" (i.e. the advanced basics...) can be found [here](http://adv-r.had.co.nz/)

## ## Step 1. Install and update R to the latest version
## New versions of R are released periodically and many packages require a recent version. 
## Before you come to the course, make sure you download and install the latest version of
## R. 

## You can download R from CRAN here: http://www.r-project.org/

## Note: Unfortunately, updating R may require that all of the packages you previously installed
## to be reinstalled. See step 4 below for quickly installing most of the phylogenetic packages you
## may need in a single step.

## ## Step 2. Install Rstudio
## Rstudio is a development environment for R (one of many). We will expect that you use it
## during the course, but if you prefer to remain in your own environment that is fine, just 
## know that Rstudio has a number of convenient features and works across platforms, making 
## troubleshooting workshop a bit easier. 

## Rstudio can be downloaded here: http://www.rstudio.com/products/rstudio/download/

## ## Step 3. Setting up your machine to compile from source.
## Most of the time R packages can be installed from CRAN using pre-compiled binaries. However, 
## submitting to CRAN can be a bit of a hassle, and many packages hosted on CRAN may not be the most 
## current version or powerful. Increasingly, developers are making these packages available on 
## version control websites such as github or bitbucket. However, installing these packages often 
## requires that your machine has the tools necessary to compile them. Setting up these tools ahead
## of time will ensure that these packages can be easily installed from the web.

## Instructions for each operating system can be found here: https://support.rstudio.com/hc/en-us/articles/200486498-Package-Development-Prerequisites

## Test out that it works by using the following code in R:
##+eval=FALSE
install.packages("devtools")
library(devtools)
install_github("hadley/dplyr")

## ## Step 4. Installing packages. 
## There are a rapidly growing number of R packages that perform phylogenetic comparative methods
## in R. Want to install them all at once? Thanks to Brian O'Meara's [phylogenetics taskview](http://cran.r-project.org/web/views/Phylogenetics.html) you can.
## Note that this process takes a while, so make sure you have a good internet connection and plenty of
## time to let it install all the packages. 
##+eval=FALSE
install.packages("ctv")
install.views("Phylogenetics")

## ## Step 5 (Optional). Installing git
## Git is a version control software that is used by websites such as github and bitbucket. Learn about
## how to install git [here](https://help.github.com/articles/set-up-git/). 
## Rstudio has a number of built in tools to manage projects using git, which makes it especially nice to use. 
## Once you've installed git, get a [free github account](https://github.com/join). 
