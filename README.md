
<!-- README.md is generated from README.Rmd. Please edit that file -->
rmaverick
=========

The goal of rmaverick is to infer population structure from genetic data. What makes rmaverick different from other similar programs is its ability to estimate the *evidence* for different numbers of sub-populations (K), and even different evolutionary models, through a method called thermodynamic integration. It also contains some advanced MCMC methods that ensure reliable results.

Installation
------------

rmaverick relies on the Rcpp package, which requires the following OS-specific steps:

-   Windows
    -   Download and install the appropriate version of [Rtools](https://cran.rstudio.com/bin/windows/Rtools/) for your version of R. On installation, ensure you check the box to arrange your system PATH as recommended by Rtools
-   Mac OS X
    -   Download and install [XCode](http://itunes.apple.com/us/app/xcode/id497799835?mt=12)
    -   Within XCode go to Preferences : Downloads and install the Command Line Tools
-   Linux (Debian/Ubuntu)
    -   Install the core software development utilities required for R package development as well as LaTeX by executing

            ```
            sudo apt-get install r-base-dev texlive-full
            ```

Next, in R, ensure that you have the devtools package installed by running

``` r
install.packages("devtools", repos='http://cran.us.r-project.org')
```

Then we can simply install the rmaverick package directly from GitHub by running

``` r
devtools::install_github("bobverity/rmaverick")
```

Finally, we need to load the package

``` r
library(rmaverick)
```

Example
-------

This is a basic example which shows you how to solve a common problem:

``` r
## basic example code
plot(1:10)
```

![](README-example-1.png)
