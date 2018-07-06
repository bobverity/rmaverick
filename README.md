
<!-- README.md is generated from README.Rmd. Please edit that file -->
rmaverick
=========

<!--
[![Build Status](https://travis-ci.org/bobverity/rmaverick.png?branch=develop)](https://travis-ci.org/bobverity/rmaverick)
[![Coverage status](https://codecov.io/gh/bobverity/rmaverick/branch/develop/graph/badge.svg)](https://codecov.io/github/bobverity/rmaverick?branch=develop)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/bobverity/rmaverick?branch=develop&svg=true)](https://ci.appveyor.com/project/bobverity/rmaverick)
-->
The goal of *rmaverick* is to infer population structure from genetic data. What makes *rmaverick* different from other similar programs is its ability to estimate the *evidence* for different numbers of sub-populations (K), and even different evolutionary models, through a method called generalised thermodynamic integration (GTI). *rmaverick* is an updated version of the earlier [MavericK](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4981280/) program, and comes with a number of new features that make it more powerful and easier to use.

Installation
------------

*rmaverick* relies on the [Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html) package, which requires the following OS-specific steps:

-   Windows
    -   Download and install the appropriate version of [Rtools](https://cran.rstudio.com/bin/windows/Rtools/) for your version of R. On installation, ensure you check the box to arrange your system PATH as recommended by Rtools
-   Mac OS X
    -   Download and install [XCode](http://itunes.apple.com/us/app/xcode/id497799835?mt=12)
    -   Within XCode go to Preferences : Downloads and install the Command Line Tools
-   Linux (Debian/Ubuntu)
    -   Install the core software development utilities required for R package development as well as LaTeX by executing

            sudo apt-get install r-base-dev texlive-full

Next, in R, ensure that you have the [devtools](https://www.rstudio.com/products/rpackages/devtools/) package installed by running

``` r
install.packages("devtools", repos='http://cran.us.r-project.org')
```

Then install the *rmaverick* package directly from GitHub by running

``` r
devtools::install_github("bobverity/rmaverick")
```

Finally, we need to load the package:

``` r
library(rmaverick)
```

Example analysis
----------------

This example demonstrates the complete analysis pipeline, including importing data, running the main MCMC, diagnosing good and bad MCMC behaviour, comparing different models, and producing basic outputs.

### Simulate some data

rmverick comes with built-in functions for simulating data from different evolutionary models. The models used in simulation are exactly the same as the models used in the inference step, allowing us to test the power of the program without worrying about discrepancies between the data and the assumed model. We will simulate a data set of 100 diploid individuals from the no-admixture model, each genotyped at 20 loci and originating from 3 distinct subpopulations:

``` r
mysim <- sim_data(n = 100, loci = 20, K = 3, admix_on = FALSE)
```

Running `names(mysim)` we can see that the simulated data contains several elements:

``` r
names(mysim)
#> [1] "data"         "allele_freqs" "admix_freqs"  "group"
```

As well as the raw data, we have a record of the allele frequencies, admixture frequencies and grouping that were used in generating these data. These can be useful in ground-truthing our estimated values later on.

Running `head(mysim$data)` we can see the general data format required by *rmaverick*:

``` r
head(mysim$data)
#>     ID pop ploidy locus1 locus2 locus3 locus4 locus5 locus6 locus7 locus8
#> 1 ind1   1      2      2      1      5      5      3      3      4      3
#> 2 ind1   1      2      5      2      5      2      2      1      5      3
#> 3 ind2   1      2      2      3      4      5      3      3      4      4
#> 4 ind2   1      2      4      2      1      2      4      1      2      3
#> 5 ind3   1      2      2      1      1      1      2      4      4      5
#> 6 ind3   1      2      3      3      5      2      1      4      2      4
#>   locus9 locus10 locus11 locus12 locus13 locus14 locus15 locus16 locus17
#> 1      5       1       3       5       2       1       5       3       2
#> 2      5       1       5       1       3       1       4       3       2
#> 3      5       1       3       1       3       3       4       3       2
#> 4      3       1       3       5       1       2       4       3       1
#> 5      3       1       5       4       3       5       4       2       4
#> 6      2       1       2       1       5       5       4       3       3
#>   locus18 locus19 locus20
#> 1       4       1       2
#> 2       4       1       2
#> 3       1       1       2
#> 4       1       1       2
#> 5       1       3       2
#> 6       4       1       2
```

Samples are in rows and loci are in columns, meaning for polyploid individuals multiple rows must be used for the same individual (here two rows per individual). There are also several meta-data columns, including the sample ID, the population of origin of the sample and the ploidy. These meta-data columns are optional and can be turned on or off when loading the data into a project.

### Create a project and read in data

*rmaverick* works with projects, which are essentially just simple lists containing all the inputs and outputs of a given analysis. We start by creating a project and loading in our data:

``` r
myproj <- mavproject()
myproj <- bind_data(myproj, mysim$data, ID_col = 1, pop_col = 2, ploidy_col = 3)
```

Notice the general format of the `bind_data()` function, which takes the same project as both input and output. This is the format that most *rmaverick* functions will take, as it allows a function to modify the project before overwriting the original version. Notice also that we have specified which columns are meta-data, and all other columns are assumed to contain genetic data. We can produce a summary of the project to check that the data have been loaded in correctly:

``` r
summary(myproj)
#> DATA:
#>    individuals = 100
#>    loci = 20
#>    ploidy = 2
#>    pops = 3 sampled populations
#>    missing data = 0 of 4000 gene copies (0%)
#> 
#> PARAMETER SETS:
#>    (none defined)
#> 
#> OUTPUT:
#>    (none saved)
```

If data has been input incorrectly, for example if mata-data columns have not been specified and so have been interpreted as genetic data, then this should be visible at this stage.

### Define parameters and run basic MCMC

We can define different evolutionary models by using different parameter sets. Our first parameter set will represent a simple no-admixture model (the same model that was used to generate the data).

``` r
myproj <- new_set(myproj, name = "correct model (no admixture)", admix_on = FALSE)
```

Producing a summary of the project we can now see additional properties, including which is the current active set and the key parameters of this set.

``` r
summary(myproj)
#> DATA:
#>    individuals = 100
#>    loci = 20
#>    ploidy = 2
#>    pops = 3 sampled populations
#>    missing data = 0 of 4000 gene copies (0%)
#> 
#> PARAMETER SETS:
#>  * SET1: correct model (no admixture)
#> 
#> ACTIVE SET: SET1
#>    model = no-admixture
#> 
#> OUTPUT:
#>    (none saved)
```

Now we are ready to run a basic MCMC. We will start by exploring values of K from 1 to 5, using 1000 burn-in iterations and 1000 sampling iterations. By default the MCMC has `auto_converge` turned on, meaning it will test for convergence every `burnin/10` iterations and will exit if convergence is reached. Hence, it is generally a good idea to set `burnin` to be higher than expected, as the MCMC will adjust this number down if needed. The number of sampling iterations can also be tuned. Our aim when choosing the number of sampling iterations should be to obtain enough samples that our posterior estimates are accurate to a given tolerance level, but not so many that we waste time running the MCMC for long periods past this point. We will look into this parameter again once the initial MCMC has completed. The most unfamiliar parameter for most users will be the number of "rungs". *rmaverick* runs multiple MCMC chains simultaneously, each at a different rung on a temperature ladder. The cold chain is our ordinary MCMC chain, and the hot chains serve two purposes: 1) they improve MCMC mixing by allowing it to explore the space more freely, 2) they are central to the GTI method of estimating the evidence for different models. Finally, for the sake of this document we will run with `pb_markdown=TRUE` to avoid printing large amounts of output, but ordinarily this argument should be omitted.

``` r
myproj <- run_mcmc(myproj, K = 1:5, burnin = 1e3, samples = 1e3, rungs = 10, pb_markdown =  TRUE)
#> Calculating exact solution for K = 1
#>    completed in 0.00208883 seconds
#> 
#> Running MCMC for K = 2
#> Burn-in phase
#> 
  |                                                                       
  |=================================================================| 100%
#>    converged within 100 iterations
#> Sampling phase
#> 
  |                                                                       
  |=================================================================| 100%
#>    completed in 0.810188 seconds
#> 
#> Running MCMC for K = 3
#> Burn-in phase
#> 
  |                                                                       
  |=================================================================| 100%
#>    converged within 100 iterations
#> Sampling phase
#> 
  |                                                                       
  |=================================================================| 100%
#>    completed in 1.07598 seconds
#> 
#> Running MCMC for K = 4
#> Burn-in phase
#> 
  |                                                                       
  |=================================================================| 100%
#>    converged within 100 iterations
#> Sampling phase
#> 
  |                                                                       
  |=================================================================| 100%
#>    completed in 1.31234 seconds
#> 
#> Running MCMC for K = 5
#> Burn-in phase
#> 
  |                                                                       
  |=================================================================| 100%
#>    converged within 100 iterations
#> Sampling phase
#> 
  |                                                                       
  |=================================================================| 100%
#>    completed in 1.54187 seconds
#> 
#> Processing results
#> Total run-time: 5.84 seconds
```

Notice that the solution for K=1 is almost instantaneous, as no MCMC is required in this case. In the example above we also saw that the MCMC converged within far fewer than the full 1000 burn-in iterations.

Next, we need to check that we are happy with the behaviour of our MCMC. This is an important step, as without it our results may not be meaningful despite having thounsands of posterior samples. We will check the MCMC in three ways:

1.  By looking at levels of auto-correlation
2.  By measuring the standard error of posterior estimates
3.  By looking at the GTI "path"

#### 1. Measuring auto-correlation

TODO!

#### 2. Posterior accuracy

TODO!

#### 3. GTI path

The GTI method estimates the evidence for a given model by combining information across multiple temperature rungs. These rungs provide a series of point estimates that together make a "path", and the final evidence estimate is computed from the area between this path and the zero-line. We can visualise this path using the `plot_GTI_path()` function:

``` r
plot_GTI_path(myproj, K = 3)
```

![](README-unnamed-chunk-14-1.png)

In order for our evidence estimate to be unbiased it is important that this path is relatively smooth. We can modify the smoothness of the path in two ways: 1) by increasing the number of `rungs` used in the MCMC, 2) by changing the value of `GTI_pow` which controls the curvature of the path - higher values lead to more steep curvature. Ideally we want a straight path, i.e. we want as little curvature as possible. In the example above we have a good number of rungs, and we might consider using a lower `GTI_pow` to reduce the curvature of the path slightly, but overall this path is good enough so there is no need to re-run the MCMC. This check should be performed on every value of K.

Once we are happy with our GTI paths we can look at our evidence estimates, first of all in log space:

``` r
plot_logevidence_K(myproj)
```

![](README-unnamed-chunk-15-1.png)

We can see a clear signal for K=3 or higher, and the 95% credible intervals are nice and tight. If we needed tighter credible intervals at this stage we could re-run the MCMC with a larger number of `samples`.

We can also plot the full posterior distribution of K, which is obtained by transforming these values out of log space and normalising to sum to one:

``` r
plot_posterior_K(myproj)
```

![](README-unnamed-chunk-16-1.png)

This second plot is usually more straightfoward to interpret, as it is in linear space and so can be understood in terms of ordinary probability. In this example we can see strong evidence for K=3, with a posterior probability of &gt;0.99. Again, if we had seen wide credible intervals at this stage then it would have be worth repeating the MCMC with a larger number of `samples`, but in this case the result is clear and so there is no need.

### Plot results

The main result of interest is usually the posterior allocation or "structure" plot. This plot contains one bar for each individual, with the proportion of each colour giving the posterior probability of belonging to each of the K subpopulations. We can use the `plot_qmatrix_ind()` function to produce posterior allocation plots for different values of K:

``` r
for (i in 2:5) {
  plot_qmatrix_ind(myproj, K = i)
}
```

![](README-unnamed-chunk-17-1.png)![](README-unnamed-chunk-17-2.png)![](README-unnamed-chunk-17-3.png)![](README-unnamed-chunk-17-4.png)

We can see that for K=3 there is a clear split into three distinct subpopulations (we can verify that this is the correct grouping by looking at `mysim$group`). When reporting and publishing results it is a good idea to produce posterior allocation plots for a range of values of K so that the reader has the option visualise structure at multiple levels, but ideally this should be backed up by a plot of the model evidence to give some idea of the model fit at each level. At this stage it is worth stressing the point made by many previous authors - that the model used by *rmaverick* and similar programs is just a cartoon of reality, and that there is no strict K in the real world. Instead, each K captures a different level of population structure, and while the evidence can help guide us towards values of K that capture the data well, it is just a guide and should be taken alongside other biological considerations.

### Admixture model

We will also run the simulated data through the admixture model, partly to explore how this model differs in implementation from the simpler no-admixture model, and partly to test whether *rmaverick* can detect that this model is incorrect based on the model evidence.

We start by creating a new parameter set, this time under the admixture model and with the parameter alpha free to be estimated by the MCMC. Alpha controls the level of admixture, with larger values indicating a greater probability of inheriting genes from multiple populations (we arrive back at the no-admixture model when alpha is zero).

``` r
myproj <- new_set(myproj, name = "incorrect model (admixture)", admix_on = TRUE, estimate_alpha = TRUE)
summary(myproj)
#> DATA:
#>    individuals = 100
#>    loci = 20
#>    ploidy = 2
#>    pops = 3 sampled populations
#>    missing data = 0 of 4000 gene copies (0%)
#> 
#> PARAMETER SETS:
#>    SET1: correct model (no admixture)
#>  * SET2: incorrect model (admixture)
#> 
#> ACTIVE SET: SET2
#>    model = admixture
#>    estimate alpha = TRUE
#> 
#> OUTPUT:
#>    (none saved)
```

We run the MCMC the same way as before:

``` r
myproj <- run_mcmc(myproj, K = 1:5, burnin = 1e3, samples = 1e3, rungs = 10, pb_markdown = TRUE)
#> Calculating exact solution for K = 1
#>    completed in 0.0739005 seconds
#> 
#> Running MCMC for K = 2
#> Burn-in phase
#> 
  |                                                                       
  |=================================================================| 100%
#>    Warning: convergence still not reached within 1000 iterations
#> Sampling phase
#> 
  |                                                                       
  |=================================================================| 100%
#>    completed in 27.3279 seconds
#> 
#> Running MCMC for K = 3
#> Burn-in phase
#> 
  |                                                                       
  |=================================================================| 100%
#>    converged within 1000 iterations
#> Sampling phase
#> 
  |                                                                       
  |=================================================================| 100%
#>    completed in 30.503 seconds
#> 
#> Running MCMC for K = 4
#> Burn-in phase
#> 
  |                                                                       
  |=================================================================| 100%
#>    Warning: convergence still not reached within 1000 iterations
#> Sampling phase
#> 
  |                                                                       
  |=================================================================| 100%
#>    completed in 32.8582 seconds
#> 
#> Running MCMC for K = 5
#> Burn-in phase
#> 
  |                                                                       
  |=================================================================| 100%
#>    Warning: convergence still not reached within 1000 iterations
#> Sampling phase
#> 
  |                                                                       
  |=================================================================| 100%
#>    completed in 37.0172 seconds
#> 
#> Processing results
#> Total run-time: 2.15 minutes
```

Notice that the MCMC took considerably longer this time round. This is because the MCMC now has to consider the population assignment of each gene copy seperately, and also has additional parameters to estimate. For this reason we suggest that if there is no prior biological reason to expect admixture between populations then consider running the simpler non-admixture model only.

As before, we need to check the behaviour of our MCMC for all values of K.

TODO - more checks at this point!

``` r
plot_GTI_path(myproj, K = 3)
```

![](README-unnamed-chunk-20-1.png)

As before, we can visualise the evidence for each value of K in log space, and in terms of the full posterior distribution:

``` r
plot_logevidence_K(myproj)
```

![](README-unnamed-chunk-21-1.png)

``` r
plot_posterior_K(myproj)
```

![](README-unnamed-chunk-21-2.png)

Again, there is clear evidence for K=3 under this model, despite the model being incorrect for the data.

Producing posterior allocation plots for a range of values of K, we see results similar to the no-admixture model, with most individuals assigned to just a single subpopulation. This is a further clue that the no-admixture model is more appropriate, as even the admixture model is attempting to converge to this simpler state.

``` r
for (i in 2:5) {
  plot_qmatrix_ind(myproj, K = i)
}
```

![](README-unnamed-chunk-22-1.png)![](README-unnamed-chunk-22-2.png)![](README-unnamed-chunk-22-3.png)![](README-unnamed-chunk-22-4.png)

### Comparing evolutionary models

Finally, we can make use of one of the major advantages of the model evidence - the ability to compare between different evolutionary models. The overall evidence of a model is taken as the average over all values of K under that model. We can plot the model evidence in log space, and in the form of a posterior distribution after transforming to linear space and normalising to sum to one:

``` r
plot_logevidence_model(myproj)
```

![](README-unnamed-chunk-23-1.png)

``` r
plot_posterior_model(myproj)
```

![](README-unnamed-chunk-23-2.png)

It is clear that there is far greater support for the no-admixture model, with a posterior probability of &gt;0.99. Therefore we would be justified in ignoring the output from the admixture model, and only reporting results of the no-admixture model.

### Final notes

-   If we are not interested in estimating K or comparing models then there is technically no need to use multiple temperature rungs, and so the MCMC can be made to run faster by setting `rungs=1`.

-   Can be run on a cluster (notes to follow soon)
