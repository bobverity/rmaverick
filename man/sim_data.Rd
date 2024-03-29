% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\name{sim_data}
\alias{sim_data}
\title{Simulate data}
\usage{
sim_data(
  n = 10,
  ploidy = 2,
  loci = 10,
  allele_num = 5,
  lambda = 1,
  admix_on = FALSE,
  alpha = 1,
  K = 3,
  pop_col_on = TRUE,
  ploidy_col_on = TRUE
)
}
\arguments{
\item{n}{number of samples}

\item{ploidy}{ploidy of every sample. If scalar then the same value applies 
to all samples}

\item{loci}{number of loci}

\item{allele_num}{number of alleles per locus}

\item{lambda}{prior parameter on allele frequency distribution}

\item{admix_on}{whether to simulate from the admixture model}

\item{alpha}{the admixture parameter}

\item{K}{the number of subpopulations}

\item{pop_col_on}{whether to simulate data with a population column}

\item{ploidy_col_on}{whether to simulate data with a ploidy column}
}
\description{
Simulate data from the same model(s) used in the rmaverick
  inference step. Data can be generated from models with or without
  admixture.
}
