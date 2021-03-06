% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{get_pids_in_pedigree}
\alias{get_pids_in_pedigree}
\title{Get pids in pedigree}
\usage{
get_pids_in_pedigree(ped)
}
\arguments{
\item{ped}{Pedigree}
}
\description{
Get pids in pedigree
}
\examples{
sim <- sample_geneology(100, 10)
peds <- build_pedigrees(sim$population)
get_pids_in_pedigree(peds[[1]])

}
