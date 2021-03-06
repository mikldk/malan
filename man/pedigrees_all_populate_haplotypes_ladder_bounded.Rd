% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{pedigrees_all_populate_haplotypes_ladder_bounded}
\alias{pedigrees_all_populate_haplotypes_ladder_bounded}
\title{Populate haplotypes in pedigrees (custom founder/bounded).}
\usage{
pedigrees_all_populate_haplotypes_ladder_bounded(
  pedigrees,
  mutation_rates,
  ladder_min,
  ladder_max,
  get_founder_haplotype = NULL,
  prob_two_step = 0,
  prob_genealogical_error = 0,
  progress = TRUE
)
}
\arguments{
\item{pedigrees}{Pedigree list in which to populate haplotypes}

\item{mutation_rates}{Vector with mutation rates}

\item{ladder_min}{Lower bounds for haplotypes, same length as \code{mutation_rates}}

\item{ladder_max}{Upper bounds for haplotypes, same length as \code{mutation_rates}; all entries must be strictly greater than \code{ladder_min}}

\item{get_founder_haplotype}{Function taking no arguments returning a haplotype of \code{length(mutation_rates)}}

\item{prob_two_step}{Given a mutation happens, this is the probability that the mutation is a two-step mutation; refer to details for information about behaviour around ladder boundaries}

\item{prob_genealogical_error}{Probability that a genealogical error happens: if so, give individual haplotype \code{get_founder_haplotype()} instead of father's}

\item{progress}{Show progress}
}
\description{
Populate haplotypes from founder and down in all pedigrees.
Note, that haplotypes are bounded by \code{ladder_min} and \code{ladder_max}.
All founders get a haplotype from calling the user
provided function \code{get_founder_haplotype()}.
}
\details{
Given that a two step mutation should happen (probability specified by \code{prob_two_step}):
With distances >= 2 to ladder bounds, mutations happen as usual.
At distance = 0 or 1 to a ladder bound, the mutation is forced to move away from the boundary.

Note, that pedigrees must first have been inferred by \code{\link[=build_pedigrees]{build_pedigrees()}}.
}
\examples{
sim <- sample_geneology(100, 10)
peds <- build_pedigrees(sim$population)
pedigrees_all_populate_haplotypes_ladder_bounded(
  peds, c(1, 1), c(0L, 0L), c(10L, 10L), 
  function(x) c(10, 10))
get_haplotype(sim$end_generation_individuals[[1]])

}
\seealso{
\code{\link[=pedigrees_all_populate_haplotypes]{pedigrees_all_populate_haplotypes()}} and
\code{\link[=pedigrees_all_populate_haplotypes_custom_founders]{pedigrees_all_populate_haplotypes_custom_founders()}}.
}
