% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{pedigrees_all_populate_haplotypes}
\alias{pedigrees_all_populate_haplotypes}
\title{Populate haplotypes in pedigrees (0-founder/unbounded).}
\usage{
pedigrees_all_populate_haplotypes(
  pedigrees,
  loci,
  mutation_rates,
  prob_two_step = 0,
  prob_genealogical_error = 0,
  progress = TRUE
)
}
\arguments{
\item{pedigrees}{Pedigree list in which to populate haplotypes}

\item{loci}{Number of loci}

\item{mutation_rates}{Vector with mutation rates, length \code{loci}}

\item{prob_two_step}{Given a mutation happens, this is the probability that the mutation is a two-step mutation}

\item{prob_genealogical_error}{Probability that a genealogical error happens: if so, give individual haplotype \code{rep(0L, loci)} instead of father's}

\item{progress}{Show progress}
}
\description{
Populate haplotypes from founder and down in all pedigrees.
Note, that haplotypes are unbounded and
that all founders get haplotype \code{rep(0L, loci)}.
}
\details{
Note, that pedigrees must first have been inferred by \code{\link[=build_pedigrees]{build_pedigrees()}}.
}
\examples{
sim <- sample_geneology(100, 10)
peds <- build_pedigrees(sim$population)
pedigrees_all_populate_haplotypes(peds, 2, c(1, 1))
get_haplotype(sim$end_generation_individuals[[1]])

}
\seealso{
\code{\link[=pedigrees_all_populate_haplotypes_custom_founders]{pedigrees_all_populate_haplotypes_custom_founders()}} and
\code{\link[=pedigrees_all_populate_haplotypes_ladder_bounded]{pedigrees_all_populate_haplotypes_ladder_bounded()}}.
}
