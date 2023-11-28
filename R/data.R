#' Mutational information about Y-STR markers
#' 
#' A dataset from yhrd.org (and their sources) containing mutational information about
#' Y chromosomal short tandem repeat (Y-STR) markers used in forensic genetics.
#' 
#' Note, that loci with duplications (DYS385a/b as well as 
#' DYF387S1a/b have been split into two loci).
#' 
#' @format A data frame with 29 rows and 5 variables:
#' \describe{
#'   \item{Marker}{name of Y-STR marker}
#'   \item{Meioses}{number of meioses observed}
#'   \item{Mutations}{number of mutations observed in the corresponding number of Meioses}
#'   \item{MutProb}{point estimate of mutation probability, MutProb = Mutations/Meioses}
#'   \item{Alleles}{observed alleles}
#' }
#' @source \url{https://www.yhrd.org}
"ystr_markers"

#' Kit information about Y-STR markers
#' 
#' A dataset containing information about the 
#' Y chromosomal short tandem repeat (Y-STR) markers that are present in the kit.
#' 
#' @format A data frame with 88 rows and 2 variables:
#' \describe{
#'   \item{Marker}{name of Y-STR marker}
#'   \item{Kit}{name of Y-STR kit}
#' }
#' @source \url{https://www.yhrd.org}
"ystr_kits"
