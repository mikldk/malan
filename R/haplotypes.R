#' Generate a function to simulate pedigree founder haplotype based on ladder information
#' 
#' @param ladder_min vector of minimum alleles; `ladder_min[i]` is the minimum allele at locus i
#' @param ladder_max vector of minimum alleles; `ladder_max[i]` is the maximum allele at locus i
#' @export
generate_get_founder_haplotype_ladder <- function(ladder_min, ladder_max) {
  stopifnot(is.vector(ladder_min))
  stopifnot(is.vector(ladder_max))
  stopifnot(length(ladder_min) == length(ladder_max))
  
  function() {
    unlist(lapply(seq_along(ladder_min), function(loc) {
      sample(seq(ladder_min[loc], ladder_max[loc]), 1)
    }))
  }
}

#' Generate a function to simulate pedigree founder haplotype based on a haplotype databasep
#' 
#' @param db data frame or matrix with haplotypes from which the founder is randomly simulated
#' @export
generate_get_founder_haplotype_db <- function(db) {
  stopifnot(is.data.frame(db) | is.matrix(db))

  function() {
    db[sample(1L:nrow(db), 1), ]
  }
}

#' Generate a function to generate the zero haplotype
#' 
#' @param loci Number of loci
#' @export
get_zero_haplotype_generator <- function(loci) {
  f <- function() {
    rep(0L, loci)
  }
  return(f)
}
