#' Infer meiotic distribution
#' 
#' @param poi Person of interest
#' @param indvs individuals in suspect population
#' 
#' @export
meiotic_distribution <- function(poi, indvs) {
  
  transfer_events <- sapply(indvs, meiotic_dist, ind2 = poi)
  
  # ...And not suspect (nor anybody not in same pedigree):
  in_ped <- which(transfer_events >= 1L)
  
  transfer_events_suspect_pop <- transfer_events[in_ped]
  
  M_set <- as.data.frame(table(transfer_events_suspect_pop))
  colnames(M_set) <- c("m", "n")
  M_set$m <- as.integer(as.character(M_set$m))
  N <- length(in_ped)
  M_set$p <- M_set$n / N
  M_set$cum_p <- cumsum(M_set$p)
  
  stopifnot(isTRUE(all.equal(sum(M_set$p), 1)))
  
  class(M_set) <- c("malan_meiotic_dist", class(M_set))
  
  return(M_set)
}

#' Calculate match probability information
#' 
#' Following [1], 
#'   $P(\text{match}) \approx \sum_{m = 1}^\infty 
#'   (1-\mu)^m \frac{P(M = m)}{P(M > 0)}$
#' for mutation probability $\mu$.
#'  
#' [1]: A. Caliebe and M. Krawczak (2018): 
#' "Match probabilities for Y-chromosomal profiles: A paradigm shift"
#' 
#' @param meiotic_dist Meiotic distribution (from [meiotic_distribution()])
#' @param mut_rate Overall profile mutation rate
#' 
#' @export
matchprob <- function(meiotic_dist, mut_rate) {
  stopifnot(is(meiotic_dist, "malan_meiotic_dist"))
  
  stopifnot(length(mut_rate) == 1L)
  stopifnot(mut_rate >= 0)
  stopifnot(mut_rate <= 1)
  
  # P(match) = \sum{m = 1}^\infty (1-\mu)^m \frac{P(M = m)}{P(M > 0)}

  matchprob <- sum(((1 - mut_rate)^meiotic_dist$m) * meiotic_dist$p)
  return(matchprob)
}


# FIXME: Speed up multiple? Distance matrix instead?

