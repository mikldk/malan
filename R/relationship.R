#' Construct M matrix
#' 
#' @param meioses number of meioses separating the two individuals
#' @param mu_dw mutation rate for 1-step down-mutation
#' @param mu_up mutation rate for 1-step up-mutation
#' 
#' @export
construct_M <- function(meioses, mu_dw, mu_up) {
  stopifnot(length(meioses) == 1L)
  stopifnot(isTRUE(all.equal(round(meioses), meioses)))
  meioses <- as.integer(meioses)
  stopifnot(meioses >= 0L)
  
  stopifnot(length(mu_dw) == 1L)
  stopifnot(mu_dw >= 0 && mu_dw <= 1)
  
  stopifnot(length(mu_up) == 1L)
  stopifnot(mu_up >= 0 && mu_up <= 1)
  
  stopifnot(mu_dw + mu_up <= 1)
  
  
  M <- diag(x = 1-(mu_up + mu_dw), nrow = 2*meioses + 1)
  
  for (i in seq_len(2*meioses)) {
    M[i, i+1] <- mu_dw
    M[i+1, i] <- mu_up
  }
  
  return(M)
}

#' Calculate distribution of allele difference
#' 
#' Calculate distribution of allele difference after `m` meioses.
#' 
#' @inheritParams construct_M
#' 
#' @return `data.frame` with columns `d` (allele difference) and `p` (prob)
#' 
#' @export
relationship_allele_diff_dist <- function(meioses, mu_dw, mu_up) {
  # parameters are validated in construct_M:
  M <- construct_M(meioses, mu_dw, mu_up)
  
  # NOTE: Exploit tridiaginal Toepliz matrix insted
  M_eig <- eigen(M, symmetric = TRUE)
  
  # M^k = P D^k P^-
  D <- diag(M_eig$values^meioses)
  P <- M_eig$vectors
  #Pinv <- solve(P)
  # P is orthogonal because all eigenvalues are distinct and then symmetric matrices have orthogonal eigenvectors
  Pinv <- t(P) 
  
  M_m <- P %*% D %*% Pinv
  p_0 <- matrix(0, 2*meioses + 1, 1)
  p_0[meioses+1] <- 1
  
  p_m <- M_m %*% p_0
  
  res <- data.frame(d = (-meioses):meioses, p = p_m)
  return(res)
}
