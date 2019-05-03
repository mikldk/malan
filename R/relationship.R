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
  
  
  M <- diag(x = 1 - (mu_up + mu_dw), nrow = 2*meioses + 1)
  
  for (i in seq_len(2*meioses)) {
    M[i, i+1] <- mu_dw
    M[i+1, i] <- mu_up
  }
  
  return(M)
}

relationship_allele_diff_dist_sym_worker <- function(meioses, P, D, Pinv) {
  M_m <- P %*% D %*% Pinv
  p_0 <- matrix(0, 2*meioses + 1, 1)
  p_0[meioses + 1] <- 1
  
  p_m <- M_m %*% p_0

  # FIXME: Numerics?  
  if (max(abs(Im(p_m))) > 1e-8) {
    print(p_m[which(abs(Im(p_m)) > 1e-8)])
    stop("Imaginary probability occured. ", 
         "Maybe the difference between mu_dw and mu_up causes numerical problems.")
  }
  
  # FIXME: If imaginary:
  p_m <- Mod(p_m)

  # Numerics can cause this deviates slightly from 1, normalise:
  if (!isTRUE(all.equal(sum(p_m), 1))) {
    warning("Probabilities summed to ", sum(p_m), " instead of 1, normalising was applied.")
  }
  p_m <- p_m / sum(p_m)
  
  res <- data.frame(d = (-meioses):meioses, p = p_m)
  
  return(res)
}

#' Calculate distribution of allele difference for symmetric mutation rates
#' 
#' Calculate distribution of allele difference after `m` meioses.
#' 
#' @inheritParams construct_M
#' @param mu_updw mutation rate for 1-step down- and up-mutations, i.e. total mutation rate is `2*mu_updw`
#' 
#' @return `data.frame` with columns `d` (allele difference) and `p` (prob)
#' 
#' @export
relationship_allele_diff_dist_sym <- function(meioses, mu_updw) {
  # parameters are validated in construct_M:
  M <- construct_M(meioses, mu_updw, mu_updw)

  # @param use_R_eigen Use `R`'s [eigen()] function or explicit calculation. Mostly for debugging.  
  use_R_eigen <- TRUE
  
  eigvals <- NULL
  P <- NULL
  
  if (use_R_eigen) {
    M_eig <- eigen(M, symmetric = TRUE)
    eigvals <- M_eig$values
    P <- M_eig$vectors
  } else {
    delta <- 1 - (mu_updw + mu_updw)
    sigma <- mu_updw
    tau <- mu_updw
    n <- 2*meioses + 1
    
    # (4) in http://www.math.kent.edu/~reichel/publications/toep3.pdf
    eigvals <- delta + 2*sqrt(sigma*tau)*cos(pi*seq_len(n)/(n+1))
    
    if (FALSE) {
      max(abs(eigvals - eigen(M, symmetric = TRUE)$values))
    }
    
    # Eigen vectors
    n <- 2*meioses + 1
    x <- seq_len(n)
    A <- outer(x, x, FUN = "*")
    B <- sin(pi*A / (n + 1))
    D <- apply(B, 2, function(x) x / norm(x, type = "2"))
    #D
    
    if (FALSE) {
      apply(B, 2, function(x) norm(x, type = "2"))
      max(abs(D - eigen(M, symmetric = TRUE)$vectors))
    }
    P <- D
  }
  # NOTE: Exploit tridiaginal Toepliz matrix instead
  
  D <- diag(eigvals^meioses)
  
  Pinv <- t(P)
  # Pinv <- t(P) only if mu_dw = mu_up because then is P orthogonal 
  # as all eigenvalues are distinct and then symmetric matrices 
  # have orthogonal eigenvectors
  
  res <- relationship_allele_diff_dist_sym_worker(meioses, P, D, Pinv)
  
  return(res)
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

  # @param use_R_eigen Use `R`'s [eigen()] function or explicit calculation. Mostly for debugging.  
  use_R_eigen <- TRUE

    
  eigvals <- NULL
  P <- NULL
  
  if (use_R_eigen) {
    M_eig <- eigen(M)
    eigvals <- M_eig$values
    P <- M_eig$vectors
  } else {
    delta <- 1 - (mu_dw + mu_up)
    sigma <- mu_up
    tau <- mu_dw
    n <- 2*meioses + 1
    
    # (4) in http://www.math.kent.edu/~reichel/publications/toep3.pdf
    eigvals <- delta + 2*sqrt(sigma*tau)*cos(pi*seq_len(n)/(n+1))
    
    # Eigen vectors?
  }
  # NOTE: Exploit tridiaginal Toepliz matrix instead
  
  D <- diag(eigvals^meioses)
  
  Pinv <- solve(P)
  # Pinv <- t(P) only if mu_dw = mu_up because then is P orthogonal 
  # as all eigenvalues are distinct and then symmetric matrices 
  # have orthogonal eigenvectors
  
  res <- relationship_allele_diff_dist_sym_worker(meioses, P, D, Pinv)

  return(res)
}
