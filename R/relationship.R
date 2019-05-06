verify_input_meioses_mu_dw_mu_up <- function(meioses, mu_dw, mu_up) {
  stopifnot(length(meioses) == 1L)
  stopifnot(isTRUE(all.equal(round(meioses), meioses)))
  meioses <- as.integer(meioses)
  stopifnot(meioses >= 0L)
  
  stopifnot(length(mu_dw) == 1L)
  stopifnot(mu_dw > 0 && mu_dw <= 1)
  
  stopifnot(length(mu_up) == 1L)
  stopifnot(mu_up > 0 && mu_up <= 1)
  
  stopifnot(mu_dw + mu_up <= 1)
}

#' Construct M matrix
#' 
#' @param meioses number of meioses separating the two individuals
#' @param mu_dw mutation rate for 1-step down-mutation
#' @param mu_up mutation rate for 1-step up-mutation
#' 
#' @export
construct_M <- function(meioses, mu_dw, mu_up) {
  verify_input_meioses_mu_dw_mu_up(meioses, mu_dw, mu_up)
  
  meioses <- as.integer(meioses)
  
  M <- diag(x = 1 - (mu_up + mu_dw), nrow = 2*meioses + 1)
  
  for (i in seq_len(2*meioses)) {
    M[i, i+1] <- mu_dw
    M[i+1, i] <- mu_up
  }
  
  return(M)
}

construct_M_mpfr <- function(meioses, mu_dw, mu_up) {
  if (!requireNamespace("Rmpfr", quietly = TRUE)) {
    stop("Rmpfr is required to use matmult_mpfr method")
  }
  
  if (!is(mu_dw, "mpfr")) {
    stop("mu_dw must be a mpfr")
  }
  
  if (!is(mu_up, "mpfr")) {
    stop("mu_up must be a mpfr")
  }
  
  verify_input_meioses_mu_dw_mu_up(meioses, mu_dw, mu_up)
  meioses <- as.integer(meioses)
  
  not_mut <- 1 - (mu_up + mu_dw)
  
  M <- Rmpfr::mpfr(x = rep(0, (2*meioses + 1)*(2*meioses + 1)), 
                   precBits = Rmpfr::getPrec(not_mut))
  dim(M) <- c(2*meioses + 1, 2*meioses + 1)
  diag(M) <- not_mut

  for (i in seq_len(2*meioses)) {
    M[i, i+1] <- mu_dw
    M[i+1, i] <- mu_up
  }
  
  return(M)
}

relationship_allele_diff_dist_matmult_worker <- function(meioses, M, p_0) {
  p_m <- p_0
  
  if (meioses > 0L) {
    M_m <- M
    
    if (meioses >= 2L) {
      # FIXME: Exploit squaring etc.: x^5 = x*x^4 = x*(x^2)^2
      for (i in 2L:meioses) {
        M_m <- M_m %*% M
      }
    }
    
    p_m <- M_m %*% p_0
  }
  
  return(p_m)
}

relationship_allele_diff_dist_matmult <- function(meioses, M) {
  p_0 <- matrix(0, nrow = 2*meioses + 1, ncol = 1)
  p_0[meioses + 1] <- 1
  
  p_m <- relationship_allele_diff_dist_matmult_worker(meioses, M, p_0)
  
  res <- data.frame(d = (-meioses):meioses, p = p_m)
  return(res)
}

relationship_allele_diff_dist_matmult_Rmpfr <- function(meioses, M, prec_bits) {
  p_0 <- Rmpfr::mpfr(rep(0, 2*meioses + 1), precBits = prec_bits)
  dim(p_0) <- c(2*meioses + 1, 1)
  p_0[meioses + 1] <- 1
  
  p_m <- relationship_allele_diff_dist_matmult_worker(meioses, M, p_0)
  dim(p_m) <- NULL
  
  d <- (-meioses):meioses
  res <- list(d = d, p = p_m)
  
  return(res)
}

relationship_allele_diff_dist_worker <- function(meioses, P, D, Pinv) {
  M_m <- P %*% D %*% Pinv
  p_0 <- matrix(0, nrow = 2*meioses + 1, ncol = 1)
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
#' @param method `"explicit"` (default): use known formulas for eigenvalues and eigenvectors. 
#' Can cause numerical problems.
#' `"matmult"`: do matrix multiplication instead of diagonalisation.
#' `"matmult_mpfr"`: as `"matmult"` but with the `Rmpfr` library (note that this returns `list` instead of `data.frame`).
#' `"r_eigen"`: use `R`'s [eigen()] function to find eigen values. Mostly for debugging.
#' 
#' @return `data.frame` with columns `d` (allele difference) and `p` (prob)
#' 
#' @export
relationship_allele_diff_dist_sym <- function(meioses, mu_updw, 
                                              method = "explicit") {
  if (!(method %in% c("explicit", "matmult", "matmult_mpfr", "r_eigen"))) {
    stop("Invalid method")
  }
  
  verify_input_meioses_mu_dw_mu_up(meioses, mu_updw, mu_updw)
  meioses <- as.integer(meioses)

  if (method == "matmult") {
    M <- construct_M(meioses, mu_updw, mu_updw)
    res <- relationship_allele_diff_dist_matmult(meioses = meioses, M = M)
    return(res)
  } else if (method == "matmult_mpfr") {
    M_mfpr <- construct_M_mpfr(meioses = meioses, mu_dw = mu_updw, mu_up = mu_updw)
    
    res <- relationship_allele_diff_dist_matmult_Rmpfr(meioses = meioses, 
                                                       M = M_mfpr, 
                                                       prec_bits = Rmpfr::getPrec(mu_updw))
    return(res)
  }
  
  M <- construct_M(meioses, mu_updw, mu_updw)
  
  eigvals <- NULL
  P <- NULL
  
  if (method == "r_eigen") {
    M_eig <- eigen(M, symmetric = TRUE)
    eigvals <- M_eig$values
    P <- M_eig$vectors
  } else if (method == "explicit") {
    delta <- 1 - (mu_updw + mu_updw)
    n <- 2*meioses + 1
    
    # (4)/(7) in http://www.math.kent.edu/~reichel/publications/toep3.pdf
    eigvals <- delta + 2*mu_updw*cos(pi*seq_len(n)/(n+1))
    
    if (FALSE) {
      max(abs(eigvals - eigen(M, symmetric = TRUE)$values))
    }
    
    # Eigen vectors
    n <- 2*meioses + 1
    D <- matrix(0, n, n)
    for (h in seq_len(n)) {
      for (k in seq_len(n)) {
        D[k, h] <- sin(h*k*pi / (n + 1))
      }
      D[, h] <- D[, h] / norm(D[, h], type = "2")
    }
    
    if (FALSE) {
      diffs <- unlist(lapply(seq_len(n), function(i) {
        M %*% D[, i] - eigvals[i] * D[, i]
      }))
      max(abs(diffs))
    }
    
    P <- D
  } else {
    stop("Unexpected")
  }
  
  D <- diag(eigvals^meioses)
  
  Pinv <- t(P)
  # Pinv <- t(P) only if mu_dw = mu_up because then is P orthogonal 
  # as all eigenvalues are distinct and then symmetric matrices 
  # have orthogonal eigenvectors
  
  res <- relationship_allele_diff_dist_worker(meioses, P, D, Pinv)
  
  return(res)
}


#' Calculate distribution of allele difference
#' 
#' Calculate distribution of allele difference after `m` meioses.
#' 
#' @inheritParams construct_M
#' @param method `"explicit"` (default): use known formulas for eigenvalues and eigenvectors. 
#' Can cause numerical problems.
#' `"matmult"`: do matrix multiplication instead of diagonalisation.
#' `"matmult_mpfr"`: as `"matmult"` but with the `Rmpfr` library (note that this returns `list` instead of `data.frame`).
#' `"r_eigen"`: use `R`'s [eigen()] function to find eigen values. Mostly for debugging.
#' 
#' @return `data.frame` with columns `d` (allele difference) and `p` (prob)
#' 
#' @export
relationship_allele_diff_dist <- function(meioses, mu_dw, mu_up,
                                          method = "explicit") {
  if (!(method %in% c("explicit", "matmult", "matmult_mpfr", "r_eigen"))) {
    stop("Invalid method")
  }
  
  verify_input_meioses_mu_dw_mu_up(meioses, mu_dw, mu_up)
  meioses <- as.integer(meioses)
  
  if (method == "matmult") {
    M <- construct_M(meioses, mu_dw, mu_up)
    res <- relationship_allele_diff_dist_matmult(meioses = meioses, M = M)
    return(res)
  } else if (method == "matmult_mpfr") {
    if (!requireNamespace("Rmpfr", quietly = TRUE)) {
      stop("Rmpfr is required to use matmult_mpfr method")
    }
    
    if (Rmpfr::getPrec(mu_dw) != Rmpfr::getPrec(mu_up)) {
      stop("Please choose same precision for both mu_dw and mu_up")
    }
    
    # In principle stupid to construct again, but really cheap
    # and construct_M() verifies arguments, so doesn't matter
    M_mfpr <- construct_M_mpfr(meioses = meioses, mu_dw = mu_dw, mu_up = mu_up)
    
    res <- relationship_allele_diff_dist_matmult_Rmpfr(meioses = meioses, 
                                                       M = M_mfpr, 
                                                       prec_bits = Rmpfr::getPrec(mu_dw))
    return(res)
  }

  M <- construct_M(meioses, mu_dw, mu_up)
  
  eigvals <- NULL
  P <- NULL
  
  if (method == "r_eigen") {
    M_eig <- eigen(M)
    eigvals <- M_eig$values
    P <- M_eig$vectors
  } else if (method == "explicit") {
    delta <- 1 - (mu_dw + mu_up)
    tau <- mu_dw
    sigma <- mu_up
    n <- 2*meioses + 1
    
    # (4)/(7) in http://www.math.kent.edu/~reichel/publications/toep3.pdf
    eigvals <- delta + 2*sqrt(sigma*tau)*cos(pi*seq_len(n)/(n+1))
    
    if (FALSE) {
      eigvals
      eigen(M)$values
      max(abs(eigvals - eigen(M)$values))
    }
    
    # Eigen vectors
    n <- 2*meioses + 1
    D <- matrix(0, n, n)
    for (h in seq_len(n)) {
      for (k in seq_len(n)) {
        D[k, h] <- ((sigma / tau)^(k/2)) * sin(h*k*pi / (n + 1))
      }
      D[, h] <- D[, h] / norm(D[, h], type = "2")
    }
    
    if (FALSE) {
      diffs <- unlist(lapply(seq_len(n), function(i) {
        M %*% D[, i] - eigvals[i] * D[, i]
      }))
      max(abs(diffs))
    }
    
    P <- D
  } else {
    stop("Unexpected")
  }

  D <- diag(eigvals^meioses)
  
  Pinv <- solve(P)
  # Pinv <- t(P) only if mu_dw = mu_up because then is P orthogonal 
  # as all eigenvalues are distinct and then symmetric matrices 
  # have orthogonal eigenvectors
  
  res <- relationship_allele_diff_dist_worker(meioses, P, D, Pinv)

  return(res)
}
