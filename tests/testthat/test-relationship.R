context("Relationship")

test_that("construct_M()", {
  expect_error(construct_M(3, 0, 0))
  expect_error(construct_M(1, 1, 1))
  
  expect_equal(construct_M(1, 0.1, 0.2),
               matrix(c(0.7, 0.2, 0, 0.1, 0.7, 0.2, 0, 0.1, 0.7), 3, 3))
})


test_that("relationship_allele_diff_dist_sym()", {
  expect_error(relationship_allele_diff_dist_sym(1, 1))
  
  expect_equal(relationship_allele_diff_dist_sym(0, 0.01), data.frame(d = 0, p = 1))
  expect_equal(relationship_allele_diff_dist_sym(1, 0.01), 
               data.frame(d = (-1):1, p = c(0.01, 0.98, 0.01)))
  
  
  # Manually calculated by hand
  mu <- 3e-3
  p_0 <- ((1-mu)^2 + mu^2/2)
  p_1 <- (1-mu)*mu
  p_2 <- (1-(p_0+2*p_1))/2
  
  expect_equal(relationship_allele_diff_dist_sym(2, mu/2), 
               data.frame(d = (-2):2, p = c(p_2, p_1, p_0, p_1, p_2)))
  
  # Compare the two implementations
  for (m in c(0, 1, 2, 15, 20)) {
    expect_equal(
      relationship_allele_diff_dist_sym(m, mu/2), 
      relationship_allele_diff_dist_sym(m, mu/2, method = "explicit"),
      info = paste0("m = ", m))
    
    expect_equal(
      relationship_allele_diff_dist_sym(m, mu/2, method = "explicit"), 
      relationship_allele_diff_dist_sym(m, mu/2, method = "r_eigen"),
      info = paste0("m = ", m))
  }
})

test_that("relationship_allele_diff_dist()", {
  expect_error(relationship_allele_diff_dist(1, 1, 1))
  
  expect_equal(relationship_allele_diff_dist(0, 0.01, 0.01), data.frame(d = 0, p = 1))
  expect_equal(relationship_allele_diff_dist(1, 0.01, 0.01), 
               data.frame(d = (-1):1, p = c(0.01, 0.98, 0.01)))
  
  
  # Manually calculated by hand
  mu <- 3e-3
  p_0 <- ((1-mu)^2 + mu^2/2)
  p_1 <- (1-mu)*mu
  p_2 <- (1-(p_0+2*p_1))/2
  
  expect_equal(relationship_allele_diff_dist(2, mu/2, mu/2), 
               data.frame(d = (-2):2, p = c(p_2, p_1, p_0, p_1, p_2)))
  
  expect_equal(relationship_allele_diff_dist(2, mu/2, mu/2, method = "explicit"), 
               data.frame(d = (-2):2, p = c(p_2, p_1, p_0, p_1, p_2)))

  expect_equal(relationship_allele_diff_dist(2, mu/2, mu/2, method = "r_eigen"), 
               data.frame(d = (-2):2, p = c(p_2, p_1, p_0, p_1, p_2)))
  
  # Verify non-symmetry for mu_dw != mu_up:
  mu <- 3e-3
  p_10 <- relationship_allele_diff_dist(10, mu/4, 3*mu/4)
  p_10_1s <- subset(p_10, d == -1L | d == 1L)
  expect_true(!isTRUE(all.equal(p_10_1s$p[1], p_10_1s$p[2])))
  
  # Compare the two implementations
  for (m in c(0, 1, 2, 15, 20)) {
    expect_equal(
      relationship_allele_diff_dist(m, mu/4, 3*mu/4), 
      relationship_allele_diff_dist(m, mu/4, 3*mu/4, method = "explicit"),
      info = paste0("m = ", m))
    
    expect_equal(
      relationship_allele_diff_dist(m, mu/4, 3*mu/4, method = "explicit"), 
      relationship_allele_diff_dist(m, mu/4, 3*mu/4, method = "r_eigen"),
      info = paste0("m = ", m))
  }
})


test_that("relationship_allele_diff_dist_sym() methods", {
  meioses <- 5
  
  mu_events <- 3
  mu_meioses <- 1000
  mu_R <- mu_events / mu_meioses
  
  mu_dw1_R <- mu_R/2
  mu_up1_R <- mu_R/2
  
  expect_equal(relationship_allele_diff_dist_sym(meioses, mu_updw = mu_dw1_R), 
               relationship_allele_diff_dist_sym(meioses, mu_updw = mu_dw1_R, method = "explicit"))
  expect_equal(relationship_allele_diff_dist_sym(meioses, mu_updw = mu_dw1_R), 
               relationship_allele_diff_dist_sym(meioses, mu_updw = mu_dw1_R, method = "r_eigen"))
  expect_equal(relationship_allele_diff_dist_sym(meioses, mu_updw = mu_dw1_R), 
               relationship_allele_diff_dist_sym(meioses, mu_updw = mu_dw1_R, method = "matmult"))
  
  
  if (requireNamespace("Rmpfr", quietly = TRUE)) {
    prec <- 100
    mu_Rmpfr_events <- Rmpfr::mpfr(mu_events, precBits = prec)
    mu_Rmpfr_meioses <- Rmpfr::mpfr(mu_meioses, precBits = prec)
    mu_Rmpfr <- mu_Rmpfr_events / mu_Rmpfr_meioses
    
    mu_dw1_Rmpf <- mu_Rmpfr/2
    mu_up1_Rmpf <- mu_Rmpfr/2
    
    res <- relationship_allele_diff_dist_sym(meioses, mu_updw = mu_dw1_Rmpf, method = "matmult_mpfr")
    
    expect_equal(as.numeric(res$p), 
                 relationship_allele_diff_dist_sym(meioses, mu_updw = mu_dw1_R)$p)
  }
})


test_that("relationship_allele_diff_dist() methods", {
  meioses <- 5
  
  mu_events <- 3
  mu_meioses <- 1000
  mu_R <- mu_events / mu_meioses
  
  mu_dw1_R <- 3*mu_R/4
  mu_up1_R <- mu_R/4
  
  expect_equal(relationship_allele_diff_dist(meioses, mu_dw = mu_dw1_R, mu_up = mu_up1_R), 
               relationship_allele_diff_dist(meioses, mu_dw = mu_dw1_R, mu_up = mu_up1_R, method = "explicit"))
  expect_equal(relationship_allele_diff_dist(meioses, mu_dw = mu_dw1_R, mu_up = mu_up1_R), 
               relationship_allele_diff_dist(meioses, mu_dw = mu_dw1_R, mu_up = mu_up1_R, method = "r_eigen"))
  expect_equal(relationship_allele_diff_dist(meioses, mu_dw = mu_dw1_R, mu_up = mu_up1_R), 
               relationship_allele_diff_dist(meioses, mu_dw = mu_dw1_R, mu_up = mu_up1_R, method = "matmult"))
  
  expect_true(!isTRUE(all.equal(relationship_allele_diff_dist(meioses, mu_dw = mu_dw1_R, mu_up = mu_up1_R),
                                relationship_allele_diff_dist_sym(meioses, mu_updw = mu_dw1_R))))
  
  if (requireNamespace("Rmpfr", quietly = TRUE)) {
    prec <- 100
    mu_Rmpfr_events <- Rmpfr::mpfr(mu_events, precBits = prec)
    mu_Rmpfr_meioses <- Rmpfr::mpfr(mu_meioses, precBits = prec)
    mu_Rmpfr <- mu_Rmpfr_events / mu_Rmpfr_meioses
    
    mu_dw1_Rmpf <- 3*mu_Rmpfr/4
    mu_up1_Rmpf <- mu_Rmpfr/4
    
    res <- relationship_allele_diff_dist(meioses, 
                                         mu_dw = mu_dw1_Rmpf, 
                                         mu_up = mu_up1_Rmpf, 
                                         method = "matmult_mpfr")
    
    expect_equal(as.numeric(res$p), 
                 relationship_allele_diff_dist(meioses, mu_dw = mu_dw1_R, mu_up = mu_up1_R)$p)
  }
})
