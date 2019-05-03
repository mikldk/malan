context("Relationship")

test_that("construct_M()", {
  expect_equal(construct_M(3, 0, 0), diag(2*3 + 1))
  expect_error(construct_M(1, 1, 1))
  
  expect_equal(construct_M(1, 0.1, 0.2),
               matrix(c(0.7, 0.2, 0, 0.1, 0.7, 0.2, 0, 0.1, 0.7), 3, 3))
})


test_that("relationship_allele_diff_dist_sym()", {
  expect_error(relationship_allele_diff_dist_sym(1, 1))
  
  expect_equal(relationship_allele_diff_dist_sym(0, 0), data.frame(d = 0, p = 1))
  expect_equal(relationship_allele_diff_dist_sym(1, 0), 
               data.frame(d = (-1):1, p = c(0, 1, 0)))
  
  
  # Manually calculated by hand
  mu <- 3e-3
  p_0 <- ((1-mu)^2 + mu^2/2)
  p_1 <- (1-mu)*mu
  p_2 <- (1-(p_0+2*p_1))/2
  
  expect_equal(relationship_allele_diff_dist_sym(2, mu/2), 
               data.frame(d = (-2):2, p = c(p_2, p_1, p_0, p_1, p_2)))
})

test_that("relationship_allele_diff_dist()", {
  expect_error(relationship_allele_diff_dist(1, 1, 1))
  
  expect_equal(relationship_allele_diff_dist(0, 0, 0), data.frame(d = 0, p = 1))
  expect_equal(relationship_allele_diff_dist(1, 0, 0), 
               data.frame(d = (-1):1, p = c(0, 1, 0)))
  
  
  # Manually calculated by hand
  mu <- 3e-3
  p_0 <- ((1-mu)^2 + mu^2/2)
  p_1 <- (1-mu)*mu
  p_2 <- (1-(p_0+2*p_1))/2
  
  expect_equal(relationship_allele_diff_dist(2, mu/2, mu/2), 
               data.frame(d = (-2):2, p = c(p_2, p_1, p_0, p_1, p_2)))
  
  
  # Verify non-symmetry for mu_dw != mu_up:
  mu <- 3e-3
  p_100 <- relationship_allele_diff_dist(100, mu/4, 3*mu/4)
  p_100_1s <- subset(p_100, d == -1L | d == 1L)
  expect_true(!isTRUE(all.equal(p_100_1s$p[1], p_100_1s$p[2])))
})
