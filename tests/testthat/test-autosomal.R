context("Autosomal")

ESTIMATION_TOL_DUE_TO_SAMPLING <- 0.01

set.seed(1)
sim_res_fixed <- sample_geneology(population_size = 1e3, 
                                  generations = 100, 
                                  extra_generations_full = 2,
                                  individuals_generations_return = 2, # default value
                                  progress = FALSE)

peds <- build_pedigrees(sim_res_fixed$population, progress = FALSE)

# alleles         A    B    C
allele_prob <- c(0.7, 0.2, 0.1)
theta <- 0.1
L <- length(allele_prob)

x <- replicate(10000, sample_autosomal_genotype(allele_prob, theta))
#prop.table(table(apply(x, 2, paste0, collapse = ";")))

x_p_tab <- prop.table(table(c(x)))
x_p <- as.numeric(x_p_tab)

test_that("sample_autosomal_genotype works", {
  expect_equal(x_p, allele_prob, tol = ESTIMATION_TOL_DUE_TO_SAMPLING)
})


geno_probs_R_mat <- matrix(0, nrow = L, ncol = L)  
for (i in 1L:L) {
  for (j in i:L) {
    #P_AA = F*p_A + (1 - F) * p_A^2
    #X != A: P_AX = (1 - F) * 2 * p_A*p_X
    prob <- if (i == j) # homzyg:
      theta*allele_prob[i] + (1-theta)*allele_prob[i]^2
    else
      (1-theta)*2*allele_prob[i]*allele_prob[j]
    
    # only assign (i, j);
    # if (j, i) should also be assigned, the factor 2 above in
    # the else case must be removed
    geno_probs_R_mat[i, j] <- prob
  }
}
geno_probs_R <- geno_probs_R_mat[upper.tri(geno_probs_R_mat, diag = TRUE)]
geno_probs_rcpp <- calc_autosomal_genotype_probs(allele_prob, theta)
test_that("calc_autosomal_genotype_probs works", {
  expect_equal(geno_probs_R, geno_probs_rcpp)
})

pedigrees_all_populate_autosomal(peds, allele_prob, theta, mutation_rate = 0, FALSE)

geno_mat <- do.call(rbind, lapply(seq_len(pedigrees_count(peds)), function(i) {
  hs <- do.call(rbind, get_haplotypes_in_pedigree(peds[[i]]))
  hs
}))
geno_vec <- c(geno_mat)
geno_vec_p <- as.numeric(prop.table(table(geno_vec)))

test_that("pedigrees_all_populate_autosomal works", {
  expect_equal(geno_vec_p, allele_prob, tol = ESTIMATION_TOL_DUE_TO_SAMPLING)
})



get_ols_quantities <- function(y) {
  x_unique_geno <- y[!duplicated(y), ]
  ols_quantities <- apply(x_unique_geno, 1,
                          function(as) {
                            geno_prob <- mean(apply(y, 1, function(as2) all(as == as2)))
                            
                            if (length(unique(as)) == 1L) {
                              p <- x_p_tab[as.character(as[1])]
                              
                              x <- p - p^2
                              y <- geno_prob - p^2
                              return(c(x, y))
                            } else {
                              p <- x_p_tab[as.character(as[1])]
                              q <- x_p_tab[as.character(as[2])]
                              
                              x <- -2*p*q
                              y <- geno_prob - 2*p*q
                              return(c(x, y))
                            }
                          })
  ols_x <- as.matrix(ols_quantities[1L, ])
  ols_y <- ols_quantities[2L, ]
  return(list(x = ols_x, y = ols_y))
}

y <- t(x)
ols_res_sample <- get_ols_quantities(y)
theta_res <- estimate_theta_1subpop_sample(y)
test_that("estimate_theta_1subpop_sample works", {
  expect_true(theta_res$error == FALSE)
  expect_equal(qr.solve(ols_res_sample$x, ols_res_sample$y), 
               theta_res$estimate)
})

theta_res_info <- estimate_theta_1subpop_sample(y, return_estimation_info = TRUE)
test_that("estimate_theta_1subpop_sample return_estimation_info works", {
  expect_true(theta_res_info$error == FALSE)
  expect_equal(qr.solve(ols_res_sample$x, ols_res_sample$y), 
               theta_res_info$estimate)
  expect_equal(sort(ols_res_sample$x), sort(theta_res_info$estimation_info$X))
})

# bootstrap:
theta_boot <- replicate(100, {
  yboot <- y[sample(nrow(y), replace = TRUE), ]
  estimate_theta_1subpop_sample(yboot)$estimate
})
# Expanding a bit...
theta_boot_rng <- c(0.9, 1.1) * range(theta_boot)
#theta_boot_rng
test_that("estimate_theta_1subpop_sample boot contains true", {
  expect_true(all(theta_boot >= 0 & theta_boot <= 1))
  expect_true(theta >= theta_boot_rng[1L] & theta <= theta_boot_rng[2L])
})


p_cumdists <- calc_autosomal_genotype_conditional_cumdist(allele_dist = allele_prob, 
                                                          theta = theta)

cumdist_r <- geno_probs_R_mat
cumdist_r[lower.tri(cumdist_r)] <- cumdist_r[upper.tri(cumdist_r)] / 2
cumdist_r[upper.tri(cumdist_r)] <- cumdist_r[upper.tri(cumdist_r)] / 2
cumdist_r <- t(apply(cumdist_r / rowSums(cumdist_r), 1, cumsum))
#cumdist_r

test_that("calc_autosomal_genotype_conditional_cumdist works", {
  expect_equal(p_cumdists, cumdist_r)
})


livepop <- sim_res_fixed$individuals_generations
y_livepop <- get_haplotypes_individuals(livepop)
ols_res_livepop <- get_ols_quantities(y_livepop)

test_that("estimate_theta_1subpop_individuals works", {
  expect_equal(qr.solve(ols_res_livepop$x, ols_res_livepop$y),
               estimate_theta_1subpop_individuals(livepop)$estimate, 
               tol = 10*ESTIMATION_TOL_DUE_TO_SAMPLING
               )
})

test_that("livepop: haplotypes/individual same theta", {
  expect_equal(estimate_theta_1subpop_sample(y_livepop)$estimate,
               estimate_theta_1subpop_individuals(livepop)$estimate)
})


test_that("geneaology independent sample and population same theta", {
  expect_equal(estimate_theta_1subpop_sample(y)$estimate,
               estimate_theta_1subpop_sample(y_livepop)$estimate,
               tol = ESTIMATION_TOL_DUE_TO_SAMPLING)
  
  expect_equal(estimate_theta_1subpop_sample(y)$estimate,
               estimate_theta_1subpop_individuals(livepop)$estimate, 
               tol = ESTIMATION_TOL_DUE_TO_SAMPLING)
})


if (FALSE) {
  prop.table(table(c(x)))
  prop.table(table(c(y_livepop)))
  
  prop.table(table(apply(x, 2, paste0, collapse = ";")))
  prop.table(table(apply(y_livepop, 1, paste0, collapse = ";")))
  
  estimate_theta_1subpop_sample(y)
  qr.solve(ols_res_sample$x, ols_res_sample$y)
  
  estimate_theta_1subpop_sample(y_livepop)
  estimate_theta_1subpop_individuals(livepop)
  qr.solve(ols_res_livepop$x, ols_res_livepop$y)
}

