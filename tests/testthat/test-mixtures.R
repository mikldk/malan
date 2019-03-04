context("Mixtures")

set.seed(1)
sim_res_growth <- sample_geneology_varying_size(population_sizes = rep(1e3, 200),
                                                generations_full = 3,
                                                generations_return = 3,
                                                enable_gamma_variance_extension = TRUE,
                                                gamma_parameter_shape = 5,
                                                gamma_parameter_scale = 1/5,
                                                progress = FALSE)

pedigrees <- build_pedigrees(sim_res_growth$population, progress = FALSE)

mutrts <- rep(0.001, 20)
pedigrees_all_populate_haplotypes(pedigrees = pedigrees, 
                                  loci = length(mutrts), 
                                  mutation_rates = mutrts, progress = FALSE)
live_individuals <- sim_res_growth$individuals_generations

U_indices <- sample.int(n = length(live_individuals), size = 2, replace = FALSE)
U1 <- live_individuals[[U_indices[1]]]
U2 <- live_individuals[[U_indices[2]]]
H1 <- get_haplotype(U1)
H2 <- get_haplotype(U2)

# Max 100 redraws
for (i in 1:100) {
  if (any(H1 != H2)) {
    break
  }
  
  # If H1 == H2, try again:
  U_indices <- sample.int(n = length(live_individuals), size = 2, replace = FALSE)
  U1 <- live_individuals[[U_indices[1]]]
  U2 <- live_individuals[[U_indices[2]]]
  H1 <- get_haplotype(U1)
  H2 <- get_haplotype(U2)
}

test_that("contributors different", {
  expect_true(any(H1 != H2))
})


mixres <- mixture_info_by_individuals_2pers(live_individuals, U1, U2)
others_haps <- get_haplotypes_pids(sim_res_growth$population, mixres$pids_others_included)

others_indv <- lapply(mixres$pids_others_included, function(pid) {
  get_individual(sim_res_growth$population, pid)
})
others_haps_2 <- get_haplotypes_individuals(individuals = others_indv)

test_that("mixture_info_by_individuals_2pers works", {
  expect_equal(mixres$donor1_profile, H1)
  expect_equal(mixres$donor2_profile, H2)
  expect_equal(mixres$loci_not_matching, sum(H1 != H2))
  expect_equal(length(mixres$pids_matching_donor1),
               count_haplotype_occurrences_individuals(individuals = live_individuals, haplotype = H1))
  expect_equal(length(mixres$pids_matching_donor2),
               count_haplotype_occurrences_individuals(individuals = live_individuals, haplotype = H2))
  expect_equal(others_haps, others_haps_2)
})

others_haps_unique <- others_haps[!duplicated(others_haps), ]
others_haps_counts <- unlist(lapply(seq_len(nrow(others_haps_unique)), function(hap_i) {
  count_haplotype_occurrences_individuals(individuals = live_individuals, 
                                          haplotype = others_haps_unique[hap_i, ])
}))

test_that("mixture others included haplotypes works", {
  expect_equal(sum(others_haps_counts), length(mixres$pids_others_included))
})

##############################################################################

set.seed(1)
sim_res_growth <- sample_geneology_varying_size(population_sizes = rep(1e3, 200),
                                                generations_full = 3,
                                                generations_return = 3,
                                                enable_gamma_variance_extension = FALSE,
                                                progress = FALSE)

pedigrees <- build_pedigrees(sim_res_growth$population, progress = FALSE)

mutrts <- rep(0.1, 2)
pedigrees_all_populate_haplotypes_ladder_bounded(pedigrees = pedigrees, 
                                                 mutation_rates = mutrts, 
                                                 ladder_min = c(1, 1),
                                                 ladder_max = c(2, 3),
                                                 get_founder_haplotype = function() c(1, 1),
                                                 progress = FALSE)

live_individuals <- sim_res_growth$individuals_generations
indvs_haps <- get_haplotypes_individuals(live_individuals)

U1 <- live_individuals[[ which(apply(indvs_haps, 1, function(x) identical(x, c(1L, 1L))))[1L] ]]
U2 <- live_individuals[[ which(apply(indvs_haps, 1, function(x) identical(x, c(2L, 2L))))[1L] ]]
U3 <- live_individuals[[ which(apply(indvs_haps, 1, function(x) identical(x, c(2L, 3L))))[1L] ]]


###############################################################################################
# 2 persons
###############################################################################################
# mix = U1 (1, 1) + U2 (2, 2) = {1, 2} x {1, 2}
mixres2 <- mixture_info_by_individuals_2pers(live_individuals, U1, U2)

# Haps included in mixture
haps2_unique <- list(c(1L, 1L),  # <- U1 = SUSPECT / KNOWN
                     c(1L, 2L), 
                     c(2L, 1L), 
                     c(2L, 2L))
haps2_unique_counts <- c(11, 12, 21, 22) # counts refering to types
res2 <- analyse_mixture_result(mixres2, haps2_unique, haps2_unique_counts)

test_that("mixture 2 persons", {
  # Hp: only (2, 2) with count 22 must be explained:
  expect_equal(res2$terms_Hp, list(22))
  
  # Hd: mixture either (1, 1)+(2, 2) or (1, 2)+(2, 1) [(2, 1) + (1, 2) is same but considering order]
  # Two random males would give:
  #   expect_equal(res2$terms_Hd, list(c(11, 22), c(12, 21)))
  # But now we say two random different from the PoI, so substract by 1 for PoI's profile (11):
  expect_equal(res2$terms_Hd, list(c(11-1, 22), c(12, 21)))
})


###############################################################################################
# 3 persons
###############################################################################################
# mix = U1 (1, 1) + U2 (2, 2) + U3 (2, 3) = {1, 2} x {1, 2, 3}
mixres3 <- mixture_info_by_individuals_3pers(live_individuals, U1, U2, U3)

if (FALSE) {
  haps_in_mixture_3pers <- get_haplotypes_pids(sim_res_growth$population, mixres3$pids_included_in_mixture)
  haps_in_mixture_list_3pers <- lapply(seq_len(nrow(haps_in_mixture_3pers)), function(i) haps_in_mixture_3pers[i, , drop = FALSE])
  unique_haps_in_mixture_3pers <- haps_in_mixture_list_3pers[!duplicated(haps_in_mixture_3pers)]
  hapcount_pop_3gen_3pers <- unlist(lapply(seq_along(unique_haps_in_mixture_3pers), function(hap_i) {
    count_haplotype_occurrences_individuals(individuals = live_individuals, 
                                            haplotype = unique_haps_in_mixture_3pers[[hap_i]])
  }))
  do.call(rbind, unique_haps_in_mixture_3pers)
}

# Haps included in mixture, pop made such that this is {1, 2} x {1, 2, 3}
haps3_unique <- list(c(1L, 1L),  # <- U1 = SUSPECT / KNOWN
                     c(1L, 2L), 
                     c(1L, 3L),
                     c(2L, 1L), 
                     c(2L, 2L),
                     c(2L, 3L))
haps3_unique_counts <- c(11, 12, 13, 21, 22, 23) # counts refering to types
res3 <- analyse_mixture_result(mixres3, haps3_unique, haps3_unique_counts)

test_that("mixture 3 persons", {
  # Hp: allele 2 at locus 1 and alleles {2, 3} at locus 2 must be explained.
  #     thus unknowns are bound to have 2/3 at locus 2, but first locus can change, 
  #     at least one must have 2 at locus 1, but both can potentially:
  #     (L1, L2)
  #     (1, 2) + (2, 3)
  #     (2, 2) + (1, 3)
  #     (2, 2) + (2, 3)
  expect_equal(res3$terms_Hp, list(
    c(12, 23), 
    c(13, 22),
    c(22, 23)
    ))
  
  # Hd: 
  #     (L1, L2)
  #     (1, 1) + (1, 2) + (2, 3)
  #     (1, 1) + (2, 2) + (2, 3)
  #     (1, 1) + (2, 2) + (1, 3) -> [enumeration order] -> (1, 1) + (1, 3) + (2, 2)
  # 
  #     (1, 2) + (1, 3) + (2, 1)
  #     (1, 2) + (2, 1) + (2, 3)
  # 
  #     (1, 3) + (2, 1) + (2, 2)
  
  # Three random males would give:
  #   expect_equal(res3$terms_Hd, list(
  #     c(11, 12, 23), 
  #     c(11, 13, 22), # (1, 1) + (2, 2) + (1, 3) -> [enumeration order] -> (1, 1) + (1, 3) + (2, 2)
  #     c(11, 22, 23),
  #   
  #     c(12, 13, 21),
  #     c(12, 21, 23),
  #   
  #     c(13, 21, 22)
  #   ))
  #   expect_equal(res2$terms_Hd, list(c(11, 22), c(12, 21)))
  # But now we say three random different from the PoI, so substract by 1 for PoI's profile (11):
  expect_equal(res3$terms_Hd, list(
    c(11-1, 12, 23), 
    c(11-1, 13, 22), # (1, 1) + (2, 2) + (1, 3) -> [enumeration order] -> (1, 1) + (1, 3) + (2, 2)
    c(11-1, 22, 23),
    
    c(12, 13, 21),
    c(12, 21, 23),
    
    c(13, 21, 22)
  ))
})


