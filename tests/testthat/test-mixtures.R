context("Mixtures")

set.seed(1)
sim_res_growth <- sample_geneology_varying_size(population_sizes = rep(1e3, 200),
                                                enable_gamma_variance_extension = TRUE,
                                                gamma_parameter_shape = 5,
                                                gamma_parameter_scale = 1/5,
                                                generations_full = 3,
                                                individuals_generations_return = 2, # default value
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

test_that("contributors differet", {
  expect_true(any(H1 != H2))
})


mixres <- mixture_info_by_individuals(live_individuals, U1, U2)
others_haps <- get_haplotypes_pids(sim_res_growth$population, mixres$pids_others_included)

others_indv <- lapply(mixres$pids_others_included, function(pid) {
  get_individual(sim_res_growth$population, pid)
})
others_haps_2 <- get_haplotypes_individuals(individuals = others_indv)

test_that("mixture_info_by_individuals works", {
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
