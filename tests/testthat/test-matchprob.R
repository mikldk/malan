context("Match probability")

# A. Caliebe and M. Krawczak (2018): 
# "Match probabilities for Y-chromosomal profiles: A paradigm shift"
# 
# P(match) \approx \sum_{m = 1}^\infty 
# (1-\mu)^m \frac{P(M = m)}{P(M > 0)}


set.seed(1)
sim_res_growth <- sample_geneology_varying_size(population_sizes = rep(1e3, 200),
                                                generations_full = 3,
                                                generations_return = 3,
                                                enable_gamma_variance_extension = TRUE,
                                                gamma_parameter_shape = 5,
                                                gamma_parameter_scale = 1/5,
                                                progress = FALSE)

peds <- build_pedigrees(sim_res_growth$population, progress = FALSE)
#peds

# # ped index 0 is largest:
# ped_biggest <- get_pedigree(peds, 0L)
# # sapply(seq_len(pedigrees_count(peds))-1L, function(i) {
# #   pedigree_size(get_pedigree(peds, i))
# # })
# 
# 
# ped_indvs <- pedi)

# FIXME: Get only pedigree?

live_individuals <- sim_res_growth$individuals_generations
indv_poi <- live_individuals[[sample.int(n = length(live_individuals), size = 1)]]
transfer_events <- sapply(live_individuals, meiotic_dist, ind2 = indv_poi)

# ...And not suspect:
in_ped <- which(transfer_events >= 1L)
suspect_pop <- live_individuals[in_ped]

M_set <- as.data.frame(table(transfer_events[in_ped]))
M_set$Var1 <- as.integer(as.character(M_set$Var1))
N <- length(in_ped)
M_set$p <- M_set$Freq / N 

test_that("suspect pop", {
  expect_equal(sum(M_set$p), 1)
})

############# 

# P(match) = \sum{m = 1}^\infty (1-\mu)^m \frac{P(M = m)}{P(M > 0)}

mutrts <- rep(0.001, 20)
mu <- 1 - prod(1 - mutrts) # approx 20 * 0.001

matchprob <- sum(((1 - mu)^M_set$Var1) * M_set$p)
matchprob
lr <- 1/matchprob
lr

