## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(fig.width = 7)

## ---- message=FALSE------------------------------------------------------
library(malan)

## ------------------------------------------------------------------------
set.seed(1)

## ------------------------------------------------------------------------
sim_res <- sample_geneology(population_size = 10, generations = 10, progress = FALSE)

## ------------------------------------------------------------------------
pedigrees <- build_pedigrees(sim_res$population, progress = FALSE)
pedigrees
pedigrees_count(pedigrees)
pedigrees_table(pedigrees)
pedigree_size(pedigrees[[1]])
pedigree_size(pedigrees[[2]])
#pedigree_size(pedigrees[[3]]) # error as there are only 2 pedigrees

## ------------------------------------------------------------------------
plot(pedigrees)

## ------------------------------------------------------------------------
plot(pedigrees[[1]])

## ------------------------------------------------------------------------
plot(pedigrees[[2]])

## ------------------------------------------------------------------------
str(sim_res, 1)
live_individuals <- sim_res$end_generation_individuals

## ------------------------------------------------------------------------
print_individual(live_individuals[[1]])

## ------------------------------------------------------------------------
indv <- get_individual(sim_res$population, 22)
print_individual(indv)

## ------------------------------------------------------------------------
set.seed(1)

mutrts <- c(0.5, 0.5)
pedigrees_all_populate_haplotypes(pedigrees = pedigrees, 
                                  loci = length(mutrts), 
                                  mutation_rates = mutrts, progress = FALSE)

## ------------------------------------------------------------------------
plot(pedigrees[[1]], haplotypes = TRUE)

## ------------------------------------------------------------------------
plot(pedigrees[[1]], ids = FALSE, haplotypes = TRUE)

## ------------------------------------------------------------------------
plot(pedigrees[[1]], ids = TRUE, haplotypes = TRUE, mark_pids = c(14, 16))

## ------------------------------------------------------------------------
set.seed(1)
sim_res <- sample_geneology(population_size = 10, 
                            generations = 5, 
                            extra_generations_full = 2,
                            progress = FALSE)
pedigrees <- build_pedigrees(sim_res$population, progress = FALSE)
plot(pedigrees)

## ------------------------------------------------------------------------
set.seed(1)
sim_res <- sample_geneology(population_size = 10, 
                            generations = 5, 
                            extra_generations_full = 10,
                            progress = FALSE)
pedigrees <- build_pedigrees(sim_res$population, progress = FALSE)
plot(pedigrees)

## ------------------------------------------------------------------------
set.seed(1)
sim_res <- sample_geneology(population_size = 10, 
                            generations = -1, 
                            progress = FALSE)
pedigrees <- build_pedigrees(sim_res$population, progress = FALSE)
plot(pedigrees)

## ------------------------------------------------------------------------
sim_res$generations

## ------------------------------------------------------------------------
set.seed(1)
sim_res <- sample_geneology(population_size = 1e3, 
                            generations = 200, 
                            extra_generations_full = 2,
                            individuals_generations_return = 2, # default value
                            progress = FALSE)

## ------------------------------------------------------------------------
pedigrees <- build_pedigrees(sim_res$population, progress = FALSE)
pedigrees_table(pedigrees)
pedigrees_count(pedigrees)

## ------------------------------------------------------------------------
ped_sizes <- sapply(1L:pedigrees_count(pedigrees), function(i) pedigree_size(pedigrees[[i]]))
ped_sizes
largest_i <- which.max(ped_sizes)
plot(pedigrees[[largest_i]])

## ------------------------------------------------------------------------
set.seed(10)
mutrts <- rep(0.001, 20)
pedigrees_all_populate_haplotypes(pedigrees = pedigrees, 
                                  loci = length(mutrts), 
                                  mutation_rates = mutrts, progress = FALSE)

## ------------------------------------------------------------------------
live_individuals <- sim_res$individuals_generations
length(live_individuals)

## ------------------------------------------------------------------------
haps <- get_haplotypes_individuals(individuals = live_individuals)

## ------------------------------------------------------------------------
head(haps)

## ------------------------------------------------------------------------
haps_str <- apply(haps, 1, paste0, collapse = ";")
haps_tab <- table(haps_str)
sort(haps_tab, decreasing = TRUE)[1:10]
spectrum <- table(haps_tab)
spectrum

## ------------------------------------------------------------------------
set.seed(100)
Q_index <- sample.int(n = length(live_individuals), size = 1)
Q <- live_individuals[[Q_index]]
Q_hap <- get_haplotype(Q)
Q_hap

## ------------------------------------------------------------------------
Q_ped <- get_pedigree_from_individual(Q)

## ------------------------------------------------------------------------
count_haplotype_occurrences_pedigree(pedigree = Q_ped, haplotype = Q_hap, generation_upper_bound_in_result = 2)
count_haplotype_occurrences_individuals(individuals = live_individuals, haplotype = Q_hap)

## ------------------------------------------------------------------------
path_details <- pedigree_haplotype_matches_in_pedigree_meiosis_L1_dists(suspect = Q, 
                                                                        generation_upper_bound_in_result = 2)

## ------------------------------------------------------------------------
nrow(path_details)
head(path_details)

## ------------------------------------------------------------------------
meioses <- path_details[, 1]
hist(meioses)

## ------------------------------------------------------------------------
L1_max <- path_details[, 2]
table(L1_max)
mean(L1_max == 0)

## ------------------------------------------------------------------------
set.seed(100)
U_indices <- sample.int(n = length(live_individuals), size = 2, replace = FALSE)
U1 <- live_individuals[[U_indices[1]]]
U2 <- live_individuals[[U_indices[2]]]
H1 <- get_haplotype(U1)
H2 <- get_haplotype(U2)

## ------------------------------------------------------------------------
rbind(H1, H2)

## ------------------------------------------------------------------------
#mixres <- indices_in_mixture_by_haplotype_matrix(haplotypes = haps, H1 = H1, H2 = H2)
mixres <- mixture_info_by_individuals(live_individuals, U1, U2)
str(mixres, 1)

## ------------------------------------------------------------------------
length(mixres$pids_matching_donor1)
count_haplotype_occurrences_individuals(individuals = live_individuals, haplotype = H1)

length(mixres$pids_matching_donor2)
count_haplotype_occurrences_individuals(individuals = live_individuals, haplotype = H2)

## ------------------------------------------------------------------------
length(mixres$pids_others_included)

## ------------------------------------------------------------------------
others_haps <- get_haplotypes_pids(sim_res$population, mixres$pids_others_included)
others_haps <- others_haps[!duplicated(others_haps), ]
others_haps

## ------------------------------------------------------------------------
others_haps_counts <- unlist(lapply(seq_len(nrow(others_haps)), function(hap_i) {
            count_haplotype_occurrences_individuals(individuals = live_individuals, 
                                                    haplotype = others_haps[hap_i, ])
          }))
sum(others_haps_counts)
length(mixres$pids_others_included)

## ------------------------------------------------------------------------
rbind(H1, H2)

## ------------------------------------------------------------------------
dirichlet_alpha <- 5

set.seed(1)
sim_res <- sample_geneology(population_size = 10, 
                            generations = 10, 
                            enable_gamma_variance_extension = TRUE,
                            gamma_parameter_shape = dirichlet_alpha,
                            gamma_parameter_scale = 1 / dirichlet_alpha,
                            progress = FALSE)
pedigrees <- build_pedigrees(sim_res$population, progress = FALSE)

## ------------------------------------------------------------------------
plot(pedigrees)

## ------------------------------------------------------------------------
N <- 10000

set.seed(1)
sim_res <- sample_geneology(population_size = N, 
                            generations = 2, 
                            extra_generations_full = 2,
                            enable_gamma_variance_extension = TRUE,
                            gamma_parameter_shape = dirichlet_alpha,
                            gamma_parameter_scale = 1/dirichlet_alpha,
                            progress = FALSE, verbose_result = TRUE)

tbl_fathers_with_children <- table(sim_res$father_pids[, 1])
tbl_fathers_no_children <- rep(0, N - length(tbl_fathers_with_children))

number_of_children <- c(tbl_fathers_with_children, tbl_fathers_no_children)
number_of_children <- as.numeric(number_of_children)

mean(number_of_children)
sd(number_of_children)

## ------------------------------------------------------------------------
get_number_children <- function(N) {
  sim_res <- sample_geneology(population_size = N, 
                              generations = 2, 
                              extra_generations_full = 2,
                              enable_gamma_variance_extension = TRUE,
                              gamma_parameter_shape = dirichlet_alpha,
                              gamma_parameter_scale = 1 / dirichlet_alpha,
                              progress = FALSE, verbose_result = TRUE)
  
  tbl_fathers_with_children <- table(sim_res$father_pids[, 1])
  tbl_fathers_no_children <- rep(0, N - length(tbl_fathers_with_children))
  
  number_of_children <- c(tbl_fathers_with_children, tbl_fathers_no_children)
  number_of_children <- as.numeric(number_of_children)
  
  return(number_of_children)
}

library(parallel)
options(mc.cores = 2)
RNGkind("L'Ecuyer-CMRG") # for mclapply
set.seed(1)
x <- mclapply(1:100, function(i) get_number_children(100))
sds <- unlist(lapply(x, sd))
mean(sds)

## ------------------------------------------------------------------------
set.seed(1)
sim_res_growth <- sample_geneology_varying_size(population_sizes = c(10, 20, 10), extra_generations_full = 3, progress = FALSE)

## ------------------------------------------------------------------------
pedigrees_growth <- build_pedigrees(sim_res_growth$population, progress = FALSE)

## ------------------------------------------------------------------------
plot(pedigrees_growth)

