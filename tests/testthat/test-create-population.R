context("Pedigrees and haplotypes")

test_pop <- test_create_population()

test_that("test_create_population works", {
  expect_failure(expect_null(test_pop))
  expect_output(print(test_pop), regexp = "^Population with 12 individuals$")
  expect_equal(pop_size(test_pop), 12L)
})

indvs <- get_individuals(test_pop)
test_that("get_individuals works", {
  expect_failure(expect_null(indvs))
  expect_equal(length(indvs), 12L)
})

peds <- build_pedigrees(test_pop, progress = FALSE)
test_that("build_pedigrees works", {
  expect_output(print(peds), regexp = "^List of 2 pedigrees \\(of size 11, 1\\)$")
  expect_equal(pedigrees_count(peds), 2L)
})
ped <- peds[[1L]]
pids <- sort(get_pids_in_pedigree(ped))

test_that("pedigree pids works", {
  expect_equal(length(pids), 11L)
  expect_true(all(pids == 1L:11L))
  
  expect_equal(length(get_pids_in_pedigree(peds[[2L]])), 1L)
})

test_that("meiotic_dist works", {
  expect_equal(0L, meiotic_dist(get_individual(test_pop, pid = 1L), 
                                get_individual(test_pop, pid = 1L)))
  
  expect_equal(1L, meiotic_dist(get_individual(test_pop, pid = 1L), 
                                get_individual(test_pop, pid = 6L)))
  
  expect_equal(4L, meiotic_dist(get_individual(test_pop, pid = 1L), 
                                get_individual(test_pop, pid = 10L)))
  
  expect_equal(-1L, meiotic_dist(get_individual(test_pop, pid = 1L), 
                                 get_individual(test_pop, pid = 12L)))
})



LOCI <- 5L

pedigrees_all_populate_haplotypes(peds, loci = LOCI, mutation_rates = rep(0L, LOCI), progress = FALSE)
test_that("pedigrees_all_populate_haplotypes works", {
  expect_output(print(peds), regexp = "^List of 2 pedigrees \\(of size 11, 1\\)$")
})

test_that("haplotype_matches_individuals works", {
  expect_equal(length(indvs), length(haplotype_matches_individuals(indvs, rep(0L, LOCI))))
  expect_equal(lapply(indvs, get_pid), lapply(haplotype_matches_individuals(indvs, rep(0L, LOCI)), get_pid))
})

test_that("count_haplotype_occurrences_individuals works", {
  expect_equal(12L, count_haplotype_occurrences_individuals(indvs, rep(0L, LOCI)))
  expect_equal(0L, count_haplotype_occurrences_individuals(indvs, rep(1L, LOCI)))
  
  expect_equal(11L, count_haplotype_occurrences_pedigree(ped, rep(0L, LOCI)))
  expect_equal(5L, count_haplotype_occurrences_pedigree(ped, 
                                                        rep(0L, LOCI), 
                                                        generation_upper_bound_in_result = 0L))
  expect_equal(5L+3L, count_haplotype_occurrences_pedigree(ped, 
                                                           rep(0L, LOCI), 
                                                           generation_upper_bound_in_result = 1L))
  expect_equal(5L+3L+2L, count_haplotype_occurrences_pedigree(ped, 
                                                              rep(0L, LOCI), 
                                                              generation_upper_bound_in_result = 2L))
})

mei_res <- pedigree_haplotype_matches_in_pedigree_meiosis_L1_dists(
  suspect = get_individual(test_pop, pid = 1L), 
  generation_upper_bound_in_result = -1L)
mei_res <- mei_res[order(mei_res[, 3L]), ] # order by pid
meioses <- mei_res[, 1L]
test_that("pedigree_haplotype_matches_in_pedigree_meiosis_L1_dists works", {
  expect_equal(mei_res[, 3L], 1L:11L) # pids ordered
  expect_true(all(mei_res[, 2L] == 0L)) # max L1 == 0
  
  # meioses in meioses[pid]
  expect_equal(length(meioses), 11L)
  expect_equal(meioses[1L], 0L) # no. meioses between pid = 1 and pid = 1 is 0...
  expect_equal(meioses[2L], 4L) # no. meioses between pid = 1 and pid = 2 is 4...
  expect_equal(meioses[3L], 4L)
  expect_equal(meioses[4L], 6L)
  expect_equal(meioses[5L], 6L)
  expect_equal(meioses[6L], 1L)
  expect_equal(meioses[7L], 3L)
  expect_equal(meioses[8L], 5L)
  expect_equal(meioses[9L], 2L)
  expect_equal(meioses[10L], 4L)
  expect_equal(meioses[11L], 3L)
})


haps_from_ped <- get_haplotypes_in_pedigree(ped)
haps_from_pids <- get_haplotypes_pids(test_pop, pids)
haps_from_indvs <- get_haplotypes_individuals(indvs)
hap_from_indv <- lapply(pids, function(pid) get_haplotype(get_individual(test_pop, pid)))

test_that("pedigrees_all_populate_haplotypes haplotypes works", {
  #haps_from_ped
  expect_true(is.list(haps_from_ped))
  expect_equal(length(haps_from_ped), 11L)
  expect_equal(length(haps_from_ped[[1L]]), LOCI)
  expect_true(all(unlist(haps_from_ped) == 0L))
  
  #haps_from_pids
  expect_true(is.matrix(haps_from_pids))
  expect_equal(nrow(haps_from_pids), 11L)
  expect_equal(ncol(haps_from_pids), LOCI)
  
  #haps_from_indvs
  expect_equal(nrow(haps_from_indvs), 12L)
  expect_equal(unique(c(haps_from_indvs)), 0L)
  
  #hap_from_indv
  expect_equal(haps_from_ped, hap_from_indv)
})




test_that("pedigrees_all_populate_haplotypes_ladder_bounded error works", {
  #haps_from_ped
  expect_error(
    pedigrees_all_populate_haplotypes_ladder_bounded(peds, 
                                                     mutation_rates = rep(1, LOCI), 
                                                     ladder_min = rep(10L, LOCI), 
                                                     ladder_max = rep(10L, LOCI), 
                                                     get_founder_haplotype = function() rep(10L, LOCI),
                                                     progress = FALSE)
  )
})


pedigrees_all_populate_haplotypes_ladder_bounded(peds, 
                                                 mutation_rates = rep(0, LOCI), 
                                                 ladder_min = rep(10L, LOCI), 
                                                 ladder_max = rep(11L, LOCI), 
                                                 get_founder_haplotype = function() rep(10L, LOCI),
                                                 progress = FALSE)

haps_from_ped <- get_haplotypes_in_pedigree(ped)
haps_from_pids <- get_haplotypes_pids(test_pop, pids)
haps_from_indvs <- get_haplotypes_individuals(indvs)
hap_from_indv <- lapply(pids, function(pid) get_haplotype(get_individual(test_pop, pid)))

test_that("pedigrees_all_populate_haplotypes_ladder_bounded haplotypes works", {
  #haps_from_ped
  expect_true(is.list(haps_from_ped))
  expect_equal(length(haps_from_ped), 11L)
  expect_equal(length(haps_from_ped[[1L]]), LOCI)
  expect_true(all(unlist(haps_from_ped) == 10L))
  
  #haps_from_pids
  expect_true(is.matrix(haps_from_pids))
  expect_equal(nrow(haps_from_pids), 11L)
  expect_equal(ncol(haps_from_pids), LOCI)
  
  #haps_from_indvs
  expect_equal(nrow(haps_from_indvs), 12L)
  expect_equal(unique(c(haps_from_indvs)), 10L)
  
  #hap_from_indv
  expect_equal(haps_from_ped, hap_from_indv)
})



#####################################

set.seed(1)
pedigrees_all_populate_haplotypes_ladder_bounded(peds, 
                                                 mutation_rates = rep(1, LOCI), 
                                                 ladder_min = rep(10L, LOCI), 
                                                 ladder_max = rep(20L, LOCI), 
                                                 get_founder_haplotype = function() rep(10L, LOCI),
                                                 progress = FALSE)

haps_from_pids <- get_haplotypes_pids(test_pop, pids)
hap_fac <- factor(apply(haps_from_pids, 1, paste0, collapse = ";"))
hap_hash <- as.integer(hap_fac)
hashes <- haplotypes_to_hashes(test_pop, pids)

test_that("hashes works", {
  expect_equal(length(hap_fac), length(hashes))
  for (j in seq_along(hap_hash)) {
    #j <- 1
    #j <- 2
    is <- which(hap_hash == hap_hash[j])
    ks <- which(hashes == hashes[j])
    expect_equal(is, ks)
  }
})

pids_split_hap <- lapply(split_by_haplotypes(test_pop, pids), sort)
x <- lapply(split(pids, hap_fac), sort)
names(x) <- NULL

# Order list in some way that make them comparable:
#  Here, first by size and then by their max element; ensures unique ordering 
#  because all elements are pids and are hence unique integers.
ord_1 <- order(sapply(pids_split_hap, length), sapply(pids_split_hap, max))
ord_2 <- order(sapply(x, length), sapply(x, max))
pids_split_hap <- pids_split_hap[ord_1]
x <- x[ord_2]

test_that("split pid by haplotype works", {
  expect_equal(length(hap_fac), length(hashes))
})
