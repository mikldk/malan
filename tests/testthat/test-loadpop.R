context("Loading data")

# Pedigree 1:
#     
# G3           11      
#           /     \
# G2      9        10
#        /  \       |
# G1    6    7      8
#       |   /\     /\
# G0    1   2 3   4  5
#
#
# Pedigree 2:
#     
# G3       39  _______   
#         /   \       \  
# G1    36     37       38
#       |\     / | \
# G0   31 32  33 34 35
loaded_data <- structure(list(pid = c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 10L, 11L,
                                      31L, 32L, 33L, 34L, 35L, 36L, 37L, 38L, 39L), 
                              pid_dad = c(6L, 7L, 7L, 8L, 8L, 9L, 9L, 10L, 11L, 11L, 0L, 
                                          36L, 36L, 37L, 37L, 37L, 39L, 39L, 39L, 0L)), 
                         class = c("tbl_df", "tbl", "data.frame"), row.names = c(NA, -20L))

test_pop <- load_individuals(pid = loaded_data$pid, 
                             pid_dad = loaded_data$pid_dad,
                             progress = FALSE)

test_that("test_create_population works", {
  expect_failure(expect_null(test_pop))
  expect_output(print(test_pop), regexp = "^Population with 20 individuals$")
  expect_equal(pop_size(test_pop), 20L)
})

indvs <- get_individuals(test_pop)
test_that("get_individuals works", {
  expect_failure(expect_null(indvs))
  expect_equal(length(indvs), 20L)
})

test_that("generation assignment works", {
  
  # 0, last generation:
  expect_equal(0L, 
               unique(unlist(lapply(1L:5L, function(pid) get_generation(get_individual(test_pop, pid)))))
  )
  expect_equal(0L, 
               unique(unlist(lapply(31L:35L, function(pid) get_generation(get_individual(test_pop, pid)))))
  )
  
  # 1, second last generation:
  expect_equal(1L, 
               unique(unlist(lapply(6L:8L, function(pid) get_generation(get_individual(test_pop, pid)))))
  )
  expect_equal(1L, 
               unique(unlist(lapply(36L:38L, function(pid) get_generation(get_individual(test_pop, pid)))))
  )
  
  # 2, third last generation:
  expect_equal(2L, 
               unique(unlist(lapply(9L:10L, function(pid) get_generation(get_individual(test_pop, pid)))))
  )
  expect_equal(2L, 
               unique(unlist(lapply(39L, function(pid) get_generation(get_individual(test_pop, pid)))))
  )
  
  # 3, fourth last generation:
  expect_equal(3L, 
               unique(unlist(lapply(11L, function(pid) get_generation(get_individual(test_pop, pid)))))
  )
})

peds <- build_pedigrees(test_pop, progress = FALSE)
test_that("build_pedigrees works", {
  expect_output(print(peds), regexp = "^List of 2 pedigrees \\(of size 11, 9\\)$")
  expect_equal(pedigrees_count(peds), 2L)
})
ped <- peds[[1L]]

pids <- sort(get_pids_in_pedigree(ped))
test_that("pedigree pids works", {
  expect_equal(length(pids), 11L)
  expect_true(all(pids == 1L:11L))
})

# ped 1:
pids <- sort(get_pids_in_pedigree(peds[[1L]]))
test_that("pedigree pids works", {
  expect_equal(length(pids), 11L)
  expect_true(all(pids == 1L:11L))
})

# ped 2:
pids <- sort(get_pids_in_pedigree(peds[[2L]]))
test_that("pedigree pids works", {
  expect_equal(length(pids), 9L)
  expect_true(all(pids == 31L:39L))
})


test_that("meiotic_dist works", {
  expect_equal(4L, meiotic_dist(get_individual(test_pop, pid = 1L), 
                                get_individual(test_pop, pid = 10L)))
})

LOCI <- 5L

pedigrees_all_populate_haplotypes(peds, loci = LOCI, mutation_rates = rep(0, LOCI), progress = FALSE)
test_that("haplotype_matches_individuals works", {
  expect_equal(length(indvs), length(haplotype_matches_individuals(indvs, rep(0L, LOCI))))
  expect_equal(lapply(indvs, get_pid), lapply(haplotype_matches_individuals(indvs, rep(0L, LOCI)), get_pid))
})

################################################



indvs <- get_individuals(test_pop)

set.seed(1)
pedigrees_all_populate_haplotypes(peds, loci = LOCI, mutation_rates = rep(0.1, LOCI), progress = FALSE)

test_that("count_haplotype_occurrences_individuals works", {
  expect_true(0L < count_haplotype_occurrences_individuals(indvs, rep(0L, LOCI)))
})

test_that("get_matching_pids_from_hashmap works for (0, 0, ..., 0)", {
  hashmap <- build_haplotype_hashmap(indvs, progress = FALSE)
  pids <- get_matching_pids_from_hashmap(hashmap, rep(0L, LOCI))
  
  expect_equal(length(pids), count_haplotype_occurrences_individuals(indvs, rep(0L, LOCI)))
  
  delete_haplotypeids_hashmap(hashmap)
})


test_that("get_matching_pids_from_hashmap works for all", {
  set.seed(1)
  pedigrees_all_populate_haplotypes(peds, loci = LOCI, mutation_rates = rep(0.5, LOCI), progress = FALSE)
  all_hs <- get_haplotypes_individuals(indvs)
  all_hs <- all_hs[!duplicated(all_hs), ]
  
  hashmap <- build_haplotype_hashmap(indvs, progress = FALSE)
  
  for (i in seq_len(nrow(all_hs))) {
    h <- all_hs[i, ]
    
    pids <- get_matching_pids_from_hashmap(hashmap, h)
    
    expect_equal(length(pids), count_haplotype_occurrences_individuals(indvs, h), info = paste0("i = ", i, "; h = (", paste0(h, collapse = ", "), ")"))
  }

  delete_haplotypeids_hashmap(hashmap)
})


