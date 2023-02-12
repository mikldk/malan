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
  if (FALSE) {
    plot(peds)
  }
  expect_equal(4L, meiotic_dist(get_individual(test_pop, pid = 1L), 
                                get_individual(test_pop, pid = 10L)))
  
  expect_equal(6L, meiotic_dist(get_individual(test_pop, pid = 1L), 
                                get_individual(test_pop, pid = 5L)))
  
  expect_equal(-1L, meiotic_dist(get_individual(test_pop, pid = 37L), 
                                 get_individual(test_pop, pid = 10L)))
})

test_that("meiotic_dist_threshold works", {
  expect_equal(4L, meiotic_dist_threshold(get_individual(test_pop, pid = 1L), 
                                          get_individual(test_pop, pid = 10L),
                                          threshold = 5))
  expect_equal(4L, meiotic_dist_threshold(get_individual(test_pop, pid = 1L), 
                                          get_individual(test_pop, pid = 10L),
                                          threshold = 4))
  expect_equal(-1L, meiotic_dist_threshold(get_individual(test_pop, pid = 1L), 
                                          get_individual(test_pop, pid = 10L),
                                          threshold = 3))
  expect_equal(-1L, meiotic_dist_threshold(get_individual(test_pop, pid = 1L), 
                                           get_individual(test_pop, pid = 10L),
                                           threshold = 2))
  
  
  expect_equal(6L, meiotic_dist_threshold(get_individual(test_pop, pid = 1L), 
                                          get_individual(test_pop, pid = 5L),
                                          threshold = 7))
  expect_equal(6L, meiotic_dist_threshold(get_individual(test_pop, pid = 1L), 
                                          get_individual(test_pop, pid = 5L),
                                          threshold = 6))
  expect_equal(-1L, meiotic_dist_threshold(get_individual(test_pop, pid = 1L), 
                                          get_individual(test_pop, pid = 5L),
                                          threshold = 5))
  
  expect_equal(-1L, meiotic_dist_threshold(get_individual(test_pop, pid = 37L), 
                                           get_individual(test_pop, pid = 10L),
                                           threshold = 1))
  expect_equal(-1L, meiotic_dist_threshold(get_individual(test_pop, pid = 37L), 
                                           get_individual(test_pop, pid = 10L),
                                           threshold = 0))
})

test_that("radius works 1", {
  if (FALSE) {
    plot(peds)
  }
  
  ans <- as.data.frame(meiotic_radius(get_individual(test_pop, pid = 1L), radius = 1))
  expect_equal(nrow(ans), 2)
  expect_equal(subset(ans, pid == 1, dist, drop = TRUE), 0L)
  expect_equal(subset(ans, pid == 6, dist, drop = TRUE), 1L)
  
  ans <- as.data.frame(meiotic_radius(get_individual(test_pop, pid = 1L), radius = 2))
  expect_equal(nrow(ans), 3)
  expect_equal(subset(ans, pid == 1, dist, drop = TRUE), 0L)
  expect_equal(subset(ans, pid == 6, dist, drop = TRUE), 1L)
  expect_equal(subset(ans, pid == 9, dist, drop = TRUE), 2L)
  
  ans <- as.data.frame(meiotic_radius(get_individual(test_pop, pid = 1L), radius = 3))
  expect_equal(nrow(ans), 5)
  expect_equal(subset(ans, pid == 1, dist, drop = TRUE), 0L)
  expect_equal(subset(ans, pid == 6, dist, drop = TRUE), 1L)
  expect_equal(subset(ans, pid == 9, dist, drop = TRUE), 2L)
  expect_equal(subset(ans, pid == 11, dist, drop = TRUE), 3L)
  expect_equal(subset(ans, pid == 7, dist, drop = TRUE), 3L)
  
  ans <- as.data.frame(meiotic_radius(get_individual(test_pop, pid = 1L), radius = 4))
  expect_equal(nrow(ans), 8)
  expect_equal(subset(ans, pid == 1, dist, drop = TRUE), 0L)
  expect_equal(subset(ans, pid == 6, dist, drop = TRUE), 1L)
  expect_equal(subset(ans, pid == 9, dist, drop = TRUE), 2L)
  expect_equal(subset(ans, pid == 11, dist, drop = TRUE), 3L)
  expect_equal(subset(ans, pid == 7, dist, drop = TRUE), 3L)
  expect_equal(subset(ans, pid == 2, dist, drop = TRUE), 4L)
  expect_equal(subset(ans, pid == 3, dist, drop = TRUE), 4L)
  expect_equal(subset(ans, pid == 10, dist, drop = TRUE), 4L)
  
  ans <- as.data.frame(meiotic_radius(get_individual(test_pop, pid = 1L), radius = 5))
  expect_equal(nrow(ans), 9)
  expect_equal(subset(ans, pid == 1, dist, drop = TRUE), 0L)
  expect_equal(subset(ans, pid == 6, dist, drop = TRUE), 1L)
  expect_equal(subset(ans, pid == 9, dist, drop = TRUE), 2L)
  expect_equal(subset(ans, pid == 11, dist, drop = TRUE), 3L)
  expect_equal(subset(ans, pid == 7, dist, drop = TRUE), 3L)
  expect_equal(subset(ans, pid == 2, dist, drop = TRUE), 4L)
  expect_equal(subset(ans, pid == 3, dist, drop = TRUE), 4L)
  expect_equal(subset(ans, pid == 10, dist, drop = TRUE), 4L)
  expect_equal(subset(ans, pid == 8, dist, drop = TRUE), 5L)
  
  #####
  
  ans <- as.data.frame(meiotic_radius(get_individual(test_pop, pid = 37L), radius = 1))
  expect_equal(nrow(ans), 5)
  expect_equal(subset(ans, pid == 39, dist, drop = TRUE), 1L)
  expect_equal(subset(ans, pid == 33, dist, drop = TRUE), 1L)
  expect_equal(subset(ans, pid == 34, dist, drop = TRUE), 1L)
  expect_equal(subset(ans, pid == 35, dist, drop = TRUE), 1L)
  
  ans <- as.data.frame(meiotic_radius(get_individual(test_pop, pid = 37L), radius = 2))
  expect_equal(nrow(ans), 7)
  expect_equal(subset(ans, pid == 39, dist, drop = TRUE), 1L)
  expect_equal(subset(ans, pid == 33, dist, drop = TRUE), 1L)
  expect_equal(subset(ans, pid == 34, dist, drop = TRUE), 1L)
  expect_equal(subset(ans, pid == 35, dist, drop = TRUE), 1L)
  expect_equal(subset(ans, pid == 36, dist, drop = TRUE), 2L)
  expect_equal(subset(ans, pid == 38, dist, drop = TRUE), 2L)
  
  ans <- as.data.frame(meiotic_radius(get_individual(test_pop, pid = 37L), radius = 3))
  expect_equal(nrow(ans), 9)
  expect_equal(subset(ans, pid == 39, dist, drop = TRUE), 1L)
  expect_equal(subset(ans, pid == 33, dist, drop = TRUE), 1L)
  expect_equal(subset(ans, pid == 34, dist, drop = TRUE), 1L)
  expect_equal(subset(ans, pid == 35, dist, drop = TRUE), 1L)
  expect_equal(subset(ans, pid == 36, dist, drop = TRUE), 2L)
  expect_equal(subset(ans, pid == 38, dist, drop = TRUE), 2L)
  expect_equal(subset(ans, pid == 31, dist, drop = TRUE), 3L)
  expect_equal(subset(ans, pid == 32, dist, drop = TRUE), 3L)
})

test_that("radius works 2", {
  if (FALSE) {
    plot(peds)
    
    ans <- meiotic_radius(get_individual(test_pop, pid = 1L), radius = 1)
    ans

    ans <- meiotic_radius(get_individual(test_pop, pid = 1L), radius = 2)
    ans
    
    ans <- meiotic_radius(get_individual(test_pop, pid = 1L), radius = 3)
    ans
    
    ans <- meiotic_radius(get_individual(test_pop, pid = 1L), radius = 4)
    ans
  }
  
  verify_distances <- function(pid_PoI, radius) {
    ans <- meiotic_radius(get_individual(test_pop, pid = pid_PoI), radius = radius)
    ans
    expect_true(isTRUE(all(ans[, 2] <= radius)))
    
    for (i in seq_len(nrow(ans))) {
      expect_equal(as.integer(ans[i, 2L]), 
                   meiotic_dist(get_individual(test_pop, pid = pid_PoI), 
                                get_individual(test_pop, pid = as.integer(ans[i, 1L]))), 
                   label = paste0("pid_PoI = ", pid_PoI, " / radius = ", radius, " / pid ", as.integer(ans[i, 1L])))
    }
  }
  
  for (radius in seq_len(10)) {
    verify_distances(pid_PoI = 1L, radius = radius)
  }
  
  for (radius in seq_len(10)) {
    verify_distances(pid_PoI = 34L, radius = radius)
  }
})

test_that("radius works 3", {
  set.seed(20230212)
  radius_pop <- sample_geneology(population_size = 1e2, 
                                 generations = 30,
                                 generations_full = 3,
                                 generations_return = 3,
                                 progress = FALSE)
  radius_peds <- build_pedigrees(radius_pop$population, progress = FALSE)
  
  if (FALSE) {
    plot(radius_peds)
  }
  
  radius_verify_distances <- function(pid_PoI, radius) {
    #cat("pid_PoI = ", pid_PoI, " / radius = ", radius, "\n", sep = "")
           
    ans <- meiotic_radius(get_individual(radius_pop$population, pid = pid_PoI), radius = radius)
    ans
    expect_true(isTRUE(all(ans[, 2] <= radius)))
    
    for (i in seq_len(nrow(ans))) {
      expect_equal(as.integer(ans[i, 2L]), 
                   meiotic_dist(get_individual(radius_pop$population, pid = pid_PoI), 
                                get_individual(radius_pop$population, pid = as.integer(ans[i, 1L]))), 
                   label = paste0("pid_PoI = ", pid_PoI, " / radius = ", radius, " / pid ", as.integer(ans[i, 1L])))
    }
  }
  
  radius_pop_live_pids <- lapply(radius_pop$individuals_generations, get_pid)
  pid_PoIs <- sample(x = radius_pop_live_pids, size = 5)
  
  for (pid_PoI in pid_PoIs) {
    for (radius in seq_len(15)) {
      radius_verify_distances(pid_PoI = pid_PoI, radius = radius)
    }
  }
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


