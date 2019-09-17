context("Match probability")

test_that("Test simple", {
  set.seed(1)
  sim_res_growth <- sample_geneology_varying_size(population_sizes = rep(1e3, 200),
                                                  generations_full = 3,
                                                  generations_return = 3,
                                                  enable_gamma_variance_extension = TRUE,
                                                  gamma_parameter_shape = 5,
                                                  gamma_parameter_scale = 1/5,
                                                  progress = FALSE)
  
  peds <- build_pedigrees(sim_res_growth$population, progress = FALSE)
  
  live_individuals <- sim_res_growth$individuals_generations
  indv_poi <- live_individuals[[sample.int(n = length(live_individuals), size = 1)]]
  
  mdist <- meiotic_distribution(poi = indv_poi, indvs = live_individuals)
  
  expect_equal(mdist, structure(list(m = c(1L, 4L, 5L, 6L, 7L, 8L, 16L, 18L, 30L, 31L, 
                                           32L, 48L, 49L, 50L, 52L, 53L, 54L, 55L, 56L, 70L, 71L, 72L, 86L, 
                                           87L, 88L, 204L, 205L, 206L, 338L, 339L, 340L, 352L, 353L, 354L
  ), n = c(1L, 3L, 4L, 9L, 3L, 2L, 2L, 1L, 6L, 8L, 15L, 12L, 14L, 
           8L, 47L, 60L, 144L, 69L, 77L, 41L, 35L, 35L, 20L, 20L, 21L, 36L, 
           29L, 37L, 111L, 103L, 112L, 35L, 42L, 44L), p = c(0.000829187396351575, 
                                                             0.00248756218905473, 0.0033167495854063, 0.00746268656716418, 
                                                             0.00248756218905473, 0.00165837479270315, 0.00165837479270315, 
                                                             0.000829187396351575, 0.00497512437810945, 0.0066334991708126, 
                                                             0.0124378109452736, 0.00995024875621891, 0.0116086235489221, 
                                                             0.0066334991708126, 0.038971807628524, 0.0497512437810945, 0.119402985074627, 
                                                             0.0572139303482587, 0.0638474295190713, 0.0339966832504146, 0.0290215588723051, 
                                                             0.0290215588723051, 0.0165837479270315, 0.0165837479270315, 0.0174129353233831, 
                                                             0.0298507462686567, 0.0240464344941957, 0.0306799336650083, 0.0920398009950249, 
                                                             0.0854063018242123, 0.0928689883913764, 0.0290215588723051, 0.0348258706467662, 
                                                             0.0364842454394693), cum_p = c(0.000829187396351575, 0.0033167495854063, 
                                                                                            0.0066334991708126, 0.0140961857379768, 0.0165837479270315, 0.0182421227197347, 
                                                                                            0.0199004975124378, 0.0207296849087894, 0.0257048092868988, 0.0323383084577114, 
                                                                                            0.0447761194029851, 0.054726368159204, 0.066334991708126, 0.0729684908789386, 
                                                                                            0.111940298507463, 0.161691542288557, 0.281094527363184, 0.338308457711443, 
                                                                                            0.402155887230514, 0.436152570480929, 0.465174129353234, 0.494195688225539, 
                                                                                            0.51077943615257, 0.527363184079602, 0.544776119402985, 0.574626865671642, 
                                                                                            0.598673300165838, 0.629353233830846, 0.721393034825871, 0.806799336650083, 
                                                                                            0.899668325041459, 0.928689883913765, 0.963515754560531, 1)), row.names = c(NA, 
                                                                                                                                                                        -34L), class = c("malan_meiotic_dist", "data.frame")))
  
  mp <- matchprob(meiotic_dist = mdist, 
                  mut_rate = 1 - prod(1 - rep(0.001, 20)))
  
  expect_equal(mp, 0.18583573467031)
})



test_that("Test multiple", {
  set.seed(1)
  
  n_pops <- 10
  n_pois <- 5
  
  res <- lapply(seq_len(n_pops), function(i) {
    sim_res_growth <- sample_geneology_varying_size(population_sizes = rep(1e3, 200),
                                                    generations_full = 3,
                                                    generations_return = 3,
                                                    enable_gamma_variance_extension = TRUE,
                                                    gamma_parameter_shape = 5,
                                                    gamma_parameter_scale = 1/5,
                                                    progress = FALSE)
    
    peds <- build_pedigrees(sim_res_growth$population, progress = FALSE)
    
    live_individuals <- sim_res_growth$individuals_generations
    
    res_pois <- sapply(seq_len(n_pois), function(j) {
      indv_poi <- live_individuals[[sample.int(n = length(live_individuals), size = 1)]]
      mdist <- meiotic_distribution(poi = indv_poi, indvs = live_individuals)
      mp <- matchprob(meiotic_dist = mdist, 
                      mut_rate = 1 - prod(1 - rep(0.001, 20)))
      return(mp)
    })
  })
  
  expect_equal(length(res), n_pops)
  expect_equal(unique(sapply(res, length)), n_pois)
  
  res_vec <- unlist(res)
  expect_equal(mean(res_vec), 0.195726438053142)
})



###############################################
# 

test_that("meiotis_dist_all", {
  set.seed(1)
  sim_res_growth <- sample_geneology_varying_size(population_sizes = rep(100, 100),
                                                  generations_return = 100,
                                                  progress = FALSE)
  
  peds <- build_pedigrees(sim_res_growth$population, progress = FALSE)
  #pedigrees_count(peds)
  
  all_individuals <- sim_res_growth$individuals_generations
  indv_poi <- all_individuals[[sample.int(n = length(all_individuals), size = 1)]]
  
  indv_poi_dists_all <- meiotis_dist_all(indv_poi)
  mdist_all <- as.data.frame(table(as.numeric(indv_poi_dists_all)))
  colnames(mdist_all) <- c("m", "n")
  mdist_all$m <- as.integer(as.character(mdist_all$m))
  expect_equal(1, as.numeric(subset(mdist_all, m == 0, n)))
  
  mdist_all <- subset(mdist_all, m > 0L)
  rownames(mdist_all) <- NULL
  mdist <- meiotic_distribution(poi = indv_poi, indvs = all_individuals)[, c("m", "n")]
  rownames(mdist) <- NULL
  class(mdist) <- "data.frame"
  expect_equal(mdist_all, mdist)
  
  if (FALSE) {
    microbenchmark::microbenchmark(
      meiotis_dist_all(indv_poi),
      meiotic_distribution(poi = indv_poi, indvs = all_individuals),
      times = 10
    )
  }
  
  if (FALSE) {
    # meiotis_dist_all_lookup
    # FIXME
    # TODO
  }
})



