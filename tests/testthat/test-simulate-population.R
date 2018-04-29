context("Simulate population")

set.seed(1)
sim_res_fixed <- sample_geneology(population_size = 1e3, 
                                  generations = 20, 
                                  generations_full = 3,
                                  generations_return = 3,
                                  progress = FALSE)

test_that("sample_geneology works", {
  expect_failure(expect_null(sim_res_fixed))
  expect_output(print(sim_res_fixed$population), regexp = "^Population with .* individuals$")
  expect_equal(length(sim_res_fixed$end_generation_individuals), 1000L)
  expect_equal(length(sim_res_fixed$individuals_generations), 3L*1000L)
})


set.seed(2)
sim_res_growth <- sample_geneology_varying_size(population_sizes = rep(1e3, 20),
                                                enable_gamma_variance_extension = TRUE,
                                                gamma_parameter_shape = 5,
                                                gamma_parameter_scale = 1/5,
                                                generations_full = 3,
                                                generations_return = 3, # default value
                                                progress = FALSE)

test_that("sample_geneology works", {
  expect_failure(expect_null(sim_res_growth))
  expect_output(print(sim_res_growth$population), regexp = "^Population with .* individuals$")
  expect_equal(length(sim_res_growth$end_generation_individuals), 1000L)
  expect_equal(length(sim_res_growth$individuals_generations), 3L*1000L)
  expect_equal(sim_res_growth$sdo_type, "GammaVariation")
})

