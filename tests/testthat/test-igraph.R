context("igraph")

test_that("igraph import works for tree", {
  g1 <- graph_from_literal( 2 +- 1 -+ 3 )
  pop1 <- from_igraph(g1)
  peds1 <- build_pedigrees(pop1, progress = FALSE)
  
  expect_error(get_generation(get_individual(pop1, 1)))
  expect_error(get_generation(get_individual(pop1, 2)))
  expect_error(get_generation(get_individual(pop1, 3)))

  infer_generations(peds1)
  
  expect_equal(get_generation(get_individual(pop1, 1)), 1L)
  expect_equal(get_generation(get_individual(pop1, 2)), 0L)
  expect_equal(get_generation(get_individual(pop1, 3)), 0L)
})


test_that("igraph import works for forest", {
  g2 <- graph_from_literal( 2 +- 1 -+ 3, 4 -+ 5 )
  pop2 <- from_igraph(g2)
  peds2 <- build_pedigrees(pop2, progress = FALSE)
  
  expect_error(get_generation(get_individual(pop2, 1)))
  expect_error(get_generation(get_individual(pop2, 2)))
  expect_error(get_generation(get_individual(pop2, 3)))
  expect_error(get_generation(get_individual(pop2, 4)))
  expect_error(get_generation(get_individual(pop2, 5)))
  
  infer_generations(peds2)
  
  expect_equal(get_generation(get_individual(pop2, 1)), 1L)
  expect_equal(get_generation(get_individual(pop2, 2)), 0L)
  expect_equal(get_generation(get_individual(pop2, 3)), 0L)
  expect_equal(get_generation(get_individual(pop2, 4)), 1L)
  expect_equal(get_generation(get_individual(pop2, 5)), 0L)
})
