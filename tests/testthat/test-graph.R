context("Graph")

test_pop <- test_create_population()
indvs <- get_individuals(test_pop)
peds <- build_pedigrees(test_pop, progress = FALSE)
ped <- peds[[1L]]

if (FALSE) {
  plot(ped)
}

pid1_meiotic_dist_calc <- meioses_generation_distribution(get_individual(test_pop, pid = 1L))
lines <- "generation meioses count
          0          0       1
          0          4       2
          0          6       2
          1          1       1
          1          3       1
          1          5       1
          2          2       1
          2          4       1
          3          3       1"
lines <- gsub("[ ]+", " ", lines)
lines <- gsub("\n ", "\n", lines)
con <- textConnection(lines)
pid1_meiotic_dist_known <- as.matrix(read.csv(con, sep = " "))

test_that("meioses_generation_distribution works", {
  expect_equal(pid1_meiotic_dist_calc, pid1_meiotic_dist_known)
})


library(igraph)
g <- pedigree_as_igraph(ped)
el <- igraph::as_edgelist(g)
el2 <- as.matrix(get_nodes_edges(peds)$edges)
dimnames(el2) <- NULL

test_that("igraph interface works", {
  expect_equal(sort(as.integer(V(g))), 1L:11L)
  expect_equal(el, el2)
  
  expect_equal(0L, c(igraph::distances(g, v = "1", to = "1")))
  
  expect_equal(1L, c(igraph::distances(g, v = "1", to = "6")))
  
  expect_equal(4L, c(igraph::distances(g, v = "1", to = "10")))
  
  for (v1 in 1L:11L) {
    for (v2 in v1:11L) {
      expect_equal(c(igraph::distances(g, v = as.character(v1), to = as.character(v2))), 
                   meiotic_dist(get_individual(test_pop, pid = v1), 
                                get_individual(test_pop, pid = v2)))
    }
  }
})
