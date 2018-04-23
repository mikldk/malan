context("Family")

test_pop <- test_create_population()
indvs <- get_individuals(test_pop)
peds <- build_pedigrees(test_pop, progress = FALSE)
ped <- peds[[1L]]

if (FALSE) {
  plot(ped)
}

i1 <- get_individual(test_pop, 1L)
i2 <- get_individual(test_pop, 2L)
i4 <- get_individual(test_pop, 4L)
i6 <- get_individual(test_pop, 6L)
i7 <- get_individual(test_pop, 7L)
i8 <- get_individual(test_pop, 8L)

test_that("children works", {
  expect_equal(0L, length(get_children(i1)))
  expect_equal(0L, length(get_children(i4)))
  expect_equal(1L, length(get_children(i6)))
  expect_equal(2L, length(get_children(i7)))
  
  expect_equal(sort(unlist(lapply(get_children(i7), get_pid))), c(2L, 3L))
  
  expect_equal(sort(unlist(lapply(get_children(i6), get_pid))), 1L)
  expect_s3_class(get_children(i6)[[1]], "malan_individual")
  expect_s3_class(get_children(i6)[[1]], "externalptr")
})

test_that("brothers works", {
  # pid 1
  expect_equal(0L, count_brothers(i1)) 
  expect_true(is.list(get_brothers(i1)))
  expect_equal(length(get_brothers(i1)), count_brothers(i1)) 
  
  # pid 4
  expect_equal(1L, count_brothers(i4)) # 5
  expect_equal(length(get_brothers(i4)), count_brothers(i4)) # 7
  expect_equal(unlist(lapply(get_brothers(i4), get_pid)), 5L)
  expect_s3_class(get_brothers(i4)[[1]], "malan_individual")
  expect_s3_class(get_brothers(i4)[[1]], "externalptr")
})

test_that("uncles works", {
  # pid 1
  expect_equal(1L, count_uncles(i1)) # 7
  expect_true(is.list(get_uncles(i1)))
  expect_equal(length(get_uncles(i1)), count_uncles(i1)) # 7
  expect_equal(unlist(lapply(get_uncles(i1), get_pid)), 7L)
  expect_s3_class(get_uncles(i1)[[1]], "malan_individual")
  expect_s3_class(get_uncles(i1)[[1]], "externalptr")
  
  # pid 4
  expect_equal(0L, count_uncles(i4)) # 7
  expect_equal(length(get_uncles(i4)), count_uncles(i4)) # 7
})

test_that("cousins works", {
  expect_equal(0L, length(get_cousins(i4)))
  expect_equal(1L, length(get_cousins(i6)))
  
  expect_equal(2L, length(get_cousins(i1)))
  expect_equal(sort(unlist(lapply(get_cousins(i1), get_pid))), c(2L, 3L))
  
  expect_equal(1L, length(get_cousins(i2)))
  expect_equal(sort(unlist(lapply(get_cousins(i2), get_pid))), 1L)
  
  expect_equal(2L, length(get_cousins(i8)))
  expect_equal(sort(unlist(lapply(get_cousins(i8), get_pid))), c(6L, 7L))
  
  expect_equal(sort(unlist(lapply(get_cousins(i6), get_pid))), 8L)
  expect_s3_class(get_cousins(i6)[[1]], "malan_individual")
  expect_s3_class(get_cousins(i6)[[1]], "externalptr")
})
