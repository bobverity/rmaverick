context("test-misc.R")

#------------------------------------------------
test_that("define_default working correctly", {
  
  # tests with NULL input
  expect_equal(define_default(NULL, "string"), "string")
  expect_equal(define_default(NULL, 1), 1)
  expect_equal(define_default(NULL, TRUE), TRUE)
  
  # tests with non-NULL input
  expect_equal(define_default("string", NULL), "string")
  expect_equal(define_default(1, NULL), 1)
  expect_equal(define_default(TRUE, NULL), TRUE)
})

#------------------------------------------------
test_that("C++ Hungarian algorithm working correctly", {
  
  # test with simple identity matching
  m <- 5 - diag(5,5)
  sol <- 1:nrow(m)-1
  expect_equal(call_hungarian(m)$best_matching, sol)
  
  # test with positive values
  m <- 5 - diag(5,10)
  set.seed(1)
  for (i in 1:1e2) {
    s <- sample(nrow(m))
    sol <- s-1
    expect_equal(call_hungarian(m[s,])$best_matching, sol)
  }
  
  # test with negative values
  m <- -2 - diag(5,10)
  for (i in 1:1e2) {
    s <- sample(nrow(m))
    sol <- s-1
    expect_equal(call_hungarian(m[s,])$best_matching, sol)
  }
  
  # test with more difficult problem that requires updating cost matrix through
  # both subtraction and addition
  m <- matrix(0, 3, 3)
  m[1:2, 2:3] <- 1
  sol <- 1:nrow(m)-1
  expect_equal(call_hungarian(m)$best_matching, sol)
  
})

#------------------------------------------------
test_that("log_sum working correctly", {
  
  # test finite and infinite inputs
  expect_equal(log_sum(c(3,4)), log(exp(3)+exp(4)))  # both finite
  expect_equal(log_sum(c(-Inf,4)), log(exp(-Inf)+exp(4)))  # first value infinite
  expect_equal(log_sum(c(3,-Inf)), log(exp(3)+exp(-Inf)))  # second value infinite
  expect_equal(log_sum(1:10), log(sum(exp(1:10))))  # more than two values
})
