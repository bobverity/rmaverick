context("test-misc.R")

test_that("define_default returns correct output", {
  # tests with NULL input
  expect_equal(define_default(NULL, "string"), "string")
  expect_equal(define_default(NULL, 1), 1)
  expect_equal(define_default(NULL, TRUE), TRUE)
  
  # tests with non-NULL input
  expect_equal(define_default("string", NULL), "string")
  expect_equal(define_default(1, NULL), 1)
  expect_equal(define_default(TRUE, NULL), TRUE)
})