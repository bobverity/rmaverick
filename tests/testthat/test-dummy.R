context("test-dummy.R")

test_that("dummy1 returns correct output", {
  expect_equal(dummy1(), "dummy1")
  expect_equal(dummy1(5), "dummy5")
})