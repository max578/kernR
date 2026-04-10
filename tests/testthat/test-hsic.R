test_that("HSIC detects dependence", {
  set.seed(42)
  n <- 200
  x <- rnorm(n)
  y <- x + rnorm(n, sd = 0.5)

  res <- hsic_test(x, y, n_permutations = 200, seed = 1)
  expect_s3_class(res, "kernel_test_result")
  expect_equal(res$method, "HSIC")
  expect_true(res$p_value < 0.05)
})

test_that("HSIC does not reject independence", {
  set.seed(42)
  n <- 200
  x <- rnorm(n)
  y <- rnorm(n)

  res <- hsic_test(x, y, n_permutations = 200, seed = 1)
  expect_true(res$p_value > 0.05)
})

test_that("HSIC detects non-linear dependence", {
  set.seed(42)
  n <- 300
  x <- rnorm(n)
  y <- x^2 + rnorm(n, sd = 0.3)

  res <- hsic_test(x, y, n_permutations = 200, seed = 1)
  expect_true(res$p_value < 0.05)
})

test_that("HSIC validates inputs", {
  expect_error(hsic_test(1:5, 1:10), "same number")
  expect_error(hsic_test(1:5, 1:5), "At least 10")
})

test_that("HSIC result prints without error", {
  set.seed(1)
  res <- hsic_test(rnorm(50), rnorm(50), n_permutations = 50, seed = 1)
  expect_output(print(res), "HSIC")
  expect_output(summary(res), "HSIC")
})
