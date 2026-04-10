test_that("MMD detects different distributions", {
  set.seed(42)
  x <- matrix(rnorm(200), 100, 2)
  y <- matrix(rnorm(200, mean = 1), 100, 2)

  res <- mmd_test(x, y, n_permutations = 200, seed = 1)
  expect_s3_class(res, "kernel_test_result")
  expect_equal(res$method, "MMD")
  expect_true(res$p_value < 0.05)
})

test_that("MMD does not reject same distribution", {
  set.seed(42)
  x <- matrix(rnorm(200), 100, 2)
  y <- matrix(rnorm(200), 100, 2)

  res <- mmd_test(x, y, n_permutations = 200, seed = 1)
  expect_true(res$p_value > 0.05)
})

test_that("MMD detects variance difference", {
  set.seed(42)
  x <- matrix(rnorm(200, sd = 1), 100, 2)
  y <- matrix(rnorm(200, sd = 3), 100, 2)

  res <- mmd_test(x, y, n_permutations = 200, seed = 1)
  expect_true(res$p_value < 0.05)
})

test_that("MMD validates inputs", {
  expect_error(mmd_test(matrix(1, 10, 2), matrix(1, 10, 3)), "same number of columns")
  expect_error(mmd_test(matrix(1, 3, 2), matrix(1, 10, 2)), "at least 5")
})
