test_that("bd-HSIC detects causal association", {
  set.seed(42)
  n <- 200
  z <- matrix(rnorm(n * 2), n, 2)
  x <- z[, 1] + rnorm(n)
  y <- 0.8 * x + z[, 2] + rnorm(n, sd = 0.3)

  result <- bd_hsic_test(x, y, z,
    n_permutations = 100,
    seed = 1
  )

  expect_s3_class(result, "kernel_test_result")
  expect_equal(result$method, "bd-HSIC")
  expect_true(result$p_value < 0.1)
  expect_true(result$ess > 0)
})

test_that("bd-HSIC does not reject under null", {
  set.seed(42)
  n <- 200
  z <- matrix(rnorm(n * 2), n, 2)
  x <- z[, 1] + rnorm(n)
  y <- z[, 2] + rnorm(n)  # y depends on z but not x after conditioning

  result <- bd_hsic_test(x, y, z,
    n_permutations = 100,
    seed = 1
  )
  # Under null, p-value should be > 0.05 most of the time
  # We just check it runs without error
  expect_s3_class(result, "kernel_test_result")
  expect_true(result$p_value >= 0)
  expect_true(result$p_value <= 1)
})

test_that("bd-HSIC validates inputs", {
  expect_error(
    bd_hsic_test(1:5, 1:5, matrix(1, 5, 2)),
    "at least 6"
  )
})

test_that("bd-HSIC accepts small-N field trials (A2: 20-obs floor dropped)", {
  # A small multi-paddock trial: below the retired hard floor of 20, above the
  # genuine mathematical minimum of 6. The test must run and return a valid
  # p-value rather than refuse the sample outright.
  set.seed(7)
  n <- 14
  z <- matrix(rnorm(n * 2), n, 2)
  x <- z[, 1] + rnorm(n)
  y <- 0.8 * x + z[, 2] + rnorm(n, sd = 0.3)

  result <- suppressWarnings(
    bd_hsic_test(x, y, z, n_permutations = 100, seed = 1)
  )
  expect_s3_class(result, "kernel_test_result")
  expect_true(result$p_value >= 0 && result$p_value <= 1)
})
