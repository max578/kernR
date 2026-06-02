# Tests for the kernel Stein discrepancy goodness-of-fit test (ksd_test) and
# the gaussian_score() factory. Covers Type-I behaviour, power against mean /
# variance / tail mis-specification, the score factory, seed reproducibility,
# both base kernels, the result contract, and input validation.

# Well-specified behaviour ---------------------------------------------------

test_that("a well-specified sample is not rejected (IMQ and RBF)", {
  set.seed(10)
  x <- matrix(stats::rnorm(600L), ncol = 2L)

  fit_imq <- ksd_test(x, n_boot = 299L, seed = 1L)
  fit_rbf <- ksd_test(x, kernel = "rbf", n_boot = 299L, seed = 1L)

  expect_s3_class(fit_imq, c("ksd_test", "kernel_test_result"))
  expect_gt(fit_imq$p_value, 0.05)
  expect_gt(fit_rbf$p_value, 0.05)
  expect_false(fit_imq$reject)

  # KSD U-statistic sits near zero under the null
  expect_lt(abs(fit_imq$statistic), 0.05)
})

# Power ----------------------------------------------------------------------

test_that("a mean-shifted sample is rejected", {
  set.seed(11)
  x <- matrix(stats::rnorm(600L), ncol = 2L) + 0.7
  fit <- ksd_test(x, n_boot = 299L, seed = 1L)

  expect_lt(fit$p_value, 0.05)
  expect_true(fit$reject)
  expect_gt(fit$statistic, 0)
})

test_that("a variance-inflated sample is rejected", {
  set.seed(12)
  x <- matrix(stats::rnorm(600L, sd = 1.8), ncol = 2L)
  fit <- ksd_test(x, n_boot = 299L, seed = 1L)

  expect_lt(fit$p_value, 0.05)
  expect_true(fit$reject)
})

test_that("a heavy-tailed (t3) sample is rejected against a normal target", {
  set.seed(13)
  x <- matrix(stats::rt(600L, df = 3), ncol = 2L)
  fit <- ksd_test(x, n_boot = 299L, seed = 1L)

  expect_lt(fit$p_value, 0.05)
})

# gaussian_score() factory ---------------------------------------------------

test_that("gaussian_score() defaults to the standard-normal score -x", {
  x <- matrix(stats::rnorm(40L), ncol = 2L)
  s <- gaussian_score()
  expect_equal(s(x), -x)
})

test_that("gaussian_score() matches a non-standard target it was built for", {
  set.seed(14)
  sig <- matrix(c(1, 0.6, 0.6, 1), nrow = 2L)
  x <- matrix(stats::rnorm(600L), ncol = 2L) %*% chol(sig)

  fit_match <- ksd_test(x, score = gaussian_score(sigma = sig),
                        n_boot = 299L, seed = 1L)
  expect_gt(fit_match$p_value, 0.05)

  # The wrong target (identity covariance) is rejected on the same sample
  fit_wrong <- ksd_test(x, n_boot = 299L, seed = 1L)
  expect_lt(fit_wrong$p_value, 0.05)
})

test_that("gaussian_score() validates mean and covariance dimensions", {
  x <- matrix(stats::rnorm(40L), ncol = 2L)
  expect_error(gaussian_score(mean = c(0, 0, 0))(x), "length 2")
  expect_error(gaussian_score(sigma = diag(3))(x), "2 x 2")
  asymmetric <- matrix(c(1, 0.3, 0.6, 1), nrow = 2L)
  expect_error(gaussian_score(sigma = asymmetric)(x), "symmetric")
  singular <- matrix(c(1, 1, 1, 1), nrow = 2L)
  expect_error(gaussian_score(sigma = singular)(x), "invertible")
})

# Reproducibility and dispatch -----------------------------------------------

test_that("the statistic is seed-independent and the p-value is reproducible", {
  set.seed(15)
  x <- matrix(stats::rnorm(400L), ncol = 2L)

  a <- ksd_test(x, n_boot = 199L, seed = 42L)
  b <- ksd_test(x, n_boot = 199L, seed = 42L)
  c <- ksd_test(x, n_boot = 199L, seed = 7L)

  # Same seed -> identical p-value
  expect_identical(a$p_value, b$p_value)
  # The observed statistic does not depend on the bootstrap seed
  expect_identical(a$statistic, c$statistic)
})

test_that("a univariate (vector) sample is handled", {
  set.seed(16)
  x <- stats::rnorm(300L)
  fit <- ksd_test(x, n_boot = 199L, seed = 1L)
  expect_equal(fit$dimension, 1L)
  expect_gt(fit$p_value, 0.05)
})

test_that("a numeric bandwidth is honoured and the print method runs", {
  set.seed(17)
  x <- matrix(stats::rnorm(200L), ncol = 2L)
  fit <- ksd_test(x, bandwidth = 1.5, n_boot = 99L, seed = 1L)
  expect_equal(fit$bandwidth, 1.5)
  expect_output(print(fit), "Goodness-of-fit verdict")
  expect_output(print(fit), "Stein kernel")
})

# Input validation -----------------------------------------------------------

test_that("ksd_test() rejects malformed input", {
  x <- matrix(stats::rnorm(200L), ncol = 2L)

  expect_error(ksd_test(matrix(1:6, ncol = 2L)), "at least 5 rows")
  bad <- x
  bad[1L, 1L] <- NA_real_
  expect_error(ksd_test(bad), "finite")
  expect_error(ksd_test(x, beta = 0), "in \\(-1, 0\\)")
  expect_error(ksd_test(x, beta = -1.5), "in \\(-1, 0\\)")
  expect_error(ksd_test(x, alpha = 0), "in \\(0, 1\\)")
  expect_error(ksd_test(x, n_boot = 0L), "positive integer")
  expect_error(ksd_test(x, bandwidth = -1), "positive number")
})

test_that("ksd_test() validates the score argument and its output", {
  x <- matrix(stats::rnorm(200L), ncol = 2L)

  expect_error(ksd_test(x, score = "not a function"), "must be `NULL` or a function")
  expect_error(
    ksd_test(x, score = function(z) z[, 1L, drop = FALSE]),
    "same dimensions"
  )
  expect_error(
    ksd_test(x, score = function(z) {
      out <- -z
      out[1L, 1L] <- Inf
      out
    }),
    "non-finite"
  )
})
