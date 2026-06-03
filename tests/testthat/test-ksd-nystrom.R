# test-ksd-nystrom.R -- low-rank (Nystrom) kernel Stein discrepancy test.
#
# Anchors ksd_test_nystrom() against the exact ksd_test(). Because the Stein
# kernel is positive semi-definite, FF^t is itself a valid Stein kernel, so the
# accelerated test is the exact procedure on a rank-m kernel: the reconstruction
# test pins the statistic, and the Type-I test confirms the degenerate-U-
# statistic wild bootstrap keeps its level under the approximation -- the one
# property the acceleration could have broken.

test_that("nystrom KSD reconstructs the exact statistic as m -> n (IMQ)", {
  set.seed(1)
  x <- matrix(stats::rnorm(200L), ncol = 2L)
  ex <- ksd_test(x, n_boot = 1L, seed = 1L)
  ap <- ksd_test_nystrom(x, m = nrow(x) - 1L, n_boot = 1L, seed = 1L)
  # KSD is a near-zero degenerate U-statistic under the null, so compare on
  # the absolute scale: the factorisation residual is O(1e-9).
  expect_lt(abs(ap$statistic - ex$statistic), 1e-6)
})

test_that("nystrom KSD reconstructs the exact statistic as m -> n (RBF)", {
  set.seed(2)
  x <- matrix(stats::rnorm(200L), ncol = 2L)
  ex <- ksd_test(x, kernel = "rbf", n_boot = 1L, seed = 1L)
  ap <- ksd_test_nystrom(x, kernel = "rbf", m = nrow(x) - 1L,
                         n_boot = 1L, seed = 1L)
  expect_lt(abs(ap$statistic - ex$statistic), 1e-6)
})

test_that("nystrom KSD is reproducible with seed", {
  set.seed(3)
  x <- matrix(stats::rnorm(800L), ncol = 2L)
  f1 <- ksd_test_nystrom(x, m = 60L, n_boot = 99L, seed = 7L)
  f2 <- ksd_test_nystrom(x, m = 60L, n_boot = 99L, seed = 7L)
  expect_identical(f1$statistic, f2$statistic)
  expect_identical(f1$p_value, f2$p_value)
  expect_identical(f1$null_distribution, f2$null_distribution)
})

test_that("a mean-shifted sample is rejected (power)", {
  set.seed(4)
  x <- matrix(stats::rnorm(600L), ncol = 2L) + 1
  fit <- ksd_test_nystrom(x, m = 60L, n_boot = 199L, seed = 4L)
  expect_true(fit$reject)
  expect_gt(fit$surprise_bits, 0)
})

test_that("a custom Gaussian-score target is honoured", {
  set.seed(5)
  sig <- matrix(c(1, 0.6, 0.6, 1), nrow = 2L)
  x_cor <- matrix(stats::rnorm(600L), ncol = 2L) %*% chol(sig)
  # Correct target: not rejected; wrong (standard-normal) target: rejected.
  fit_right <- ksd_test_nystrom(x_cor, score = gaussian_score(sigma = sig),
                                m = 60L, n_boot = 199L, seed = 5L)
  fit_wrong <- ksd_test_nystrom(x_cor, m = 60L, n_boot = 199L, seed = 5L)
  expect_false(fit_right$reject)
  expect_true(fit_wrong$reject)
})

test_that("the object class and fields match the exact test", {
  set.seed(6)
  x <- matrix(stats::rnorm(400L), ncol = 2L)
  fit <- ksd_test_nystrom(x, m = 40L, n_boot = 49L, seed = 6L)
  expect_s3_class(fit, c("ksd_test", "kernel_test_result"))
  expect_identical(fit$approximation, "nystrom")
  expect_identical(fit$stein_kernel, "imq")
  expect_true(is.finite(fit$surprise_bits))
  expect_identical(fit$m, 40L)
})

test_that("input validation guards small n, bad m, and bad method", {
  x_small <- matrix(stats::rnorm(16L), ncol = 2L)         # n = 8 < 10
  expect_error(ksd_test_nystrom(x_small, m = 5L), "at least 10")

  x <- matrix(stats::rnorm(200L), ncol = 2L)
  expect_error(ksd_test_nystrom(x, m = 1L), "integer >= 2")
  expect_error(ksd_test_nystrom(x, method = "rff", m = 20L),
               "should be")
  expect_error(
    ksd_test_nystrom(x, score = function(z) z[, 1L, drop = FALSE], m = 20L),
    "same dimensions"
  )
})

test_that("m is capped at n - 1", {
  set.seed(8)
  x <- matrix(stats::rnorm(60L), ncol = 2L)               # n = 30
  fit <- ksd_test_nystrom(x, m = 500L, n_boot = 49L, seed = 8L)
  expect_identical(fit$m, 29L)
})

test_that("Type-I error rate stays near the nominal level under the null", {
  skip_on_cran()
  set.seed(99)
  alpha <- 0.1
  n_trials <- 80L
  reject <- logical(n_trials)
  for (t in seq_len(n_trials)) {
    x <- matrix(stats::rnorm(600L), ncol = 2L)             # n = 300, ~ N(0, I)
    fit <- ksd_test_nystrom(x, m = 50L, alpha = alpha, n_boot = 99L)
    reject[t] <- fit$reject
  }
  # FF^t is a valid Stein kernel, so the degenerate-U wild bootstrap keeps its
  # level; expect the empirical size near alpha (slack for 80 Monte-Carlo runs).
  expect_lt(mean(reject), alpha + 0.12)
})
