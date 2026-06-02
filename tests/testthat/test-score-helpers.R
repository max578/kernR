# Tests for numeric_score(): the finite-difference score adapter for ksd_test.

test_that("numeric_score() matches the closed-form standard-normal score", {
  set.seed(30)
  x <- matrix(stats::rnorm(200L), ncol = 2L)
  ld <- function(z) -0.5 * rowSums(z^2)
  s <- numeric_score(ld)
  expect_lt(max(abs(s(x) - (-x))), 1e-6)
})

test_that("an additive constant in the log density cancels in the score", {
  set.seed(31)
  x <- matrix(stats::rnorm(100L), ncol = 2L)
  ld <- function(z) -0.5 * rowSums(z^2)
  ld_shift <- function(z) -0.5 * rowSums(z^2) + 42
  expect_equal(numeric_score(ld)(x), numeric_score(ld_shift)(x))
})

test_that("numeric_score() drives ksd_test on a well-specified target", {
  set.seed(32)
  x <- matrix(stats::rnorm(400L), ncol = 2L)
  ld <- function(z) -0.5 * rowSums(z^2)
  fit <- ksd_test(x, score = numeric_score(ld), n_boot = 199L, seed = 1L)
  expect_gt(fit$p_value, 0.05)
})

test_that("numeric_score() drives ksd_test to reject a mis-specified target", {
  set.seed(33)
  x <- matrix(stats::rnorm(400L), ncol = 2L) + 0.8
  ld <- function(z) -0.5 * rowSums(z^2)
  fit <- ksd_test(x, score = numeric_score(ld), n_boot = 199L, seed = 1L)
  expect_lt(fit$p_value, 0.05)
})

test_that("numeric_score() agrees with gaussian_score() for a correlated target", {
  set.seed(34)
  x <- matrix(stats::rnorm(100L), ncol = 2L)
  sig <- matrix(c(1, 0.5, 0.5, 1), nrow = 2L)
  prec <- solve(sig)
  ld <- function(z) -0.5 * rowSums((z %*% prec) * z)
  s_num <- numeric_score(ld)(x)
  s_cf <- gaussian_score(sigma = sig)(x)
  expect_lt(max(abs(s_num - s_cf)), 1e-5)
})

test_that("numeric_score() validates its arguments and the log-density output", {
  expect_error(numeric_score("not a function"), "must be a function")
  expect_error(numeric_score(function(z) z, h = 0), "positive number")
  x <- matrix(stats::rnorm(40L), ncol = 2L)
  bad_len <- numeric_score(function(z) 1)
  expect_error(bad_len(x), "length nrow")
})
