test_that("propensity estimation produces valid scores", {
  set.seed(42)
  n <- 200
  x <- matrix(rnorm(n * 2), n, 2)
  logit_p <- 0.5 * x[, 1]
  t <- rbinom(n, 1, plogis(logit_p))

  ps <- estimate_propensity(t, x, cross_fit = FALSE)
  expect_s3_class(ps, "propensity_fit")
  expect_length(ps$scores, n)
  expect_true(all(ps$scores >= 0.01 & ps$scores <= 0.99))
})

test_that("cross-fitted propensity scores work", {
  set.seed(42)
  n <- 200
  x <- matrix(rnorm(n * 2), n, 2)
  t <- rbinom(n, 1, 0.5)

  ps <- estimate_propensity(t, x, cross_fit = TRUE, n_folds = 3)
  expect_length(ps$scores, n)
})

test_that("propensity rejects non-binary treatment", {
  expect_error(estimate_propensity(c(0, 1, 2), matrix(1, 3, 1)), "binary")
})

test_that("density ratio estimation produces valid output", {
  set.seed(42)
  n <- 150
  z <- matrix(rnorm(n * 2), n, 2)
  x <- z[, 1] + rnorm(n)

  dr <- estimate_density_ratio(x, z, method = "logistic", seed = 1)
  expect_s3_class(dr, "density_ratio_fit")
  expect_length(dr$weights, n)
  expect_true(dr$ess > 0 && dr$ess <= n)
  expect_true(all(dr$weights > 0))
})

test_that("overlap assessment works", {
  set.seed(42)
  scores <- c(runif(50, 0.3, 0.7), runif(50, 0.3, 0.7))
  treatment <- c(rep(1, 50), rep(0, 50))

  ov <- assess_overlap(scores, treatment)
  expect_s3_class(ov, "overlap_diagnostic")
  expect_false(ov$overlap_warning)
})

test_that("effective sample size is correct", {
  # Uniform weights -> ESS = n
  w <- rep(1, 100)
  expect_equal(effective_sample_size(w), 100)

  # All weight on one obs -> ESS = 1
  w2 <- c(100, rep(0, 99))
  expect_equal(effective_sample_size(w2), 1)
})
