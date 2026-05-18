# Tests for aggregate_downscale() (kernR 0.0.0.9015).
# Third downscaling method alongside kernel_downscale (CME) and
# dist_regression (bag-of-points): aggregate-likelihood inversion of
# a known aggregator with a GMM latent prior. Linear-Gaussian closed
# form + non-linear importance sampling.

.toy_prior_2c <- function() {
  list(
    means       = list(c(0, 0), c(2, 2)),
    covariances = list(diag(2L) * 0.5, diag(2L) * 0.5),
    weights     = c(0.5, 0.5)
  )
}

test_that("linear-Gaussian closed form: single-component recovers Kalman update", {
  set.seed(1L)
  A <- matrix(c(1, 1), nrow = 1L)
  prior_1c <- list(
    means = list(c(0, 0)),
    covariances = list(diag(2L)),
    weights = 1
  )
  y <- 2.0
  sigma_y <- 0.1
  fit <- aggregate_downscale(y, A, prior_1c, sigma_y = sigma_y)
  expect_equal(fit$method, "linear_closed_form")
  expect_equal(fit$aggregator_type, "linear")
  # Closed form: S_y = A I A^T + sigma_y^2 = 2 + 0.01. Kalman gain
  # K = I A^T S_y^{-1} = (1, 1)^T / 2.01. Posterior mean = (0, 0) +
  # K * (2 - 0) = (2/2.01, 2/2.01).
  expect_equal(fit$posterior_mean,
               c(2, 2) / 2.01, tolerance = 1e-9)
  expect_equal(length(fit$posterior_weights), 1L)
  expect_equal(fit$posterior_weights, 1)
})

test_that("linear-Gaussian closed form: 2-component prior reweights toward likely component", {
  set.seed(2L)
  A <- matrix(c(0.5, 0.5), nrow = 1L)         # spatial average
  prior <- .toy_prior_2c()
  # y near "average of 0,0" should favour the first component
  fit_low  <- aggregate_downscale(0.0, A, prior, sigma_y = 0.2)
  # y near "average of 2,2" should favour the second
  fit_high <- aggregate_downscale(2.0, A, prior, sigma_y = 0.2)
  expect_gt(fit_low$posterior_weights[1L],
            fit_low$posterior_weights[2L])
  expect_gt(fit_high$posterior_weights[2L],
            fit_high$posterior_weights[1L])
})

test_that("non-linear IS path runs and returns valid moments", {
  set.seed(3L)
  prior <- .toy_prior_2c()
  agg_fn <- function(x) matrix(sin(rowSums(x)), ncol = 1L)
  fit <- aggregate_downscale(0.5, agg_fn, prior, sigma_y = 0.2,
                             n_samples_per_component = 400L,
                             seed = 3L)
  expect_equal(fit$method, "nonlinear_is")
  expect_equal(fit$aggregator_type, "nonlinear")
  expect_length(fit$posterior_mean, 2L)
  expect_equal(dim(fit$posterior_cov), c(2L, 2L))
  expect_equal(sum(fit$posterior_weights), 1, tolerance = 1e-10)
  # Posterior covariance must be at least p.s.d. (eigenvalues >= -eps)
  expect_true(min(eigen(fit$posterior_cov, only.values = TRUE)$values)
              > -1e-8)
  # All component ESS finite and within (0, n_samples]
  expect_true(all(is.finite(fit$ess_per_component)))
  expect_true(all(fit$ess_per_component > 0 &
                  fit$ess_per_component <= 400))
})

test_that("non-linear IS path is reproducible under seed", {
  prior <- .toy_prior_2c()
  agg_fn <- function(x) matrix(rowSums(x^2), ncol = 1L)
  # The ESS-floor warning may incidentally fire on this toy
  # design; reproducibility is what this test asserts.
  f1 <- suppressWarnings(
    aggregate_downscale(1.0, agg_fn, prior, sigma_y = 0.3,
                        n_samples_per_component = 200L, seed = 7L)
  )
  f2 <- suppressWarnings(
    aggregate_downscale(1.0, agg_fn, prior, sigma_y = 0.3,
                        n_samples_per_component = 200L, seed = 7L)
  )
  expect_identical(f1$posterior_mean,    f2$posterior_mean)
  expect_identical(f1$posterior_weights, f2$posterior_weights)
})

test_that("ESS-floor warning fires when IS collapses (very narrow likelihood)", {
  set.seed(4L)
  prior <- .toy_prior_2c()
  agg_fn <- function(x) matrix(rowSums(x), ncol = 1L)
  # Very small sigma_y + y far from prior mass -> IS weights collapse
  expect_warning(
    aggregate_downscale(10, agg_fn, prior, sigma_y = 0.01,
                        n_samples_per_component = 50L,
                        min_ess_fraction = 0.5, seed = 4L),
    "effective sample size"
  )
})

test_that("posterior_sample_aggregate draws from the posterior mixture", {
  set.seed(5L)
  A <- matrix(c(0.5, 0.5), nrow = 1L)
  prior <- .toy_prior_2c()
  fit <- aggregate_downscale(2.0, A, prior, sigma_y = 0.3)
  samp <- posterior_sample_aggregate(fit, n = 1000L, seed = 5L)
  expect_equal(dim(samp), c(1000L, 2L))
  # Sample mean approaches posterior mean (loose tolerance for 1000)
  expect_equal(colMeans(samp), fit$posterior_mean, tolerance = 0.15)
})

test_that("input validation: bad sigma_y, bad aggregator dim, bad prior", {
  prior <- .toy_prior_2c()
  expect_error(
    aggregate_downscale(0, matrix(c(1, 1), 1L), prior, sigma_y = -1),
    "sigma_y"
  )
  # Aggregator with wrong cols
  expect_error(
    aggregate_downscale(c(0, 0), matrix(0, 2L, 3L), prior),
    "expected 2x2"
  )
  # Prior missing covariances
  expect_error(
    aggregate_downscale(0, matrix(c(1, 1), 1L),
                        list(means = list(c(0, 0)), weights = 1)),
    "covariances"
  )
  # Aggregator neither matrix nor function
  expect_error(
    aggregate_downscale(0, "not_a_matrix", prior),
    "matrix .linear. or a function"
  )
})

test_that("proxymix gmm_fit prior is accepted via slot extraction", {
  skip_if_not_installed("proxymix", minimum_version = "0.3.0")
  set.seed(6L)
  n <- 80L
  # Two-cluster training data
  X_train <- rbind(
    matrix(stats::rnorm(n * 2L, mean = 0), n, 2L),
    matrix(stats::rnorm(n * 2L, mean = 3), n, 2L)
  )
  target <- proxymix::gmm_target_from_samples(X_train)
  prior_fit <- proxymix::fit_proxymix(target, N = 2L, regime = "sample")
  A <- matrix(c(0.5, 0.5), nrow = 1L)
  fit <- aggregate_downscale(1.5, A, prior_fit, sigma_y = 0.2)
  expect_equal(fit$method, "linear_closed_form")
  expect_equal(fit$n_components, 2L)
  expect_equal(sum(fit$posterior_weights), 1, tolerance = 1e-10)
})
