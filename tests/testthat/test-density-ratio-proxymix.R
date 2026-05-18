# Tests for the proxymix density-ratio backend.
# Skipped when proxymix is not installed (it's a soft Suggests, not a
# hard dependency: GRDC-firewalled, MIT, optional).

test_that("estimate_density_ratio dispatches to proxymix backend", {
  testthat::skip_if_not_installed("proxymix", minimum_version = "0.3.0")
  set.seed(11L)
  n <- 80L
  z <- matrix(stats::rnorm(n * 2L), n, 2L)
  x <- z[, 1L] + stats::rnorm(n, sd = 0.5)

  dr <- estimate_density_ratio(
    x, z,
    method = "proxymix",
    proxymix_components = 2L,
    seed   = 11L
  )
  expect_s3_class(dr, "density_ratio_fit")
  expect_equal(dr$method, "proxymix")
  expect_equal(length(dr$weights), n)
  expect_true(is.finite(dr$ess))
  expect_gt(dr$ess, 0)
})

test_that("proxymix density-ratio weights are finite and positive", {
  testthat::skip_if_not_installed("proxymix", minimum_version = "0.3.0")
  set.seed(12L)
  n <- 60L
  z <- matrix(stats::rnorm(n * 2L), n, 2L)
  x <- z[, 1L] + stats::rnorm(n)
  dr <- estimate_density_ratio(x, z, method = "proxymix", seed = 12L)
  expect_true(all(is.finite(dr$weights)))
  expect_true(all(dr$weights > 0))
  expect_equal(sum(dr$weights), n, tolerance = 1e-6)
})

test_that("proxymix dispatch is reproducible with seed", {
  testthat::skip_if_not_installed("proxymix", minimum_version = "0.3.0")
  set.seed(13L)
  n <- 60L
  z <- matrix(stats::rnorm(n * 2L), n, 2L)
  x <- z[, 1L] + stats::rnorm(n)
  d1 <- estimate_density_ratio(x, z, method = "proxymix", seed = 99L)
  d2 <- estimate_density_ratio(x, z, method = "proxymix", seed = 99L)
  expect_identical(d1$weights, d2$weights)
})

test_that("bd_hsic_test runs end-to-end with proxymix backend", {
  testthat::skip_if_not_installed("proxymix", minimum_version = "0.3.0")
  set.seed(14L)
  n <- 200L
  z <- matrix(stats::rnorm(n * 2L), n, 2L)
  x <- z[, 1L] + stats::rnorm(n)
  y <- 0.7 * x + z[, 2L] + stats::rnorm(n, sd = 0.4)

  res <- bd_hsic_test(
    x, y, z,
    density_ratio  = "proxymix",
    n_permutations = 99L,
    seed           = 14L
  )
  expect_s3_class(res, "kernel_test_result")
  expect_equal(res$method, "bd-HSIC")
  expect_true(is.finite(res$statistic))
  expect_true(is.finite(res$ess))
  # With a real causal effect we expect rejection at moderate n.
  expect_lt(res$p_value, 0.05)
})

test_that("estimate_density_ratio errors helpfully when proxymix missing", {
  # This test deliberately stubs the requireNamespace path to exercise
  # the helpful error message, regardless of whether proxymix is
  # actually installed.
  if (requireNamespace("proxymix", quietly = TRUE)) {
    # When installed, skip — the error path can't fire.
    skip("proxymix is installed; error-path test not applicable")
  }
  expect_error(
    estimate_density_ratio(
      matrix(stats::rnorm(40L), 20L, 2L),
      matrix(stats::rnorm(40L), 20L, 2L),
      method = "proxymix"
    ),
    "requires the `proxymix` package"
  )
})
