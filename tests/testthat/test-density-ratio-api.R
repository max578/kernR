# test-density-ratio-api.R -- public density-ratio API, weight diagnostics, and
# the conditional-mean-embedding print method. These exports previously had no
# direct coverage; the round-trips below pin their contracts.

# fit_density_ratio + predict_density_ratio (logistic backend) ----------------

test_that("fit_density_ratio + predict_density_ratio round-trip on the logistic backend", {
  set.seed(1L)
  n <- 200L
  z <- matrix(stats::rnorm(n * 2L), n, 2L)
  x <- stats::rbinom(n, 1L, stats::plogis(z[, 1L]))

  fit <- fit_density_ratio(x, z, method = "logistic", seed = 2L)
  expect_s3_class(fit, "density_ratio_fit")

  weights <- predict_density_ratio(fit, new_x = x, new_z = z, type = "weight")
  expect_length(weights, n)
  expect_true(all(is.finite(weights)))
  expect_true(all(weights > 0))

  log_ratio <- predict_density_ratio(fit, new_x = x, new_z = z, type = "log_ratio")
  expect_true(all(is.finite(log_ratio)))
  # The raw ratio is exp(log_ratio); the weight self-normalises to sum to n.
  ratio <- predict_density_ratio(fit, new_x = x, new_z = z, type = "ratio")
  expect_equal(ratio, exp(log_ratio), tolerance = 1e-8)
  expect_equal(sum(weights), n, tolerance = 1e-6)
})

test_that("fit_density_ratio is reproducible under a fixed seed", {
  set.seed(3L)
  n <- 150L
  z <- matrix(stats::rnorm(n * 2L), n, 2L)
  x <- stats::rbinom(n, 1L, stats::plogis(z[, 2L]))

  first <- predict_density_ratio(
    fit_density_ratio(x, z, method = "logistic", seed = 5L),
    new_x = x, new_z = z, type = "weight"
  )
  second <- predict_density_ratio(
    fit_density_ratio(x, z, method = "logistic", seed = 5L),
    new_x = x, new_z = z, type = "weight"
  )
  expect_identical(first, second)
})

# plot_weights diagnostic -----------------------------------------------------

test_that("plot_weights renders and returns the weights invisibly", {
  weights <- stats::runif(100L, 0.5, 1.5)
  pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  result <- plot_weights(weights)
  expect_identical(result, weights)
})

# print.cme_fit ---------------------------------------------------------------

test_that("print.cme_fit emits a compact summary and returns invisibly", {
  set.seed(7L)
  x <- matrix(stats::rnorm(80L), ncol = 2L)
  y <- matrix(x[, 1L] + stats::rnorm(40L, sd = 0.2), ncol = 1L)
  fit <- fit_cme(x, y, lambda = 1e-2)

  expect_output(print(fit), "Conditional mean embedding")
  expect_output(print(fit), "Training points: *40")
  expect_invisible(print(fit))
})
