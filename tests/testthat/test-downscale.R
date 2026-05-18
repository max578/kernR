test_that("kernel_downscale recovers a smooth coarse->fine map", {
  set.seed(1)
  n <- 120L
  coarse <- matrix(stats::runif(n * 2L, -1, 1), n, 2L,
                   dimnames = list(NULL, c("x1", "x2")))
  truth <- function(z) {
    2 * z[, "x1"] - 0.5 * z[, "x2"]^2 + sin(2 * z[, "x1"] * z[, "x2"])
  }
  fine <- matrix(truth(coarse) + stats::rnorm(n, sd = 0.05), ncol = 1L,
                 dimnames = list(NULL, "yield"))

  # Hold-out RMSE on a fresh sample
  n_new <- 60L
  new_coarse <- matrix(stats::runif(n_new * 2L, -1, 1), n_new, 2L,
                       dimnames = list(NULL, c("x1", "x2")))
  new_truth <- truth(new_coarse)

  fit <- kernel_downscale(coarse, fine, new_coarse)
  rmse <- sqrt(mean((fit$prediction[, 1L] - new_truth)^2))
  expect_s3_class(fit, "kernel_downscale")
  expect_equal(dim(fit$prediction), c(n_new, 1L))
  expect_lt(rmse, 0.3)  # generous; CV ridge should land well within this
})

test_that("kernel_downscale handles multi-output fine targets", {
  set.seed(2)
  n <- 80L
  coarse <- matrix(stats::runif(n * 2L), n, 2L,
                   dimnames = list(NULL, c("temp", "rain")))
  fine <- cbind(
    yield   = coarse[, "rain"] + stats::rnorm(n, sd = 0.1),
    biomass = coarse[, "temp"]^2 + stats::rnorm(n, sd = 0.1)
  )
  new_coarse <- matrix(stats::runif(20L * 2L), 20L, 2L,
                       dimnames = list(NULL, c("temp", "rain")))
  fit <- kernel_downscale(coarse, fine, new_coarse)
  expect_equal(dim(fit$prediction), c(20L, 2L))
  expect_equal(colnames(fit$prediction), c("yield", "biomass"))
})

test_that("kernel_downscale is reproducible across calls (deterministic given lambda='cv')", {
  set.seed(3)
  n <- 60L
  coarse <- matrix(stats::rnorm(n * 2L), n, 2L)
  fine <- matrix(coarse[, 1L] + stats::rnorm(n, sd = 0.1), ncol = 1L)
  new_coarse <- matrix(stats::rnorm(20L), 20L, 2L)
  f1 <- kernel_downscale(coarse, fine, new_coarse)
  f2 <- kernel_downscale(coarse, fine, new_coarse)
  expect_identical(f1$prediction, f2$prediction)
  expect_identical(f1$lambda,     f2$lambda)
})

test_that("kernel_downscale honours fixed lambda override", {
  set.seed(4)
  n <- 50L
  coarse <- matrix(stats::rnorm(n * 2L), n, 2L)
  fine <- matrix(coarse[, 1L] + stats::rnorm(n, sd = 0.1), ncol = 1L)
  fit <- kernel_downscale(coarse, fine, coarse, lambda = 0.01)
  expect_equal(fit$lambda, 0.01)
})

test_that("kernel_downscale returns weights when asked", {
  set.seed(5)
  n <- 30L
  coarse <- matrix(stats::rnorm(n * 2L), n, 2L)
  fine <- matrix(coarse[, 1L] + stats::rnorm(n, sd = 0.1), ncol = 1L)
  fit_no <- kernel_downscale(coarse, fine, coarse)
  fit_yes <- kernel_downscale(coarse, fine, coarse, return_weights = TRUE)
  expect_null(fit_no$weights)
  expect_equal(dim(fit_yes$weights), c(n, n))
  # Reconstruct prediction from weights
  expect_equal(fit_yes$weights %*% fine, fit_yes$prediction,
               ignore_attr = TRUE)
})

test_that("kernel_downscale validates inputs", {
  expect_error(
    kernel_downscale(matrix(0, 5L, 2L), matrix(0, 5L, 1L),
                     matrix(0, 5L, 2L)),
    "At least 10"
  )
  expect_error(
    kernel_downscale(matrix(0, 20L, 2L), matrix(0, 10L, 1L),
                     matrix(0, 5L, 2L)),
    "same number of rows"
  )
  expect_error(
    kernel_downscale(matrix(0, 20L, 2L), matrix(0, 20L, 1L),
                     matrix(0, 5L, 3L)),
    "same number of columns"
  )
})

test_that("print and as.data.frame work", {
  set.seed(6)
  n <- 30L
  coarse <- matrix(stats::rnorm(n * 2L), n, 2L)
  fine <- matrix(coarse[, 1L] + stats::rnorm(n, sd = 0.1), ncol = 1L,
                 dimnames = list(NULL, "y"))
  fit <- kernel_downscale(coarse, fine, coarse)
  expect_output(print(fit), "Kernel Downscaling")
  df <- as.data.frame(fit)
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), n)
})

test_that("fit_cme and predict.cme_fit are exported and round-trip", {
  set.seed(7)
  n <- 40L
  x <- matrix(stats::rnorm(n * 2L), n, 2L)
  y <- matrix(x[, 1L] + stats::rnorm(n, sd = 0.1), ncol = 1L)
  cme <- fit_cme(x, y)
  expect_s3_class(cme, "cme_fit")
  weights <- predict(cme, x[1:5, , drop = FALSE])
  expect_equal(dim(weights), c(5L, n))
})
