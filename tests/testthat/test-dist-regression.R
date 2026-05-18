make_bags_mean <- function(M, n_per = 40L, sigma = 1, seed = 1L) {
  set.seed(seed)
  mu <- stats::runif(M, -3, 3)
  bags <- lapply(mu, function(m) {
    matrix(stats::rnorm(n_per, mean = m, sd = sigma), ncol = 1L)
  })
  list(bags = bags, mu = mu)
}

make_bags_var <- function(M, n_per = 40L, seed = 1L) {
  set.seed(seed)
  s <- stats::runif(M, 0.5, 2.5)
  bags <- lapply(s, function(sigma) {
    matrix(stats::rnorm(n_per, mean = 0, sd = sigma), ncol = 1L)
  })
  list(bags = bags, s = s)
}


test_that("dist_regression with linear outer recovers Y = bag mean", {
  d <- make_bags_mean(M = 60L, n_per = 50L, seed = 2L)
  fit <- dist_regression(d$bags, y = d$mu, outer = "linear")
  expect_s3_class(fit, "dist_regression")
  expect_equal(fit$d_y, 1L)
  expect_equal(fit$outer, "linear")

  # Predict on fresh bags
  d_new <- make_bags_mean(M = 20L, n_per = 50L, seed = 3L)
  pred <- predict(fit, d_new$bags)
  rmse <- sqrt(mean((pred - d_new$mu)^2))
  expect_lt(rmse, 0.5)
})

test_that("dist_regression with rbf outer recovers Y = bag SD (distributional)", {
  # Mean-zero bags with varying SDs: Y = bag SD is invisible to a linear
  # mean-embedding outer kernel (means are all zero), but an RBF outer
  # over embeddings *can* pick up second-moment structure indirectly.
  # We test the weaker claim that RBF gives lower RMSE than linear here.
  d <- make_bags_var(M = 80L, n_per = 60L, seed = 4L)
  fit_lin <- dist_regression(d$bags, y = d$s, outer = "linear")
  fit_rbf <- dist_regression(d$bags, y = d$s, outer = "rbf")

  d_new <- make_bags_var(M = 20L, n_per = 60L, seed = 5L)
  rmse_lin <- sqrt(mean((predict(fit_lin, d_new$bags) - d_new$s)^2))
  rmse_rbf <- sqrt(mean((predict(fit_rbf, d_new$bags) - d_new$s)^2))

  expect_true(is.finite(rmse_lin) && is.finite(rmse_rbf))
  # RBF should be no worse than linear on this distributional target.
  expect_lte(rmse_rbf, rmse_lin + 0.05)
})

test_that("dist_regression is reproducible given lambda='cv'", {
  d <- make_bags_mean(M = 30L, n_per = 30L, seed = 7L)
  f1 <- dist_regression(d$bags, y = d$mu, outer = "linear")
  f2 <- dist_regression(d$bags, y = d$mu, outer = "linear")
  expect_identical(f1$alpha, f2$alpha)
  expect_identical(f1$lambda, f2$lambda)
})

test_that("dist_regression honours fixed lambda override", {
  d <- make_bags_mean(M = 20L, n_per = 25L, seed = 8L)
  fit <- dist_regression(d$bags, y = d$mu, lambda = 0.01)
  expect_equal(fit$lambda, 0.01)
})

test_that("dist_regression supports variable bag sizes", {
  set.seed(9L)
  M <- 25L
  mu <- stats::runif(M, -2, 2)
  bags <- lapply(mu, function(m) {
    n_i <- sample(15:60, 1L)
    matrix(stats::rnorm(n_i, mean = m), ncol = 1L)
  })
  fit <- dist_regression(bags, y = mu, outer = "linear")
  expect_equal(fit$M, M)

  bag_sizes <- vapply(bags, nrow, integer(1L))
  expect_true(min(bag_sizes) < max(bag_sizes))
  # Prediction works
  pred <- predict(fit, bags[1:5])
  expect_equal(length(pred), 5L)
})

test_that("dist_regression validates inputs", {
  expect_error(
    dist_regression(list(matrix(1, 5L, 1L)),
                    y = 1),
    "at least"
  )
  d <- make_bags_mean(M = 10L, n_per = 20L, seed = 10L)
  expect_error(
    dist_regression(d$bags, y = d$mu[-1L]),
    "length"
  )
  bags_mismatch <- d$bags
  bags_mismatch[[1L]] <- matrix(1, 5L, 2L)
  expect_error(
    dist_regression(bags_mismatch, y = d$mu),
    "same number of columns"
  )
  expect_error(
    dist_regression(d$bags, y = d$mu, outer = "rbf",
                    outer_bandwidth = -1),
    "outer_bandwidth"
  )
})

test_that("predict.dist_regression validates new bag dim", {
  d <- make_bags_mean(M = 20L, n_per = 30L, seed = 11L)
  fit <- dist_regression(d$bags, y = d$mu, outer = "linear")
  expect_error(
    predict(fit, list(matrix(1, 10L, 2L))),
    "same number of columns"
  )
})

test_that("dist_regression with rbf outer + median bandwidth resolves positively", {
  d <- make_bags_mean(M = 30L, n_per = 30L, seed = 12L)
  fit <- dist_regression(d$bags, y = d$mu, outer = "rbf")
  expect_true(fit$outer_bandwidth > 0)
  expect_true(is.finite(fit$outer_bandwidth))
})

test_that("print and as.data.frame work", {
  d <- make_bags_mean(M = 15L, n_per = 25L, seed = 13L)
  fit <- dist_regression(d$bags, y = d$mu, outer = "rbf")
  expect_output(print(fit), "Distribution Regression")
  df <- as.data.frame(fit)
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 15L)
  expect_setequal(colnames(df),
                  c("bag_id", "bag_size", "y_train", "y_fit"))
})

test_that("dist_regression accepts multivariate y", {
  set.seed(14L)
  M <- 25L
  mu <- stats::runif(M, -2, 2)
  bags <- lapply(mu, function(m) {
    matrix(stats::rnorm(40L, mean = m), ncol = 1L)
  })
  y_mat <- cbind(target1 = mu, target2 = mu^2)
  fit <- dist_regression(bags, y = y_mat, outer = "linear")
  expect_equal(fit$d_y, 2L)
  pred <- predict(fit, bags[1:5])
  expect_equal(dim(pred), c(5L, 2L))
})
