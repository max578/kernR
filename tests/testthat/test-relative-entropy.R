test_that("relative_entropy is zero for identical Gaussians", {
  expect_equal(relative_entropy(0, 1, 0, 1), 0)
  expect_equal(relative_entropy(c(0, 0), diag(2), c(0, 0), diag(2)), 0)
})

test_that("relative_entropy matches the scalar closed form", {
  # KL(N(1,1) || N(0,1)) = (mu diff)^2 / 2 = 0.5
  expect_equal(relative_entropy(1, 1, 0, 1), 0.5, tolerance = 1e-10)
  # KL(N(0,4) || N(0,1)) = 0.5[ 4 - 1 + log(1/4) ] = 0.8068528
  expect_equal(relative_entropy(0, 4, 0, 1), 0.5 * (4 - 1 + log(1 / 4)),
               tolerance = 1e-10)
})

test_that("relative_entropy matches a vector closed form", {
  # mean shift of 1 in one coordinate, identity covariances -> 0.5
  expect_equal(relative_entropy(c(0, 0), diag(2), c(1, 0), diag(2)), 0.5,
               tolerance = 1e-10)
})

test_that("relative_entropy is non-negative and asymmetric", {
  set.seed(1)
  for (i in 1:20) {
    sp <- runif(1, 0.2, 3); sq <- runif(1, 0.2, 3)
    mp <- rnorm(1); mq <- rnorm(1)
    expect_gte(relative_entropy(mp, sp, mq, sq), 0)
  }
  # asymmetric in p, q when covariances differ
  expect_false(isTRUE(all.equal(relative_entropy(0, 1, 0, 4),
                                relative_entropy(0, 4, 0, 1))))
})

test_that("relative_entropy accepts scalar / vector / matrix covariance forms", {
  expect_equal(relative_entropy(c(0, 0), 1, c(0, 0), 1), 0)          # scalar
  expect_equal(relative_entropy(c(0, 0), c(1, 1), c(0, 0), diag(2)), 0)  # vec vs mat
})

test_that("relative_entropy errors on a non-positive-definite reference", {
  expect_error(relative_entropy(c(0, 0), diag(2), c(0, 0),
                                matrix(c(1, 2, 2, 1), 2)),
               "positive definite")
  expect_error(relative_entropy(c(0, 0), diag(2), 0, 1),
               "same length")
})

test_that("relative_entropy_ensemble approximates the closed form", {
  set.seed(42)
  n <- 40000L
  p <- matrix(rnorm(n * 2), ncol = 2)              # ~ N(0, I)
  q <- matrix(rnorm(n * 2), ncol = 2)
  q[, 1] <- q[, 1] + 1                             # ~ N((1,0), I)
  closed <- relative_entropy(c(0, 0), diag(2), c(1, 0), diag(2))
  expect_equal(relative_entropy_ensemble(p, q), closed, tolerance = 0.05)
})

test_that("cir_objective integrates the divergence profile correctly", {
  # constant divergence M over [0, L] -> CIR = L (the full horizon)
  expect_equal(cir_objective(seq(0, 5, by = 1), rep(2, 6)), 5, tolerance = 1e-10)
  # no recoverable future information (M = 0) -> 0
  expect_equal(cir_objective(c(0, 1, 2), c(0, 0, 0)), 0)
  # exponential decay D(L) = M exp(-r L): CIR -> 1/r as the grid extends
  lag <- seq(0, 40, by = 0.1); r <- 0.6
  expect_equal(cir_objective(lag, 2 * exp(-r * lag)), 1 / r, tolerance = 1e-2)
})

test_that("cir_objective validates its inputs", {
  expect_error(cir_objective(1, 1), "equal length >= 2")
  expect_error(cir_objective(c(2, 1), c(1, 1)), "increasing")
})
