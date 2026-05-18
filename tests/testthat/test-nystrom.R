test_that("nystrom_factor reconstructs the full kernel as m -> n", {
  set.seed(1)
  n <- 60L
  x <- matrix(stats::rnorm(n * 2L), n, 2L)
  k <- kernel_spec("rbf", bandwidth = 1.0)
  K_full <- kernel_matrix(x, kernel = k)

  fit <- nystrom_factor(x, kernel = k, m = n - 1L, seed = 1L)
  K_approx <- tcrossprod(fit$F)
  # Frobenius error normalised by full norm: tight at near-full rank.
  rel_err <- sqrt(sum((K_full - K_approx)^2)) / sqrt(sum(K_full^2))
  expect_lt(rel_err, 0.05)
})

test_that("nystrom_factor caps m at n - 1 and reports effective rank", {
  set.seed(2)
  x <- matrix(stats::rnorm(40L), ncol = 2L)
  f <- nystrom_factor(x, m = 200L, seed = 2L)
  expect_equal(f$m, nrow(x) - 1L)
  expect_equal(ncol(f$F), f$m)
})

test_that("nystrom_factor is reproducible with seed", {
  set.seed(3)
  x <- matrix(stats::rnorm(100L), ncol = 2L)
  f1 <- nystrom_factor(x, m = 20L, seed = 7L)
  f2 <- nystrom_factor(x, m = 20L, seed = 7L)
  expect_identical(f1$F, f2$F)
})

test_that("nystrom_factor validates inputs", {
  x <- matrix(stats::rnorm(40L), ncol = 2L)
  expect_error(nystrom_factor(x, m = 1L), "m")
  expect_error(nystrom_factor(x, regularise = -1), "regularise")
})

test_that("rff_features approximates the RBF kernel at moderate D", {
  set.seed(4)
  n <- 80L
  x <- matrix(stats::rnorm(n * 2L), n, 2L)
  k <- kernel_spec("rbf", bandwidth = 1.0)
  K_full <- kernel_matrix(x, kernel = k)

  f <- rff_features(x, kernel = k, D = 1500L, seed = 4L)
  K_approx <- tcrossprod(f$F)
  rel_err <- sqrt(sum((K_full - K_approx)^2)) / sqrt(sum(K_full^2))
  expect_lt(rel_err, 0.1)
})

test_that("rff_features rejects non-RBF kernels", {
  x <- matrix(stats::rnorm(40L), ncol = 2L)
  expect_error(
    rff_features(x, kernel = kernel_spec("matern", nu = 2.5)),
    "RBF"
  )
})

test_that("rff_features is reproducible with seed", {
  x <- matrix(stats::rnorm(60L), ncol = 2L)
  f1 <- rff_features(x, D = 50L, seed = 11L)
  f2 <- rff_features(x, D = 50L, seed = 11L)
  expect_identical(f1$F, f2$F)
})

test_that("hsic_test_nystrom detects dependence and respects independence", {
  set.seed(5)
  n <- 400L
  xx <- stats::rnorm(n)
  y_dep <- xx + stats::rnorm(n, sd = 0.5)
  y_ind <- stats::rnorm(n)

  res_dep <- hsic_test_nystrom(xx, y_dep, m = 60L,
                               n_permutations = 99L, seed = 1L)
  res_ind <- hsic_test_nystrom(xx, y_ind, m = 60L,
                               n_permutations = 99L, seed = 1L)
  expect_s3_class(res_dep, "kernel_test_result")
  expect_equal(res_dep$approximation, "nystrom")
  expect_lt(res_dep$p_value, 0.05)
  expect_gt(res_ind$p_value, 0.05)
})

test_that("hsic_test_nystrom with method='rff' detects dependence", {
  set.seed(6)
  n <- 400L
  xx <- stats::rnorm(n)
  y  <- xx^2 + stats::rnorm(n, sd = 0.3)
  res <- hsic_test_nystrom(xx, y, method = "rff", m = 150L,
                           n_permutations = 99L, seed = 1L)
  expect_equal(res$approximation, "rff")
  expect_lt(res$p_value, 0.05)
})

test_that("hsic_test_nystrom is reproducible with seed", {
  set.seed(7)
  n <- 200L
  xx <- stats::rnorm(n); y <- xx + stats::rnorm(n, sd = 0.5)
  r1 <- hsic_test_nystrom(xx, y, m = 40L, n_permutations = 50L, seed = 9L)
  r2 <- hsic_test_nystrom(xx, y, m = 40L, n_permutations = 50L, seed = 9L)
  expect_identical(r1$statistic, r2$statistic)
  expect_identical(r1$p_value,   r2$p_value)
})

test_that("hsic_test_nystrom validates inputs", {
  xx <- stats::rnorm(50L)
  expect_error(hsic_test_nystrom(xx, stats::rnorm(40L)), "same number")
  expect_error(hsic_test_nystrom(stats::rnorm(5L), stats::rnorm(5L)),
               "At least 10")
})

test_that("hsic_test_nystrom and hsic_test agree at high rank", {
  # Sanity check: at high enough rank, Nystrom should give a similar
  # p-value verdict on the same data as the exact HSIC test.
  set.seed(8)
  n <- 150L
  xx <- stats::rnorm(n)
  y  <- xx + stats::rnorm(n, sd = 0.5)
  exact <- hsic_test(xx, y, n_permutations = 199L, seed = 1L)
  nys   <- hsic_test_nystrom(xx, y, m = 100L,
                             n_permutations = 199L, seed = 1L)
  expect_lt(exact$p_value, 0.05)
  expect_lt(nys$p_value, 0.05)
})

test_that("m_x / m_y override m per factor", {
  set.seed(9)
  n <- 200L
  xx <- matrix(stats::rnorm(n * 2L), n, 2L)
  y  <- xx[, 1L] + stats::rnorm(n, sd = 0.5)
  res <- hsic_test_nystrom(xx, y, m_x = 30L, m_y = 60L,
                           n_permutations = 50L, seed = 1L)
  expect_equal(res$m_x, 30L)
  expect_equal(res$m_y, 60L)
})
