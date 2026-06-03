test_that("weighted_hsic_stat matches the internal compiled engine", {
  set.seed(42)
  x <- matrix(rnorm(60), ncol = 2)
  y <- matrix(rnorm(60), ncol = 2)
  Kx <- kernel_matrix(x, kernel = resolve_bandwidth(kernel_spec(), x))
  Ky <- kernel_matrix(y, kernel = resolve_bandwidth(kernel_spec(), y))
  w <- runif(nrow(x))

  expect_equal(
    weighted_hsic_stat(Kx, Ky, w),
    weighted_hsic_stat_cpp(Kx, Ky, w)
  )
})

test_that("weighted_hsic_stat defaults to uniform weights", {
  set.seed(7)
  x <- matrix(rnorm(40), ncol = 2)
  Kx <- kernel_matrix(x, kernel = resolve_bandwidth(kernel_spec(), x))
  Ky <- kernel_matrix(x, kernel = resolve_bandwidth(kernel_spec(), x))

  expect_equal(
    weighted_hsic_stat(Kx, Ky),
    weighted_hsic_stat(Kx, Ky, rep(1, nrow(x)))
  )
})

test_that("weighted_hsic_stat validates its inputs", {
  Kx <- diag(4)
  expect_error(weighted_hsic_stat(matrix(1:6, nrow = 2), Kx), "square")
  expect_error(weighted_hsic_stat(Kx, diag(3)), "same dimensions")
  expect_error(weighted_hsic_stat(Kx, Kx, w = rep(1, 3)), "length")
  expect_error(weighted_hsic_stat(Kx, Kx, w = c(-1, 1, 1, 1)), "non-negative")
  expect_error(weighted_hsic_stat(Kx, Kx, w = rep(0, 4)), "positive sum")
})

test_that("resolve_bandwidth fills the median-heuristic bandwidth", {
  set.seed(1)
  x <- matrix(rnorm(40), ncol = 2)

  k <- resolve_bandwidth(kernel_spec(), x)
  expect_true(is.numeric(k$bandwidth))
  expect_gt(k$bandwidth, 0)

  # A fixed bandwidth passes through untouched.
  kf <- kernel_spec("rbf", bandwidth = 2.5)
  expect_identical(resolve_bandwidth(kf, x)$bandwidth, 2.5)
})
