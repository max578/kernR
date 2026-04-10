test_that("kernel_spec creates valid objects", {
  k <- kernel_spec()
  expect_s3_class(k, "kernel_spec")
  expect_equal(k$type, "rbf")
  expect_equal(k$bandwidth, "median")

  k2 <- kernel_spec("linear")
  expect_equal(k2$type, "linear")

  k3 <- kernel_spec("matern", nu = 1.5)
  expect_equal(k3$nu, 1.5)

  expect_error(kernel_spec("rbf", bandwidth = -1))
  expect_error(kernel_spec("matern", nu = -1))
})

test_that("kernel_matrix produces correct dimensions", {
  x <- matrix(rnorm(100), 50, 2)
  K <- kernel_matrix(x)
  expect_equal(dim(K), c(50, 50))

  y <- matrix(rnorm(60), 30, 2)
  K2 <- kernel_matrix(x, y)
  expect_equal(dim(K2), c(50, 30))
})

test_that("RBF kernel matrix is symmetric and positive semi-definite", {
  set.seed(42)
  x <- matrix(rnorm(100), 50, 2)
  K <- kernel_matrix(x, kernel = kernel_spec("rbf", bandwidth = 1.0))

  # Symmetric

  expect_equal(K, t(K), tolerance = 1e-10)

  # Diagonal is 1
  expect_equal(diag(K), rep(1, 50), tolerance = 1e-10)

  # PSD: all eigenvalues >= 0
  eigs <- eigen(K, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eigs > -1e-10))
})

test_that("linear kernel is x %*% t(y)", {
  set.seed(1)
  x <- matrix(rnorm(30), 10, 3)
  y <- matrix(rnorm(15), 5, 3)
  K <- kernel_matrix(x, y, kernel = kernel_spec("linear"))
  expect_equal(K, x %*% t(y), tolerance = 1e-10)
})

test_that("polynomial kernel computes correctly", {
  set.seed(1)
  x <- matrix(rnorm(20), 10, 2)
  K <- kernel_matrix(x, kernel = kernel_spec("polynomial", degree = 2, offset = 1))
  expected <- (x %*% t(x) + 1)^2
  expect_equal(K, expected, tolerance = 1e-10)
})

test_that("median bandwidth is positive", {
  set.seed(1)
  x <- matrix(rnorm(200), 100, 2)
  bw <- select_bandwidth(x, "median")
  expect_true(is.numeric(bw))
  expect_true(bw > 0)
})

test_that("kernel_matrix rejects mismatched dimensions", {
  x <- matrix(1, 10, 2)
  y <- matrix(1, 10, 3)
  expect_error(kernel_matrix(x, y), "same number of columns")
})
