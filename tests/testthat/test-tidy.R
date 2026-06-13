# Tidiers for kernR verdict objects (broom-style tidy() S3 methods).

test_that("tidy() on a kernel_test_result returns the canonical columns", {
  set.seed(1)
  x <- matrix(rnorm(200), 100, 2)
  y <- x[, 1] + rnorm(100)
  res <- hsic_test(x, y, n_permutations = 199, seed = 1)

  td <- tidy(res)
  expect_s3_class(td, "data.frame")
  expect_identical(nrow(td), 1L)
  expect_true(all(c("term", "statistic", "p.value", "n",
                    "n_permutations", "ess") %in% names(td)))
  # The broom-canonical dotted name carries the same value as the native field.
  expect_identical(td$p.value, res$p_value)
  expect_identical(td$statistic, res$statistic)
  expect_identical(td$term, res$method)
})

test_that("tidy() dispatches through broom and generics generics", {
  set.seed(2)
  a <- matrix(rnorm(100), 50, 2)
  b <- matrix(rnorm(100, mean = 0.5), 50, 2)
  res <- mmd_test(a, b, n_permutations = 199, seed = 2)

  td_generics <- generics::tidy(res)
  expect_s3_class(td_generics, "data.frame")
  expect_identical(td_generics$p.value, res$p_value)
})

test_that("tidy() on an mmd_ppc carries the PPC-specific columns", {
  set.seed(3)
  post <- matrix(rnorm(400), ncol = 2)
  obs  <- matrix(rnorm(40), ncol = 2)
  res <- mmd_ppc(post, obs, n_permutations = 199, seed = 3)

  td <- tidy(res)
  expect_true(all(c("surprise_bits", "reject") %in% names(td)))
  expect_identical(td$reject, res$reject)
  expect_equal(td$surprise_bits, res$surprise_bits)
})

test_that("tidy() on a taci_result threads decision and grounding", {
  set.seed(1)
  n <- 60
  nrate <- runif(n, 0, 200)
  yield <- 1.1 + 4.2 * (1 - exp(-0.018 * nrate)) + rnorm(n, 0, 0.25)
  post <- cbind(ymax = rnorm(150, 4.2, 0.30),
                rate = rnorm(150, 0.018, 0.004),
                y0   = rnorm(150, 1.1, 0.15))
  mitscherlich <- function(theta, X, t) {
    theta[3] + theta[1] * (1 - exp(-theta[2] * t))
  }

  # No provenance declared -> grounding is "[unverified]".
  res_ung <- taci_test(post, mitscherlich, X = matrix(1, n, 1),
                       treatment = nrate, outcome = yield,
                       n_perm = 80, seed = 1)
  td_ung <- tidy(res_ung)
  expect_identical(td_ung$term, "taci")
  expect_identical(td_ung$grounding, "[unverified]")
  expect_identical(td_ung$decision, res_ung$decision)
  expect_identical(td_ung$p.value, res_ung$p_h0)

  # Declaring provenance grounds the verdict.
  res_grd <- taci_test(post, mitscherlich, X = matrix(1, n, 1),
                       treatment = nrate, outcome = yield,
                       n_perm = 80, seed = 1,
                       mechanism_provenance = list(source = "test fixture"))
  expect_identical(tidy(res_grd)$grounding, "grounded")
})
