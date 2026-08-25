# test-manifest.R -- as_orchestra_manifest() emits kernR's verdicts as
# orchestra_manifest S7 objects (K8): a taci_result or a kernel_test_result
# can now enter decideR through the C2-symmetric emit contract.

test_that("a TACI verdict emits a treatment_effects manifest, reliable case", {
  set.seed(1L)
  n <- 60L
  X <- matrix(stats::rnorm(n), ncol = 1L)
  Tt <- stats::rbinom(n, 1L, 0.5)
  mech <- function(theta, X, t) theta[1] + theta[2] * t
  Y <- mech(c(1, 2), X, Tt) + stats::rnorm(n, sd = 0.3)
  Theta <- cbind(rnorm(200L, 1, 0.05), rnorm(200L, 2, 0.3))

  fit <- taci_test(posterior = Theta, mechanism = mech, X = X,
                   treatment = Tt, outcome = Y,
                   treatment_type = "binary", n_perm = 100L, seed = 1L)

  m <- as_orchestra_manifest(fit)
  expect_s3_class(m, "orchestra_manifest")
  expect_equal(m@inferential_target, "treatment_effects")
  expect_equal(m@emitter_package, "kernR")
  expect_s3_class(m@summary, "manifest_summary")
  expect_equal(m@summary$abstained, !fit$posterior_adequacy$ok)
  expect_equal(m@summary$metrics$decision, fit$decision)
  expect_true(verify_manifest(m)$ok)
})

test_that("an unreliable TACI verdict (posterior over-determined) abstains in the manifest", {
  set.seed(2L)
  n <- 40L
  X <- matrix(stats::rnorm(n), ncol = 1L)
  Tt <- stats::rbinom(n, 1L, 0.5)
  mech <- function(theta, X, t) theta[1] + theta[2] * t
  Y <- mech(c(1, 2), X, Tt) + stats::rnorm(n, sd = 0.3)
  # A near-degenerate posterior (tiny sd on the slope) trips the
  # posterior-adequacy guard (effect_cv < 0.02).
  Theta <- cbind(rnorm(200L, 1, 1e-6), rnorm(200L, 2, 1e-6))

  fit <- suppressWarnings(
    taci_test(posterior = Theta, mechanism = mech, X = X,
             treatment = Tt, outcome = Y,
             treatment_type = "binary", n_perm = 50L, seed = 2L)
  )
  expect_false(fit$posterior_adequacy$ok)

  m <- as_orchestra_manifest(fit)
  expect_true(m@summary$abstained)
  expect_false(is.na(m@summary$abstain_reason))
})

test_that("a kernel_test_result emits a treatment_effects manifest", {
  set.seed(3L)
  x <- matrix(stats::rnorm(100L), ncol = 1L)
  y <- matrix(stats::rnorm(100L), ncol = 1L)
  fit <- hsic_test(x, y, n_permutations = 100L)
  expect_s3_class(fit, "kernel_test_result")

  m <- as_orchestra_manifest(fit)
  expect_s3_class(m, "orchestra_manifest")
  expect_equal(m@inferential_target, "treatment_effects")
  expect_equal(m@method, "kernR:HSIC")
  expect_false(m@summary$abstained)
  expect_true(verify_manifest(m)$ok)
})

test_that("verify_manifest() detects a tampered payload", {
  set.seed(4L)
  x <- matrix(stats::rnorm(80L), ncol = 1L)
  y <- matrix(stats::rnorm(80L), ncol = 1L)
  fit <- hsic_test(x, y, n_permutations = 80L)
  m <- as_orchestra_manifest(fit)
  m@summary$metrics$statistic <- m@summary$metrics$statistic + 1
  expect_false(verify_manifest(m)$ok)
})
