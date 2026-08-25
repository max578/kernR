# Tests for the typed kernR_abstention marker (K4 residual / K8 close-out).
#
# The leader-node `is_orchestra_decline()` predicate
# (`ORCHESTRA_dev/integration/refusal_contract.R`) recognises either the
# shared `orchestra_refusal` marker or any class matching `_(refusal|
# abstention)$`. This file exercises that same naming convention directly
# (`grepl("_(refusal|abstention)$", class(x))`) rather than sourcing the
# contract file itself, so the package suite stays self-contained; the two
# are kept in lock-step by inspection, not by a shared path.
#
# Each of the three reliability gates named in the audit (K4) is exercised
# through the real verb, not by hand-building a result list, so a
# regression in the gate itself would also break these tests.

.is_orchestra_decline_locally <- function(x) {
  any(grepl("_(refusal|abstention)$", class(x)))
}

test_that("bd_hsic_test() stamps kernR_abstention when the ESS floor trips", {
  set.seed(42L)
  n <- 60L
  z <- matrix(stats::rnorm(n * 2L), n, 2L)
  x <- z[, 1L] + stats::rnorm(n, sd = 0.5)
  y <- 0.5 * x + stats::rnorm(n, sd = 0.5)

  res <- suppressWarnings(bd_hsic_test(
    x, y, z,
    density_ratio    = "logistic",
    n_permutations   = 49L,
    seed             = 42L,
    min_ess_fraction = 0.999
  ))

  expect_true(res$ess_warning)
  expect_s3_class(res, "kernel_test_result")
  expect_s3_class(res, "kernR_abstention")
  expect_true(.is_orchestra_decline_locally(res))
})

test_that("bd_hsic_test() stamps kernR_abstention when the density-ratio fit gate trips", {
  # The proxymix fit-quality gate is DGP-dependent to trigger on real data
  # (K7); mock the two proxymix-backend calls so the C4 fit-quality flag
  # fires deterministically, exercising bd_hsic_test()'s own gate logic
  # rather than proxymix's convergence behaviour.
  local_mocked_bindings(
    fit_density_ratio = function(...) structure(
      list(
        method      = "proxymix",
        fit_quality = list(ok = FALSE, status = "degraded",
                           reason = "mocked non-convergence for testing",
                           marg_converged = FALSE),
        ncol_x      = 1L, ncol_z = 2L
      ),
      class = "density_ratio_fit"
    ),
    predict_density_ratio = function(object, new_x, new_z, type) {
      # Non-degenerate weights so the downstream propensity clustering
      # step is well posed; the specific values are immaterial to what
      # this test asserts.
      stats::runif(nrow(new_x), 0.5, 1.5)
    }
  )

  set.seed(1L)
  n <- 20L
  z <- matrix(stats::rnorm(n * 2L), n, 2L)
  x <- z[, 1L] + stats::rnorm(n, sd = 0.5)
  y <- 0.5 * x + stats::rnorm(n, sd = 0.5)

  res <- suppressWarnings(bd_hsic_test(
    x, y, z,
    density_ratio    = "proxymix",
    n_permutations   = 19L,
    seed             = 1L,
    min_ess_fraction = 0
  ))

  expect_true(res$density_ratio_warning)
  expect_s3_class(res, "kernel_test_result")
  expect_s3_class(res, "kernR_abstention")
  expect_true(.is_orchestra_decline_locally(res))
})

test_that("taci_test() stamps kernR_abstention when the posterior-adequacy guard trips", {
  mechanism <- function(theta, X, t) theta[3] + theta[1] * (1 - exp(-theta[2] * t))
  set.seed(1L)
  n <- 60L
  nrate <- stats::runif(n, 0, 200)
  yield <- 1.1 + 4.2 * (1 - exp(-0.018 * nrate)) + stats::rnorm(n, 0, 0.25)
  # Near-degenerate posterior (tiny sd on every parameter): the model-implied
  # effect is over-determined, so effect_cv < 0.02 and the adequacy guard
  # trips.
  post <- cbind(
    ymax = stats::rnorm(100L, 4.2, 0.0001),
    rate = stats::rnorm(100L, 0.018, 0.00001),
    y0   = stats::rnorm(100L, 1.1, 0.0001)
  )

  verdict <- suppressWarnings(taci_test(
    post, mechanism, X = matrix(1, n, 1L),
    treatment = nrate, outcome = yield,
    n_perm = 49L, seed = 1L
  ))

  expect_false(verdict$posterior_adequacy$ok)
  expect_s3_class(verdict, "taci_result")
  expect_s3_class(verdict, "kernR_abstention")
  expect_true(.is_orchestra_decline_locally(verdict))
})

test_that("kernR_abstention is not stamped when every reliability gate passes", {
  set.seed(7L)
  n <- 200L
  z <- matrix(stats::rnorm(n * 2L), n, 2L)
  x <- z[, 1L] + stats::rnorm(n)
  y <- 0.7 * x + z[, 2L] + stats::rnorm(n, sd = 0.4)

  res <- bd_hsic_test(x, y, z, density_ratio = "logistic",
                      n_permutations = 49L, seed = 7L)
  expect_false(res$ess_warning)
  expect_false(res$density_ratio_warning)
  expect_false(inherits(res, "kernR_abstention"))
})

test_that("as_orchestra_manifest() sets summary$abstained from the typed marker, not the raw flag directly", {
  set.seed(42L)
  n <- 60L
  z <- matrix(stats::rnorm(n * 2L), n, 2L)
  x <- z[, 1L] + stats::rnorm(n, sd = 0.5)
  y <- 0.5 * x + stats::rnorm(n, sd = 0.5)

  res <- suppressWarnings(bd_hsic_test(
    x, y, z,
    density_ratio    = "logistic",
    n_permutations   = 49L,
    seed             = 42L,
    min_ess_fraction = 0.999
  ))
  expect_true(inherits(res, "kernR_abstention"))

  m <- as_orchestra_manifest(res)
  expect_true(m@summary$abstained)
})
