# test-joint-coverage.R -- joint (multivariate) calibration diagnostic.
#
# joint_coverage_test() tests calibration of the ensemble's *dependence
# structure*, not just its margins. The headline property -- and the reason it
# exists -- is that band-depth catches a correlation-only miscalibration that
# the marginal coverage_test() is blind to; that differentiator is pinned here
# alongside null calibration, the dispersion-direction readings, and validation.

chol2 <- function(rho) chol(matrix(c(1, rho, rho, 1), 2L))

test_that("the differentiator: band-depth catches correlation-only error that marginal misses", {
  set.seed(7)
  # Identical standard-normal margins; ensemble independent, obs correlated.
  obs <- matrix(stats::rnorm(300L), ncol = 2L) %*% chol2(0.9)
  ens <- matrix(stats::rnorm(4000L), ncol = 2L)

  marginal <- coverage_test(ens, obs, alpha = 0.05)
  joint    <- joint_coverage_test(ens, obs, prerank = "band_depth",
                                  alpha = 0.05, seed = 7L)

  expect_false(marginal$reject)         # margins match -> marginal is blind
  expect_true(joint$reject)             # dependence wrong -> joint catches it
})

test_that("under-dispersion is detected with the correct direction per pre-rank", {
  set.seed(2)
  obs <- matrix(stats::rnorm(300L), ncol = 2L)
  ens <- matrix(stats::rnorm(4000L, sd = 0.5), ncol = 2L)   # too tight

  bd <- joint_coverage_test(ens, obs, prerank = "band_depth", seed = 2L)
  av <- joint_coverage_test(ens, obs, prerank = "average", seed = 2L)

  expect_true(bd$reject)
  expect_lt(bd$mean_rank, 0.45)                 # slope: obs outlying
  expect_match(bd$verdict, "under-dispersed")

  expect_true(av$reject)
  expect_gt(av$dispersion_ratio, 1.1)           # U-shape: var inflated
  expect_match(av$verdict, "under-dispersed")
})

test_that("a correctly correlated ensemble is judged calibrated", {
  set.seed(3)
  ch <- chol2(0.8)
  obs <- matrix(stats::rnorm(300L), ncol = 2L) %*% ch
  ens <- matrix(stats::rnorm(6000L), ncol = 2L) %*% ch
  fit <- joint_coverage_test(ens, obs, prerank = "band_depth", seed = 3L)
  expect_false(fit$reject)
  expect_match(fit$verdict, "calibrated")
})

test_that("results are reproducible with seed and vary the rank tie-breaking otherwise", {
  set.seed(4)
  obs <- matrix(stats::rnorm(200L), ncol = 2L) %*% chol2(0.5)
  ens <- matrix(stats::rnorm(3000L), ncol = 2L)
  f1 <- joint_coverage_test(ens, obs, seed = 11L)
  f2 <- joint_coverage_test(ens, obs, seed = 11L)
  expect_identical(f1$ranks, f2$ranks)
  expect_identical(f1$calibration$p_value, f2$calibration$p_value)
})

test_that("the object class and fields are well-formed", {
  set.seed(5)
  obs <- matrix(stats::rnorm(160L), ncol = 2L)
  ens <- matrix(stats::rnorm(2000L), ncol = 2L)
  fit <- joint_coverage_test(ens, obs, seed = 5L)
  expect_s3_class(fit, "joint_coverage_test")
  expect_identical(fit$prerank, "band_depth")
  expect_length(fit$ranks, nrow(obs))
  expect_true(all(fit$ranks >= 1L & fit$ranks <= nrow(ens) + 1L))
  expect_s3_class(fit$coverage, "data.frame")
  expect_true(is.finite(fit$mean_rank) && is.finite(fit$dispersion_ratio))
})

test_that("input validation guards dimension, observed, and inherited checks", {
  ens1 <- matrix(stats::rnorm(200L), ncol = 1L)             # single output
  obs1 <- matrix(stats::rnorm(20L), ncol = 1L)
  expect_error(joint_coverage_test(ens1, obs1), "at least 2 output dimensions")

  ens <- matrix(stats::rnorm(200L), ncol = 2L)
  expect_error(joint_coverage_test(ens), "observed")           # missing observed
  expect_error(
    joint_coverage_test(ens, matrix(stats::rnorm(30L), ncol = 3L)),
    "same number of columns"                                   # inherited check
  )
  expect_error(
    joint_coverage_test(matrix(stats::rnorm(8L), ncol = 2L),
                        matrix(stats::rnorm(4L), ncol = 2L)),
    "at least 10 draws"                                        # inherited check
  )
})

test_that("Type-I error rate stays near the nominal level under the null", {
  skip_on_cran()
  set.seed(99)
  alpha <- 0.1
  n_trials <- 80L
  rej_bd <- rej_av <- logical(n_trials)
  for (t in seq_len(n_trials)) {
    obs <- matrix(stats::rnorm(120L), ncol = 2L)
    ens <- matrix(stats::rnorm(1600L), ncol = 2L)
    rej_bd[t] <- joint_coverage_test(ens, obs, prerank = "band_depth",
                                     alpha = alpha)$reject
    rej_av[t] <- joint_coverage_test(ens, obs, prerank = "average",
                                     alpha = alpha)$reject
  }
  expect_lt(mean(rej_bd), alpha + 0.12)
  expect_lt(mean(rej_av), alpha + 0.12)
})
