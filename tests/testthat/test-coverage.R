# Tests for coverage_test(): the PIT-based coverage / calibration diagnostic --
# calibrated vs under-/over-dispersed vs biased detection, coverage accuracy,
# the pesto_ensemble / manifest methods, and input validation.

test_that("a calibrated ensemble is not rejected and coverage matches nominal", {
  set.seed(1)
  obs <- matrix(stats::rnorm(200L), ncol = 2L)
  fit <- coverage_test(matrix(stats::rnorm(4000L), ncol = 2L), obs)

  expect_s3_class(fit, "coverage_test")
  expect_gt(fit$calibration$p_value, 0.05)
  expect_false(fit$reject)
  expect_lt(abs(fit$dispersion_ratio - 1), 0.25)          # ratio near 1
  cov90 <- fit$coverage$empirical[fit$coverage$nominal == 0.9]
  expect_gt(cov90, 0.83)                                   # ~90% covered
  expect_match(fit$verdict, "calibrated")
})

test_that("an under-dispersed (over-confident) ensemble is detected", {
  set.seed(2)
  obs <- matrix(stats::rnorm(200L), ncol = 2L)
  fit <- coverage_test(matrix(stats::rnorm(4000L, sd = 0.4), ncol = 2L), obs)

  expect_true(fit$reject)
  expect_gt(fit$dispersion_ratio, 1.1)
  cov90 <- fit$coverage$empirical[fit$coverage$nominal == 0.9]
  expect_lt(cov90, 0.75)                                   # 90% intervals miss
  expect_match(fit$verdict, "under-dispersed")
})

test_that("an over-dispersed ensemble is detected", {
  set.seed(3)
  obs <- matrix(stats::rnorm(200L), ncol = 2L)
  fit <- coverage_test(matrix(stats::rnorm(4000L, sd = 2.5), ncol = 2L), obs)

  expect_true(fit$reject)
  expect_lt(fit$dispersion_ratio, 0.9)
  expect_match(fit$verdict, "over-dispersed")
})

test_that("a biased ensemble is detected with the correct direction", {
  set.seed(4)
  obs <- matrix(stats::rnorm(200L), ncol = 2L)
  # Ensemble shifted high relative to observations
  fit <- coverage_test(matrix(stats::rnorm(4000L) + 1.2, ncol = 2L), obs)

  expect_true(fit$reject)
  expect_lt(fit$mean_pit, 0.4)                             # obs below ensemble
  expect_match(fit$verdict, "too high")
})

test_that("PIT values are in [0, 1] with the right shape", {
  set.seed(5)
  obs <- matrix(stats::rnorm(60L), ncol = 3L)
  fit <- coverage_test(matrix(stats::rnorm(900L), ncol = 3L), obs)
  expect_equal(dim(fit$pit), c(20L, 3L))
  expect_true(all(fit$pit >= 0 & fit$pit <= 1))
  expect_equal(nrow(fit$coverage_by_dim), length(fit$levels))
})

test_that("the print method runs and is reproducible", {
  set.seed(6)
  obs <- matrix(stats::rnorm(100L), ncol = 2L)
  ens <- matrix(stats::rnorm(2000L), ncol = 2L)
  a <- coverage_test(ens, obs)
  b <- coverage_test(ens, obs)
  expect_identical(a$dispersion_ratio, b$dispersion_ratio)   # deterministic
  expect_output(print(a), "Coverage / calibration diagnostic")
  expect_output(print(a), "Dispersion")
})

# pesto_ensemble / manifest methods ------------------------------------------

test_that("the pesto_ensemble method uses the bundled observations", {
  set.seed(7)
  post <- matrix(stats::rnorm(4000L), ncol = 2L)
  obs  <- matrix(stats::rnorm(200L), ncol = 2L)
  ens  <- pesto_ensemble(post, obs, metadata = list(run = "r1"))
  fit  <- coverage_test(ens)
  expect_s3_class(fit, "coverage_test")
  # Robust to the uniformity test's Type-I rate: a calibrated ensemble has a
  # dispersion ratio near 1 regardless of the occasional rank-histogram blip.
  expect_lt(abs(fit$dispersion_ratio - 1), 0.3)
  expect_equal(fit$pesto_metadata$run, "r1")
})

test_that("the manifest method threads fidelity provenance", {
  skip_if_not_installed("PESTO")
  set.seed(8)
  n <- 300L
  outputs <- data.frame(
    real_name = paste0("r", seq_len(n)),
    o1 = stats::rnorm(n), o2 = stats::rnorm(n)
  )
  params <- data.frame(real_name = paste0("r", seq_len(n)),
                       p1 = stats::rnorm(n))
  m <- PESTO::pesto_ensemble_manifest(
    run_id = "cov", params = params, outputs = outputs,
    weights = c(o1 = 1, o2 = 1), obs_target = c(o1 = 0, o2 = 0),
    data_hash = "sha256:test", pesto_version = "0.4.1",
    timestamp = Sys.time(), method = "ies_callback",
    noptmax = 3L, lambda_schedule = c(1, 1, 1),
    fidelity = list(type = "multifidelity", final_level = 1L, n_levels = 2L)
  )
  obs <- matrix(stats::rnorm(40L), ncol = 2L)
  fit <- coverage_test(m, observed = obs)
  expect_s3_class(fit, "coverage_test")
  expect_equal(fit$pesto_metadata$run_id, "cov")
  expect_false(is.null(fit$pesto_metadata$fidelity))
})

# Input validation -----------------------------------------------------------

test_that("coverage_test() rejects malformed input", {
  obs <- matrix(stats::rnorm(40L), ncol = 2L)
  good <- matrix(stats::rnorm(2000L), ncol = 2L)

  expect_error(coverage_test(good), "observed")
  expect_error(coverage_test(matrix(stats::rnorm(18L), ncol = 2L), obs),
               "at least 10 draws")
  expect_error(coverage_test(good, matrix(stats::rnorm(30L), ncol = 3L)),
               "same number of columns")
  bad <- good
  bad[1L, 1L] <- NA_real_
  expect_error(coverage_test(bad, obs), "finite")
  expect_error(coverage_test(good, obs, levels = c(0.5, 1.2)), "in \\(0, 1\\)")
  expect_error(coverage_test(good, obs, n_bins = 1L), ">= 2")
  expect_error(coverage_test(good, obs, alpha = 0), "in \\(0, 1\\)")
})
