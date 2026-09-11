# Calibration suite: does a prediction get the LEVEL right, not just the order.
#
# The oracle is construction: data are generated with a known calibration
# slope, so the estimate is graded against a value this package did not
# compute.

# A calibrated prediction is one the outcome varies AROUND -- so the prediction
# is generated first and the outcome drawn from it. `truth + noise` is the
# opposite: it carries noise the outcome does not, so its slope attenuates to
# Var(truth)/(Var(truth) + Var(noise)) by construction. Getting that backwards
# is the classic regression-dilution trap and this suite detects it, which is
# the point.
test_that("a well-calibrated prediction recovers slope 1 and intercept 0", {
  set.seed(20260912)
  pred <- rnorm(400L, mean = 5, sd = 2)
  truth <- pred + rnorm(400L, sd = 0.3)
  r <- calibration_suite(pred, truth, n_boot = 300L, seed = 1L)
  expect_lt(abs(r$slope[["estimate"]] - 1), 0.05)
  expect_lt(abs(r$citl[["estimate"]]), 0.1)
  expect_true(r$slope[["lower"]] < 1 && r$slope[["upper"]] > 1)
  expect_identical(r$verdict, "calibrated")
  expect_false(r$abstained)
})

test_that("the slope recovers a KNOWN miscalibration", {
  # Predictions spread 1.6x too wide about the mean. Regressing observed on
  # predicted must return a slope near 1/1.6 = 0.625 -- an analytic target
  # this package plays no part in setting.
  set.seed(20260912)
  truth <- rnorm(600L, mean = 5, sd = 2)
  pred <- 5 + 1.6 * (truth - 5)
  r <- calibration_suite(pred, truth, n_boot = 300L, seed = 1L)
  expect_lt(abs(r$slope[["estimate"]] - 1 / 1.6), 0.02)
  expect_identical(r$verdict, "miscalibrated")
})

test_that("calibration-in-the-large isolates a pure level shift", {
  set.seed(20260912)
  truth <- rnorm(400L, mean = 5, sd = 2)
  pred <- truth + 0.75                       # level shift, slope untouched
  r <- calibration_suite(pred, truth, n_boot = 300L, seed = 1L)
  expect_lt(abs(r$citl[["estimate"]] + 0.75), 0.05)
  expect_lt(abs(r$slope[["estimate"]] - 1), 0.05)
})

# --- the headline property ---------------------------------------------------
# A calibration test on a handful of points cannot resolve a slope of 0.8 from
# a slope of 1.0. Reporting "calibrated" there reports the absence of power,
# not the presence of calibration. It must decline instead, and the decline
# must be typed so a caller gating on calibration can route on it.

# Power depends on signal-to-noise, not on n alone: six points on a nearly
# noiseless line pin the slope tightly. The regime that matters is a small,
# NOISY sample -- which is exactly a held-out site-year.
test_that("a small noisy sample yields 'undetermined', not a pass", {
  set.seed(20260912)
  pred <- rnorm(6L, mean = 5, sd = 2)
  truth <- pred + rnorm(6L, sd = 3)
  r <- calibration_suite(pred, truth, n_boot = 300L, tolerance = 0.1, seed = 1L)
  expect_identical(r$verdict, "undetermined")
  expect_true(r$abstained)
  expect_s3_class(r, "kernR_abstention")
  expect_true(is_kernR_abstention(r))
  # and it says what it could have detected
  expect_gt(r$detectable, r$tolerance)
})

test_that("the same data certifies once the tolerance is honest about n", {
  set.seed(20260912)
  pred <- rnorm(6L, mean = 5, sd = 2)
  truth <- pred + rnorm(6L, sd = 3)
  strict <- calibration_suite(pred, truth, n_boot = 300L, tolerance = 0.05,
                              seed = 1L)
  loose <- calibration_suite(pred, truth, n_boot = 300L, tolerance = 5,
                             seed = 1L)
  expect_identical(strict$verdict, "undetermined")
  expect_identical(loose$verdict, "calibrated")
  # The estimates are identical; only the claim the caller may make changes.
  expect_equal(strict$slope[["estimate"]], loose$slope[["estimate"]])
})

test_that("a large sample can detect what a small one cannot", {
  set.seed(20260912)
  pred_big <- rnorm(2000L, mean = 5, sd = 2)
  big <- calibration_suite(pred_big, pred_big + rnorm(2000L, sd = 1),
                           n_boot = 300L, seed = 1L)
  pred_small <- rnorm(20L, mean = 5, sd = 2)
  small <- calibration_suite(pred_small, pred_small + rnorm(20L, sd = 1),
                             n_boot = 300L, seed = 1L)
  expect_lt(big$detectable, small$detectable)
  expect_identical(big$verdict, "calibrated")
})

test_that("the seed makes the bootstrap interval reproducible", {
  set.seed(20260912)
  pred <- rnorm(200L, mean = 5, sd = 2)
  truth <- pred + rnorm(200L, sd = 0.4)
  a <- calibration_suite(pred, truth, n_boot = 200L, seed = 99L)
  b <- calibration_suite(pred, truth, n_boot = 200L, seed = 99L)
  expect_equal(a$slope, b$slope)
  expect_equal(a$ici, b$ici)
})

test_that("input is validated rather than coerced into a wrong answer", {
  expect_error(calibration_suite(1:5, 1:4), "same length")
  expect_error(calibration_suite(numeric(0), numeric(0)), "non-empty")
  expect_error(calibration_suite(c(1, NA, 3), c(1, 2, 3)), "finite")
  expect_error(calibration_suite(1:10, 1:10, tolerance = 0), "tolerance")
  expect_error(calibration_suite(1:10, 1:10, conf_level = 1), "conf_level")
})

test_that("a constant prediction has no slope and does not fabricate one", {
  r <- calibration_suite(rep(3, 30L), rnorm(30L), n_boot = 50L, seed = 1L)
  expect_true(is.na(r$slope[["estimate"]]))
  expect_identical(r$verdict, "undetermined")
  expect_true(is_kernR_abstention(r))
})

test_that("regression dilution is DETECTED, not waved through", {
  # `truth + noise` is over-dispersed as a prediction of `truth`: the slope
  # attenuates to Var(truth)/(Var(truth) + Var(noise)). With sd 2 and 0.6 that
  # is 4/(4 + 0.36) = 0.917, an analytic target this package does not set.
  set.seed(20260912)
  truth <- rnorm(3000L, mean = 5, sd = 2)
  pred <- truth + rnorm(3000L, sd = 0.6)
  r <- calibration_suite(pred, truth, n_boot = 300L, seed = 1L)
  expect_lt(abs(r$slope[["estimate"]] - 4 / 4.36), 0.03)
  expect_identical(r$verdict, "miscalibrated")
})
