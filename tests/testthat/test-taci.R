.taci_mitscherlich <- function(theta, X, t) {
  theta[3] + theta[1] * (1 - exp(-theta[2] * t))
}

test_that("TACI recovers a mechanism-consistent continuous effect", {
  set.seed(1)
  n <- 120
  nrate <- runif(n, 0, 200)
  yield <- 1.1 + 4.2 * (1 - exp(-0.018 * nrate)) + rnorm(n, 0, 0.25)
  post <- cbind(ymax = rnorm(200, 4.2, 0.30),
                rate = rnorm(200, 0.018, 0.004),
                y0   = rnorm(200, 1.1, 0.15))

  res <- taci_test(post, .taci_mitscherlich, X = matrix(1, n, 1),
                   treatment = nrate, outcome = yield,
                   n_perm = 150, seed = 1)

  expect_s3_class(res, "taci_result")
  expect_identical(res$treatment_type, "continuous")
  expect_true(res$in_tail)
  expect_identical(res$decision, "mechanism_consistent_effect")
  expect_true(res$posterior_adequacy$ok)
})

test_that("TACI returns no_effect on a flat outcome", {
  set.seed(2)
  n <- 120
  nrate <- runif(n, 0, 200)
  yield <- 2.5 + rnorm(n, 0, 0.30)            # no N effect in the data
  post <- cbind(ymax = rnorm(200, 4.2, 0.30),
                rate = rnorm(200, 0.018, 0.004),
                y0   = rnorm(200, 1.1, 0.15))

  res <- taci_test(post, .taci_mitscherlich, X = matrix(1, n, 1),
                   treatment = nrate, outcome = yield,
                   n_perm = 150, seed = 1)

  expect_false(res$in_tail)
  expect_identical(res$decision, "no_effect")
})

test_that("the posterior-adequacy guard flags an over-determined effect", {
  set.seed(3)
  n <- 120
  nrate <- runif(n, 0, 200)
  yield <- 1.1 + 4.2 * (1 - exp(-0.018 * nrate)) + rnorm(n, 0, 0.25)
  # A near-degenerate posterior pins the model-implied effect too precisely.
  post <- cbind(ymax = rnorm(200, 4.2, 1e-4),
                rate = rnorm(200, 0.018, 1e-6),
                y0   = rnorm(200, 1.1, 1e-4))

  expect_warning(
    res <- taci_test(post, .taci_mitscherlich, X = matrix(1, n, 1),
                     treatment = nrate, outcome = yield,
                     n_perm = 100, seed = 1),
    "over-determined"
  )
  expect_false(res$posterior_adequacy$ok)
  expect_true(res$posterior_adequacy$effect_cv < 0.02)
})

test_that("TACI handles a binary treatment", {
  set.seed(4)
  n <- 160
  trt <- rbinom(n, 1, 0.5)
  y <- 1.0 + 1.5 * trt + rnorm(n, 0, 0.5)
  # Linear mechanism: theta = (effect, intercept).
  mech <- function(theta, X, t) theta[2] + theta[1] * t
  post <- cbind(effect = rnorm(200, 1.5, 0.25),
                intercept = rnorm(200, 1.0, 0.2))

  res <- taci_test(post, mech, X = matrix(1, n, 1),
                   treatment = trt, outcome = y,
                   n_perm = 150, seed = 1)

  expect_s3_class(res, "taci_result")
  expect_identical(res$treatment_type, "binary")
  expect_true(res$in_tail)
})

test_that("TACI validates the posterior argument", {
  mech <- function(theta, X, t) theta[1] * t
  expect_error(
    taci_test(list(1, 2), mech, X = matrix(1, 4, 1),
              treatment = c(0, 1, 0, 1), outcome = rnorm(4)),
    "matrix or data.frame"
  )
})

test_that("print.taci_result returns its input invisibly", {
  set.seed(5)
  n <- 80
  nrate <- runif(n, 0, 200)
  yield <- 1.1 + 4.2 * (1 - exp(-0.018 * nrate)) + rnorm(n, 0, 0.25)
  post <- cbind(ymax = rnorm(150, 4.2, 0.30),
                rate = rnorm(150, 0.018, 0.004),
                y0   = rnorm(150, 1.1, 0.15))
  res <- taci_test(post, .taci_mitscherlich, X = matrix(1, n, 1),
                   treatment = nrate, outcome = yield,
                   n_perm = 80, seed = 1)

  expect_invisible(print(res))
})

# --- FX-9: mechanism-provenance grounding label (Independent Oracle Principle) -

test_that("TACI labels an un-grounded verdict [unverified] by default", {
  set.seed(1)
  n <- 120
  nrate <- runif(n, 0, 200)
  yield <- 1.1 + 4.2 * (1 - exp(-0.018 * nrate)) + rnorm(n, 0, 0.25)
  post <- cbind(ymax = rnorm(200, 4.2, 0.30),
                rate = rnorm(200, 0.018, 0.004),
                y0   = rnorm(200, 1.1, 0.15))

  # No mechanism_provenance: the verdict must be labelled, not confident.
  res <- taci_test(post, .taci_mitscherlich, X = matrix(1, n, 1),
                   treatment = nrate, outcome = yield, n_perm = 150, seed = 1)
  expect_identical(res$decision, "mechanism_consistent_effect")  # enum intact
  expect_identical(res$grounding, "[unverified]")
  expect_identical(res$verdict, "mechanism_consistent_effect [unverified]")
  expect_true(is.na(res$mechanism_provenance))

  # Declaring provenance grounds the verdict; the verdict string drops the tag.
  res2 <- taci_test(post, .taci_mitscherlich, X = matrix(1, n, 1),
                    treatment = nrate, outcome = yield, n_perm = 150, seed = 1,
                    mechanism_provenance = list(run_id = "pesto-123",
                                                apsim_version = "2026.5.8046.0"))
  expect_identical(res2$grounding, "grounded")
  expect_identical(res2$verdict, "mechanism_consistent_effect")
  expect_false(isTRUE(is.na(res2$mechanism_provenance)))
  # Grounding is pure metadata: the scientific decision is unchanged either way.
  expect_identical(res$decision, res2$decision)
})
