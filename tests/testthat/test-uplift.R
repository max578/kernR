# Regression guards for the 0.3.0 correctness uplift:
#  * permutation reproducibility via `seed=` (C++ RNG now routed through R)
#  * DR-DATE / DR-DETT genuine augmented-IPW (outcome model is actually used)
#  * per-arm effective-sample-size reliability gate

test_that("HSIC permutation test is reproducible via seed", {
  set.seed(99)
  x <- matrix(rnorm(120), 60, 2)
  y <- x[, 1] + rnorm(60)
  r1 <- hsic_test(x, y, n_permutations = 200, seed = 7)
  r2 <- hsic_test(x, y, n_permutations = 200, seed = 7)
  expect_identical(r1$null_distribution, r2$null_distribution)
  expect_identical(r1$p_value, r2$p_value)

  r3 <- hsic_test(x, y, n_permutations = 200, seed = 8)
  expect_false(identical(r1$null_distribution, r3$null_distribution))
})

test_that("MMD permutation test is reproducible via seed", {
  set.seed(99)
  a <- matrix(rnorm(100), 50, 2)
  b <- matrix(rnorm(100, mean = 0.5), 50, 2)
  m1 <- mmd_test(a, b, n_permutations = 200, seed = 3)
  m2 <- mmd_test(a, b, n_permutations = 200, seed = 3)
  expect_identical(m1$null_distribution, m2$null_distribution)
})

test_that("DR-DATE is reproducible via seed", {
  set.seed(5)
  n <- 150
  x <- matrix(rnorm(n * 2), n, 2)
  t <- rbinom(n, 1, plogis(0.4 * x[, 1]))
  y <- t + x[, 1] + rnorm(n)
  d1 <- dr_date_test(y, t, x, n_permutations = 100, seed = 11)
  d2 <- dr_date_test(y, t, x, n_permutations = 100, seed = 11)
  expect_identical(d1$statistic, d2$statistic)
  expect_identical(d1$p_value, d2$p_value)
})

test_that("DR-DATE outcome model changes the statistic (genuine AIPW)", {
  set.seed(7)
  n <- 200
  x <- matrix(rnorm(n * 2), n, 2)
  t <- rbinom(n, 1, plogis(0.6 * x[, 1]))
  # strong covariate -> outcome dependence so the CME augmentation matters
  y <- t * 0.8 + 1.5 * x[, 1] + 0.8 * x[, 2] + rnorm(n, sd = 0.4)
  s_krr <- dr_date_test(y, t, x,
    outcome_model = "krr",
    n_permutations = 0, seed = 1
  )$statistic
  s_ipw <- dr_date_test(y, t, x,
    outcome_model = "zero",
    n_permutations = 0, seed = 1
  )$statistic
  # Pre-0.3.0 these were byte-identical (outcome model was discarded).
  expect_false(isTRUE(all.equal(s_krr, s_ipw)))
})

test_that("DR-DATE reports per-arm ESS and warns on overlap collapse", {
  set.seed(3)
  n <- 200
  x <- matrix(rnorm(n * 2), n, 2)
  res <- dr_date_test(y = x[, 1] + rnorm(n),
    treatment = rbinom(n, 1, plogis(0.3 * x[, 1])),
    covariates = x, outcome_model = "zero",
    n_permutations = 30, seed = 1
  )
  expect_true(is.finite(res$ess))
  expect_false(is.na(res$ess))

  # Near-deterministic assignment -> extreme weights -> tiny ESS -> warning.
  # Muffle the incidental glm separation warning so only the ESS reliability
  # warning is asserted.
  set.seed(4)
  xb <- matrix(rnorm(n * 2), n, 2)
  tb <- rbinom(n, 1, plogis(10 * xb[, 1]))
  yb <- tb + xb[, 1] + rnorm(n)
  expect_warning(
    withCallingHandlers(
      dr_date_test(yb, tb, xb,
        outcome_model = "zero",
        n_permutations = 30, seed = 1
      ),
      warning = function(w) {
        if (grepl("fitted probabilities", conditionMessage(w))) {
          invokeRestart("muffleWarning")
        }
      }
    ),
    "[Ee]ffective sample size"
  )
})

test_that("DR-DETT is reproducible and uses the control outcome model", {
  set.seed(8)
  n <- 200
  x <- matrix(rnorm(n * 2), n, 2)
  t <- rbinom(n, 1, plogis(0.5 * x[, 1]))
  y <- t * rnorm(n, sd = 1.5) + 1.2 * x[, 1] + rnorm(n, sd = 0.5)
  r1 <- dr_dett_test(y, t, x, n_permutations = 0, seed = 2)$statistic
  r2 <- dr_dett_test(y, t, x, n_permutations = 0, seed = 2)$statistic
  expect_identical(r1, r2)
  s_zero <- dr_dett_test(y, t, x,
    outcome_model = "zero",
    n_permutations = 0, seed = 2
  )$statistic
  expect_false(isTRUE(all.equal(r1, s_zero)))
})
