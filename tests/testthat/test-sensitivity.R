test_that("hsic_sensitivity returns indices in [0, 1]", {
  set.seed(1)
  n <- 80L
  theta <- matrix(stats::runif(n * 3L), nrow = n,
                  dimnames = list(NULL, c("a", "b", "c")))
  y <- 2 * theta[, "a"] + theta[, "b"]^2 +
       stats::rnorm(n, sd = 0.1)
  fit <- hsic_sensitivity(theta, y, n_permutations = 99L, seed = 1L)
  expect_s3_class(fit, "hsic_sensitivity")
  expect_equal(dim(fit$index), c(3L, 1L))
  expect_true(all(fit$index >= 0 & fit$index <= 1))
  expect_true(all(fit$total_index >= 0 & fit$total_index <= 1))
})

test_that("hsic_sensitivity ranks active parameters above inert", {
  set.seed(2)
  n <- 100L
  theta <- matrix(stats::runif(n * 3L), nrow = n,
                  dimnames = list(NULL, c("active", "weak", "inert")))
  y <- 3 * theta[, "active"] + 0.5 * theta[, "weak"] +
       stats::rnorm(n, sd = 0.1)
  fit <- hsic_sensitivity(theta, y, n_permutations = 99L, seed = 2L)
  expect_equal(fit$param_names[fit$rank[1L]], "active")
  expect_equal(fit$param_names[fit$rank[3L]], "inert")
})

test_that("hsic_sensitivity is reproducible with seed", {
  set.seed(3)
  n <- 60L
  theta <- matrix(stats::runif(n * 2L), nrow = n,
                  dimnames = list(NULL, c("p1", "p2")))
  y <- cbind(o1 = theta[, 1L] + stats::rnorm(n, sd = 0.1),
             o2 = stats::rnorm(n))
  f1 <- hsic_sensitivity(theta, y, n_permutations = 99L, seed = 42L)
  f2 <- hsic_sensitivity(theta, y, n_permutations = 99L, seed = 42L)
  expect_identical(f1$index,    f2$index)
  expect_identical(f1$p_value,  f2$p_value)
})

test_that("hsic_sensitivity supports p_value = FALSE (skip null)", {
  set.seed(4)
  n <- 60L
  theta <- matrix(stats::runif(n * 2L), nrow = n,
                  dimnames = list(NULL, c("a", "b")))
  y <- theta[, 1L] + stats::rnorm(n, sd = 0.1)
  fit <- hsic_sensitivity(theta, y, p_value = FALSE,
                          n_permutations = 99L, seed = 4L)
  expect_null(fit$p_value)
  expect_null(fit$p_value_adjusted)
  expect_true(all(is.finite(fit$index)))
})

test_that("hsic_sensitivity validates inputs", {
  expect_error(
    hsic_sensitivity(matrix(0, 5L, 2L), matrix(0, 5L, 1L)),
    "At least 10"
  )
  expect_error(
    hsic_sensitivity(matrix(0, 20L, 2L), matrix(0, 10L, 1L)),
    "same number of rows"
  )
  expect_error(
    hsic_sensitivity(matrix(0, 20L, 2L), matrix(0, 20L, 1L),
                     n_permutations = 0L),
    "n_permutations"
  )
})

test_that("print, plot and as.data.frame work; as.data.frame respects p_value=FALSE", {
  set.seed(5)
  n <- 40L
  theta <- matrix(stats::runif(n * 2L), nrow = n,
                  dimnames = list(NULL, c("a", "b")))
  y <- cbind(o1 = theta[, 1L] + stats::rnorm(n, sd = 0.2),
             o2 = stats::rnorm(n))
  fit <- hsic_sensitivity(theta, y, n_permutations = 99L, seed = 5L)
  expect_output(print(fit), "HSIC-Sensitivity")
  expect_silent({
    grDevices::pdf(tempfile(fileext = ".pdf"))
    on.exit(grDevices::dev.off(), add = TRUE)
    plot(fit)
  })
  df <- as.data.frame(fit)
  expect_equal(nrow(df), 4L)
  expect_setequal(colnames(df),
                  c("parameter", "output", "index", "statistic",
                    "p_value", "p_value_adjusted"))

  fit_no_p <- hsic_sensitivity(theta, y, p_value = FALSE,
                               n_permutations = 99L, seed = 5L)
  df_no_p <- as.data.frame(fit_no_p)
  expect_setequal(colnames(df_no_p),
                  c("parameter", "output", "index", "statistic"))
})

test_that("HSIC-SI scales like a Sobol-comparable index on a near-additive case", {
  # On Y = a*X1 + b*X2 + small noise, X1 with larger coefficient should
  # have larger index. This is qualitative -- HSIC-SI is not Sobol but
  # should preserve ordering.
  set.seed(6)
  n <- 200L
  theta <- matrix(stats::runif(n * 2L), nrow = n,
                  dimnames = list(NULL, c("big", "small")))
  y <- 5 * theta[, "big"] + 1 * theta[, "small"] +
       stats::rnorm(n, sd = 0.05)
  fit <- hsic_sensitivity(theta, y, p_value = FALSE,
                          n_permutations = 50L, seed = 6L)
  expect_gt(fit$index["big", 1L], fit$index["small", 1L])
})

# ---- Total-order indices --------------------------------------------------

test_that("total_order = TRUE returns total-order matrices in [0, 1]", {
  set.seed(11)
  n <- 100L
  theta <- matrix(stats::runif(n * 3L), nrow = n,
                  dimnames = list(NULL, c("a", "b", "c")))
  y <- 2 * theta[, "a"] +
       theta[, "a"] * theta[, "b"] +
       stats::rnorm(n, sd = 0.1)
  fit <- hsic_sensitivity(theta, y, total_order = TRUE,
                          p_value = FALSE, n_permutations = 99L,
                          seed = 11L)
  expect_true(fit$total_order)
  expect_equal(dim(fit$index_total_order), c(3L, 1L))
  expect_true(all(fit$index_total_order >= 0 &
                  fit$index_total_order <= 1))
  expect_equal(dim(fit$statistic_total_order), c(3L, 1L))
})

test_that("total_order = TRUE requires p >= 2", {
  set.seed(12)
  n <- 60L
  theta <- matrix(stats::runif(n), ncol = 1L,
                  dimnames = list(NULL, "only_one"))
  y <- theta[, 1L] + stats::rnorm(n, sd = 0.1)
  expect_error(
    hsic_sensitivity(theta, y, total_order = TRUE,
                     n_permutations = 50L, seed = 12L),
    "p >= 2"
  )
})

test_that("total_order recovers ~first-order on purely additive models", {
  # Y = X1 + X2 with no interaction => T_j ~ S_j for both j.
  set.seed(13)
  n <- 250L
  theta <- matrix(stats::runif(n * 2L), nrow = n,
                  dimnames = list(NULL, c("a", "b")))
  y <- theta[, "a"] + theta[, "b"] + stats::rnorm(n, sd = 0.05)
  fit <- hsic_sensitivity(theta, y, total_order = TRUE,
                          p_value = FALSE, n_permutations = 50L,
                          seed = 13L)
  # On a purely additive model the gap |T - S| should be small per j.
  gap <- abs(fit$index_total_order[, 1L] - fit$index[, 1L])
  expect_true(all(gap < 0.3))
})

test_that("total_order > first_order on a strongly interactive model", {
  # Y = X1 * X2 with X1, X2 ~ U[-1, 1]: pure interaction, marginal
  # E[Y | X_j] = 0 so first-order signal is weaker than total-order.
  set.seed(14)
  n <- 400L
  theta <- matrix(stats::runif(n * 2L, min = -1, max = 1), nrow = n,
                  dimnames = list(NULL, c("a", "b")))
  y <- theta[, "a"] * theta[, "b"] + stats::rnorm(n, sd = 0.05)
  fit <- hsic_sensitivity(theta, y, total_order = TRUE,
                          p_value = FALSE, n_permutations = 50L,
                          seed = 14L)
  # Total-order should pick up both parameters strongly via interaction.
  expect_gt(fit$index_total_order["a", 1L], fit$index["a", 1L])
  expect_gt(fit$index_total_order["b", 1L], fit$index["b", 1L])
})

test_that("total_order = FALSE leaves new fields NULL (backwards-compat)", {
  set.seed(15)
  n <- 60L
  theta <- matrix(stats::runif(n * 2L), nrow = n,
                  dimnames = list(NULL, c("a", "b")))
  y <- theta[, 1L] + stats::rnorm(n, sd = 0.1)
  fit <- hsic_sensitivity(theta, y, p_value = FALSE,
                          n_permutations = 50L, seed = 15L)
  expect_false(fit$total_order)
  expect_null(fit$index_total_order)
  expect_null(fit$statistic_total_order)
})

test_that("total_order is reproducible with seed", {
  set.seed(16)
  n <- 80L
  theta <- matrix(stats::runif(n * 3L), nrow = n,
                  dimnames = list(NULL, c("a", "b", "c")))
  y <- theta[, 1L] + theta[, 2L] * theta[, 3L] +
       stats::rnorm(n, sd = 0.1)
  f1 <- hsic_sensitivity(theta, y, total_order = TRUE,
                         p_value = FALSE,
                         n_permutations = 49L, seed = 7L)
  f2 <- hsic_sensitivity(theta, y, total_order = TRUE,
                         p_value = FALSE,
                         n_permutations = 49L, seed = 7L)
  expect_identical(f1$index_total_order, f2$index_total_order)
})

test_that("print emits total-order block when total_order = TRUE", {
  set.seed(17)
  n <- 60L
  theta <- matrix(stats::runif(n * 2L), nrow = n,
                  dimnames = list(NULL, c("a", "b")))
  y <- theta[, 1L] + stats::rnorm(n, sd = 0.1)
  fit <- hsic_sensitivity(theta, y, total_order = TRUE,
                          p_value = FALSE, n_permutations = 50L,
                          seed = 17L)
  expect_output(print(fit), "T_total_max")
  expect_output(print(fit), "interaction")
})

test_that("plot which='total' / 'both' work + reject if total_order=FALSE", {
  set.seed(18)
  n <- 60L
  theta <- matrix(stats::runif(n * 2L), nrow = n,
                  dimnames = list(NULL, c("a", "b")))
  y <- theta[, 1L] + theta[, 2L] + stats::rnorm(n, sd = 0.1)
  fit_first <- hsic_sensitivity(theta, y, p_value = FALSE,
                                n_permutations = 50L, seed = 18L)
  fit_both <- hsic_sensitivity(theta, y, total_order = TRUE,
                               p_value = FALSE, n_permutations = 50L,
                               seed = 18L)

  expect_error(plot(fit_first, which = "total"),
               "total_order = TRUE")
  expect_error(plot(fit_first, which = "both"),
               "total_order = TRUE")

  expect_silent({
    grDevices::pdf(tempfile(fileext = ".pdf"))
    on.exit(grDevices::dev.off(), add = TRUE)
    plot(fit_both, which = "total")
    plot(fit_both, which = "both")
  })
})

test_that("as.data.frame includes total-order columns when present", {
  set.seed(19)
  n <- 60L
  theta <- matrix(stats::runif(n * 2L), nrow = n,
                  dimnames = list(NULL, c("a", "b")))
  y <- theta[, 1L] + stats::rnorm(n, sd = 0.1)
  fit <- hsic_sensitivity(theta, y, total_order = TRUE,
                          p_value = TRUE, n_permutations = 50L,
                          seed = 19L)
  df <- as.data.frame(fit)
  expect_true("index_total_order" %in% colnames(df))
  expect_true("statistic_total_order" %in% colnames(df))
  expect_equal(nrow(df), 2L)
})

# ---- Total-order pair-bootstrap CI (0.0.0.9013+) ----------------------
# 0.0.0.9012 shipped a `total_order_p_value` field that was not
# null-calibrated (under pure-noise y all parameters got p ~= 0.01).
# Critical review 2026-05-16 retracted the claim; 0.0.0.9013 replaces
# it with `total_order_ci` (CI on the index itself, not a hypothesis
# test). Tests below reflect the retraction.

test_that("total_order_ci populates CI fields with no p-value claims", {
  set.seed(31)
  n <- 60L
  theta <- matrix(stats::runif(n * 3L), nrow = n,
                  dimnames = list(NULL, c("active", "weak", "inert")))
  y <- 2 * theta[, "active"] + 0.5 * theta[, "weak"] +
       stats::rnorm(n, sd = 0.05)
  fit <- hsic_sensitivity(theta, y,
                          total_order = TRUE,
                          total_order_ci = TRUE,
                          n_permutations = 50L,
                          n_bootstrap = 50L,
                          seed = 31L)
  expect_true(fit$total_order_ci)
  expect_equal(fit$n_bootstrap, 50L)
  expect_equal(dim(fit$ci_total_order_lower), c(3L, 1L))
  expect_equal(dim(fit$ci_total_order_upper), c(3L, 1L))
  expect_true(all(fit$ci_total_order_lower <=
                  fit$ci_total_order_upper))
  expect_true(all(fit$ci_total_order_lower >= 0))
  expect_true(all(fit$ci_total_order_upper <= 1))
  # No mis-named p-value fields lurking on the object.
  expect_null(fit$p_value_total_order)
  expect_null(fit$p_value_total_order_adjusted)
})

test_that("total-order CIs are absent by default (backwards compat)", {
  set.seed(32)
  n <- 50L
  theta <- matrix(stats::runif(n * 2L), nrow = n,
                  dimnames = list(NULL, c("a", "b")))
  y <- theta[, 1L] + stats::rnorm(n, sd = 0.1)
  fit <- hsic_sensitivity(theta, y, total_order = TRUE,
                          n_permutations = 50L, seed = 32L)
  expect_false(fit$total_order_ci)
  expect_null(fit$ci_total_order_lower)
  expect_null(fit$ci_total_order_upper)
  expect_null(fit$p_value_total_order)
})

test_that("total_order_ci requires total_order = TRUE", {
  set.seed(33)
  n <- 40L
  theta <- matrix(stats::runif(n * 2L), n, 2L)
  y <- theta[, 1L] + stats::rnorm(n, sd = 0.1)
  expect_error(
    hsic_sensitivity(theta, y, total_order = FALSE,
                     total_order_ci = TRUE,
                     n_permutations = 30L),
    "requires `total_order = TRUE`"
  )
})

test_that("pair-bootstrap CI is reproducible under seed", {
  set.seed(34)
  n <- 50L
  theta <- matrix(stats::runif(n * 2L), n, 2L,
                  dimnames = list(NULL, c("a", "b")))
  y <- theta[, 1L] + stats::rnorm(n, sd = 0.1)
  f1 <- hsic_sensitivity(theta, y, total_order = TRUE,
                         total_order_ci = TRUE,
                         n_permutations = 30L, n_bootstrap = 50L,
                         seed = 7L)
  f2 <- hsic_sensitivity(theta, y, total_order = TRUE,
                         total_order_ci = TRUE,
                         n_permutations = 30L, n_bootstrap = 50L,
                         seed = 7L)
  expect_identical(f1$ci_total_order_lower,   f2$ci_total_order_lower)
  expect_identical(f1$ci_total_order_upper,   f2$ci_total_order_upper)
})

test_that("n_bootstrap and ci_level are validated", {
  set.seed(36)
  n <- 40L
  theta <- matrix(stats::runif(n * 2L), n, 2L)
  y <- theta[, 1L] + stats::rnorm(n, sd = 0.1)
  expect_error(
    hsic_sensitivity(theta, y, total_order = TRUE,
                     total_order_ci = TRUE,
                     n_bootstrap = 0L,
                     n_permutations = 30L),
    "n_bootstrap"
  )
  expect_error(
    hsic_sensitivity(theta, y, total_order = TRUE,
                     total_order_ci = TRUE,
                     ci_level = 1.5,
                     n_permutations = 30L),
    "ci_level"
  )
})

test_that("defunct total_order_p_value arg errors with a pointer", {
  set.seed(37)
  n <- 40L
  theta <- matrix(stats::runif(n * 2L), n, 2L)
  y <- theta[, 1L] + stats::rnorm(n, sd = 0.1)
  expect_error(
    hsic_sensitivity(theta, y, total_order = TRUE,
                     total_order_p_value = TRUE,
                     n_permutations = 30L),
    "defunct"
  )
  expect_error(
    hsic_sensitivity(theta, y, total_order = TRUE,
                     total_order_p_value = TRUE,
                     n_permutations = 30L),
    "total_order_ci"
  )
})

# Critical-review regression test (2026-05-16): document the
# 0.0.0.9012 mis-calibration failure mode so it cannot silently
# return. Under pure-noise y, the bootstrap CI on T_j is dragged
# toward 1 by complement-formulation finite-sample bias -- this is
# expected behaviour for the *index*, but it is NOT a significance
# verdict, and 0.0.0.9013 no longer reports one.
test_that("under pure noise, no significance fields are returned", {
  skip_on_cran()
  set.seed(2026L)
  n <- 80L
  theta <- matrix(stats::runif(n * 3L), n, 3L,
                  dimnames = list(NULL, c("a", "b", "c")))
  y_noise <- stats::rnorm(n)
  fit <- hsic_sensitivity(theta, y_noise,
                          total_order = TRUE,
                          total_order_ci = TRUE,
                          n_permutations = 50L,
                          n_bootstrap = 80L,
                          seed = 11L)
  # No claim of significance survives the retraction.
  expect_null(fit$p_value_total_order)
  expect_null(fit$p_value_total_order_adjusted)
  # CI is allowed to be high (complement-formulation finite-sample
  # bias is a known property, not a methodological claim): the test
  # documents the behaviour rather than passing it off as 'no effect'.
  expect_true(all(fit$ci_total_order_lower >= 0))
  expect_true(all(fit$ci_total_order_upper <= 1))
})

# ---- Conditional-permutation total-order significance (0.0.0.9014+) --

test_that("cond_perm test populates p-value fields with calibration flag", {
  set.seed(41)
  n <- 60L
  theta <- matrix(stats::runif(n * 3L), n, 3L,
                  dimnames = list(NULL, c("active", "weak", "inert")))
  y <- 2 * theta[, "active"] + 0.5 * theta[, "weak"] +
       stats::rnorm(n, sd = 0.05)
  fit <- hsic_sensitivity(theta, y,
                          total_order      = TRUE,
                          total_order_test = "cond_perm",
                          n_permutations   = 50L,
                          n_clusters_cp    = 5L,
                          seed             = 41L)
  expect_equal(fit$total_order_test, "cond_perm")
  expect_equal(dim(fit$p_value_total_order), c(3L, 1L))
  expect_true(all(fit$p_value_total_order >= 0 &
                  fit$p_value_total_order <= 1))
  expect_true(all(fit$p_value_total_order_adjusted >=
                  fit$p_value_total_order - 1e-12))
})

test_that("cond_perm test requires total_order = TRUE", {
  set.seed(42)
  n <- 40L
  theta <- matrix(stats::runif(n * 2L), n, 2L)
  y <- theta[, 1L] + stats::rnorm(n, sd = 0.1)
  expect_error(
    hsic_sensitivity(theta, y, total_order = FALSE,
                     total_order_test = "cond_perm",
                     n_permutations = 30L),
    "requires `total_order = TRUE`"
  )
})

# Critical-review regression: under pure-noise Y, cond_perm p-values
# must NOT be systematically tiny across all parameters. This is the
# null calibration the 0.0.0.9012 implementation failed.
test_that("cond_perm test does not systematically reject under pure noise", {
  skip_on_cran()
  set.seed(2027L)
  n <- 80L
  theta <- matrix(stats::runif(n * 3L), n, 3L,
                  dimnames = list(NULL, c("a", "b", "c")))
  y_noise <- stats::rnorm(n)
  fit <- hsic_sensitivity(theta, y_noise,
                          total_order      = TRUE,
                          total_order_test = "cond_perm",
                          n_permutations   = 99L,
                          n_clusters_cp    = 5L,
                          p_adjust         = "none",
                          seed             = 11L)
  # Under H_0 (X_j _||_ Y | X_{~j}, here just X_j _||_ Y), p-values are
  # roughly uniform on (0, 1). With 3 parameters x B = 99 permutations
  # the strict critical-review regression is: at least one of the
  # three raw p-values must be >= 0.05 (i.e., not ALL parameters
  # spuriously rejected at alpha = 0.05). The retracted 0.0.0.9012 mode
  # failed this loudly (returned 0.01 for every parameter).
  expect_true(any(fit$p_value_total_order >= 0.05))
})

test_that("defunct arg error mentions both replacements", {
  set.seed(43)
  n <- 40L
  theta <- matrix(stats::runif(n * 2L), n, 2L)
  y <- theta[, 1L] + stats::rnorm(n, sd = 0.1)
  expect_error(
    hsic_sensitivity(theta, y, total_order = TRUE,
                     total_order_p_value = TRUE,
                     n_permutations = 30L),
    "total_order_test"
  )
})
