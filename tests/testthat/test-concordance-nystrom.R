# test-concordance-nystrom.R -- low-rank (Nystrom / RFF) concordance test.
#
# Anchors concordance_test_nystrom() against the exact concordance_test():
# reconstruction as m -> n, source localisation, reproducibility, the RFF
# path, power under a real departure, and -- gated for cost -- the Type-I
# error rate under the null. The reconstruction test is the correctness
# anchor; everything else guards behaviour the acceleration must preserve.

test_that("nystrom concordance reconstructs the exact statistic as m -> n", {
  set.seed(1)
  draws <- list(
    a = matrix(stats::rnorm(120L), ncol = 2L),
    b = matrix(stats::rnorm(120L), ncol = 2L) + 0.5,
    c = matrix(stats::rnorm(120L), ncol = 2L)
  )
  n_total <- sum(vapply(draws, nrow, integer(1L)))

  exact <- concordance_test(draws, n_permutations = 1L, seed = 1L)
  approx <- concordance_test_nystrom(draws, m = n_total - 1L,
                                     n_permutations = 1L, seed = 1L)

  # Statistic and full pairwise matrix coincide to Nystrom residual.
  expect_equal(approx$statistic, exact$statistic, tolerance = 1e-4)
  expect_equal(approx$pairwise, exact$pairwise, tolerance = 1e-4)
})

test_that("pairwise matrix localises the departing source", {
  set.seed(2)
  draws <- list(
    engine_a = matrix(stats::rnorm(600L), ncol = 2L),
    engine_b = matrix(stats::rnorm(600L), ncol = 2L),
    engine_c = matrix(stats::rnorm(600L), ncol = 2L) + 1
  )
  fit <- concordance_test_nystrom(draws, m = 60L,
                                  n_permutations = 99L, seed = 2L)

  # The shifted source carries the largest off-diagonal discrepancies.
  ac <- fit$pairwise["engine_a", "engine_c"]
  bc <- fit$pairwise["engine_b", "engine_c"]
  ab <- fit$pairwise["engine_a", "engine_b"]
  expect_gt(ac, ab)
  expect_gt(bc, ab)
  expect_true(fit$reject)
})

test_that("nystrom concordance is reproducible with seed", {
  set.seed(3)
  draws <- list(
    a = matrix(stats::rnorm(400L), ncol = 2L),
    b = matrix(stats::rnorm(400L), ncol = 2L) + 0.3
  )
  f1 <- concordance_test_nystrom(draws, m = 50L,
                                 n_permutations = 99L, seed = 7L)
  f2 <- concordance_test_nystrom(draws, m = 50L,
                                 n_permutations = 99L, seed = 7L)
  expect_identical(f1$statistic, f2$statistic)
  expect_identical(f1$p_value, f2$p_value)
  expect_identical(f1$null_distribution, f2$null_distribution)
})

test_that("the RFF path runs and localises a departure", {
  set.seed(4)
  draws <- list(
    a = matrix(stats::rnorm(500L), ncol = 2L),
    b = matrix(stats::rnorm(500L), ncol = 2L) + 0.8
  )
  fit <- concordance_test_nystrom(draws, method = "rff", m = 200L,
                                  n_permutations = 99L, seed = 4L)
  expect_identical(fit$approximation, "rff")
  expect_gt(fit$statistic, 0)
  expect_true(fit$reject)
})

test_that("RFF rejects a non-RBF kernel", {
  draws <- list(
    a = matrix(stats::rnorm(60L), ncol = 2L),
    b = matrix(stats::rnorm(60L), ncol = 2L)
  )
  expect_error(
    concordance_test_nystrom(draws, kernel = kernel_spec("linear"),
                             method = "rff", m = 20L, seed = 1L),
    "RBF"
  )
})

test_that("input validation is inherited from concordance_test", {
  expect_error(
    concordance_test_nystrom(list(matrix(stats::rnorm(20L), ncol = 2L))),
    "at least two"
  )
  expect_error(
    concordance_test_nystrom(
      list(a = matrix(stats::rnorm(8L), ncol = 2L),
           b = matrix(stats::rnorm(6L), ncol = 2L)),
      m = 10L
    ),
    "at least 5 rows"
  )
})

test_that("the object class and fields match the exact test", {
  set.seed(5)
  draws <- list(
    a = matrix(stats::rnorm(200L), ncol = 2L),
    b = matrix(stats::rnorm(200L), ncol = 2L)
  )
  fit <- concordance_test_nystrom(draws, m = 40L,
                                  n_permutations = 49L, seed = 5L)
  expect_s3_class(fit, c("concordance_test", "kernel_test_result"))
  expect_identical(fit$n_groups, 2L)
  expect_identical(fit$approximation, "nystrom")
  expect_true(is.finite(fit$surprise_bits))
})

test_that("Type-I error rate stays near the nominal level under the null", {
  skip_on_cran()
  set.seed(11)
  alpha <- 0.1
  n_trials <- 80L
  reject <- logical(n_trials)
  for (t in seq_len(n_trials)) {
    draws <- list(
      a = matrix(stats::rnorm(300L), ncol = 2L),
      b = matrix(stats::rnorm(300L), ncol = 2L),
      c = matrix(stats::rnorm(300L), ncol = 2L)
    )
    fit <- concordance_test_nystrom(draws, m = 50L, alpha = alpha,
                                    n_permutations = 99L)
    reject[t] <- fit$reject
  }
  # Permutation test => the empirical size should sit near alpha, not blow
  # up. Allow Monte-Carlo slack on 80 trials.
  expect_lt(mean(reject), alpha + 0.12)
})
