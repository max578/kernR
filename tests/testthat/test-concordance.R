# Tests for the kernel k-sample concordance test (concordance_test): mutual
# concordance of several sources, power to detect and localise a divergent
# source, the pairwise discrepancy matrix, labels, reproducibility, and input
# validation.

test_that("concordant sources are not rejected", {
  set.seed(20)
  draws <- list(
    a = matrix(stats::rnorm(400L), ncol = 2L),
    b = matrix(stats::rnorm(400L), ncol = 2L),
    c = matrix(stats::rnorm(400L), ncol = 2L)
  )
  fit <- concordance_test(draws, n_permutations = 199L, seed = 1L)

  expect_s3_class(fit, c("concordance_test", "kernel_test_result"))
  expect_gt(fit$p_value, 0.05)
  expect_false(fit$reject)
  expect_equal(fit$n_groups, 3L)
  expect_equal(unname(fit$group_sizes), c(200L, 200L, 200L))
})

test_that("a mean-divergent source is rejected and localised", {
  set.seed(21)
  draws <- list(
    a = matrix(stats::rnorm(400L), ncol = 2L),
    b = matrix(stats::rnorm(400L), ncol = 2L),
    c = matrix(stats::rnorm(400L), ncol = 2L) + 1
  )
  fit <- concordance_test(draws, n_permutations = 199L, seed = 1L)

  expect_lt(fit$p_value, 0.05)
  expect_true(fit$reject)

  # The divergent source c carries the largest pairwise discrepancies
  pw <- fit$pairwise
  expect_gt(pw["a", "c"], pw["a", "b"])
  expect_gt(pw["b", "c"], pw["a", "b"])
})

test_that("a variance-divergent source is rejected", {
  set.seed(22)
  draws <- list(
    a = matrix(stats::rnorm(400L), ncol = 2L),
    b = matrix(stats::rnorm(400L), ncol = 2L),
    c = matrix(stats::rnorm(400L, sd = 2), ncol = 2L)
  )
  fit <- concordance_test(draws, n_permutations = 199L, seed = 1L)
  expect_lt(fit$p_value, 0.05)
})

test_that("the pairwise matrix is symmetric with a zero diagonal", {
  set.seed(23)
  draws <- list(
    a = matrix(stats::rnorm(400L), ncol = 2L),
    b = matrix(stats::rnorm(400L), ncol = 2L)
  )
  fit <- concordance_test(draws, n_permutations = 99L, seed = 1L)
  pw <- fit$pairwise
  expect_equal(pw, t(pw))
  expect_equal(diag(pw), c(a = 0, b = 0))
})

test_that("unnamed lists get default source labels", {
  set.seed(24)
  draws <- list(
    matrix(stats::rnorm(200L), ncol = 2L),
    matrix(stats::rnorm(200L), ncol = 2L)
  )
  fit <- concordance_test(draws, n_permutations = 99L, seed = 1L)
  expect_equal(names(fit$group_sizes), c("Source 1", "Source 2"))
})

test_that("the verdict is reproducible under seed and the print method runs", {
  set.seed(25)
  draws <- list(
    a = matrix(stats::rnorm(300L), ncol = 3L),
    b = matrix(stats::rnorm(300L), ncol = 3L)
  )
  a <- concordance_test(draws, n_permutations = 199L, seed = 7L)
  b <- concordance_test(draws, n_permutations = 199L, seed = 7L)
  expect_identical(a$p_value, b$p_value)
  expect_identical(a$statistic, b$statistic)
  expect_output(print(a), "Concordance verdict")
  expect_output(print(a), "Pairwise")
})

test_that("concordance_test() rejects malformed input", {
  good <- matrix(stats::rnorm(200L), ncol = 2L)

  expect_error(concordance_test(list(good)), "at least two samples")
  expect_error(concordance_test(good), "list of at least two")
  expect_error(
    concordance_test(list(a = good, b = matrix(stats::rnorm(300L), ncol = 3L))),
    "same number of columns"
  )
  expect_error(
    concordance_test(list(a = good, b = matrix(stats::rnorm(6L), ncol = 2L))),
    "at least 5 rows"
  )
  bad <- good
  bad[1L, 1L] <- NA_real_
  expect_error(concordance_test(list(a = good, b = bad)), "finite")
  expect_error(concordance_test(list(a = good, b = good), alpha = 1.5),
               "in \\(0, 1\\)")
})
