# test-auto-switch.R -- the n_exact_max auto-switch on the exact tests.
#
# The exact ksd_test() / concordance_test() delegate to their low-rank
# counterparts above n_exact_max. The contract (Max's confirmed mitigation):
# the switch is announced (message), recorded (approximation / m in the object,
# so the result is reproducible from the object), and escapable (Inf forces the
# exact path). These tests pin all three, plus the no-switch default below the
# ceiling.

test_that("ksd_test auto-switches above n_exact_max: announced + recorded", {
  set.seed(1)
  x <- matrix(stats::rnorm(120L), ncol = 2L)               # n = 60
  expect_message(
    fit <- ksd_test(x, n_exact_max = 30L, n_boot = 99L, seed = 1L),
    "delegating to ksd_test_nystrom"
  )
  expect_identical(fit$approximation, "nystrom")
  expect_identical(fit$m, 60L - 1L)                        # m=100 capped at n-1
  expect_s3_class(fit, c("ksd_test", "kernel_test_result"))
})

test_that("ksd_test stays exact below the ceiling and at Inf", {
  set.seed(2)
  x <- matrix(stats::rnorm(120L), ncol = 2L)               # n = 60
  # Default ceiling (5000) -> exact, no message, no approximation field.
  expect_no_message(fit <- ksd_test(x, n_boot = 49L, seed = 2L))
  expect_null(fit$approximation)
  # A big sample with n_exact_max = Inf stays exact (no delegation message).
  x_big <- matrix(stats::rnorm(12000L), ncol = 2L)         # n = 6000 > 5000
  expect_no_message(
    fit_inf <- ksd_test(x_big, n_exact_max = Inf, n_boot = 1L, seed = 2L)
  )
  expect_null(fit_inf$approximation)
})

test_that("ksd_test validates n_exact_max", {
  x <- matrix(stats::rnorm(120L), ncol = 2L)
  expect_error(ksd_test(x, n_exact_max = 0L), "positive number")
  expect_error(ksd_test(x, n_exact_max = c(10L, 20L)), "positive number")
})

test_that("concordance_test auto-switches above n_exact_max", {
  set.seed(3)
  draws <- list(
    a = matrix(stats::rnorm(120L), ncol = 2L),             # 60 each -> 180
    b = matrix(stats::rnorm(120L), ncol = 2L),
    c = matrix(stats::rnorm(120L), ncol = 2L) + 1
  )
  expect_message(
    fit <- concordance_test(draws, n_exact_max = 100L,
                            n_permutations = 99L, seed = 3L),
    "delegating to concordance_test_nystrom"
  )
  expect_identical(fit$approximation, "nystrom")
  # Pairwise localisation survives the delegated low-rank path.
  expect_gt(fit$pairwise["a", "c"], fit$pairwise["a", "b"])
})

test_that("concordance_test stays exact below the ceiling and at Inf", {
  set.seed(4)
  draws <- list(
    a = matrix(stats::rnorm(200L), ncol = 2L),             # 100 each -> 200
    b = matrix(stats::rnorm(200L), ncol = 2L)
  )
  expect_no_message(fit <- concordance_test(draws, n_permutations = 49L,
                                            seed = 4L))
  expect_null(fit$approximation)
  expect_no_message(
    fit_inf <- concordance_test(draws, n_exact_max = Inf,
                                n_permutations = 1L, seed = 4L)
  )
  expect_null(fit_inf$approximation)
})

test_that("concordance_test validates n_exact_max", {
  draws <- list(
    a = matrix(stats::rnorm(60L), ncol = 2L),
    b = matrix(stats::rnorm(60L), ncol = 2L)
  )
  expect_error(concordance_test(draws, n_exact_max = -1), "positive number")
})
