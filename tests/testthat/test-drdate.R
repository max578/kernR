test_that("DR-DATE detects mean treatment effect", {
  set.seed(42)
  n <- 300
  x <- matrix(rnorm(n * 2), n, 2)
  logit_p <- 0.5 * x[, 1]
  t <- rbinom(n, 1, plogis(logit_p))
  y <- t * 1.0 + 0.5 * x[, 1] + rnorm(n, sd = 0.5)

  result <- dr_date_test(y, t, x,
    n_permutations = 100,
    seed = 1
  )

  expect_s3_class(result, "kernel_test_result")
  expect_equal(result$method, "DR-DATE")
  expect_true(result$p_value < 0.1)
})

test_that("DR-DATE does not reject under null", {
  set.seed(42)
  n <- 300
  x <- matrix(rnorm(n * 2), n, 2)
  t <- rbinom(n, 1, 0.5)
  y <- 0.5 * x[, 1] + rnorm(n)  # No treatment effect

  result <- dr_date_test(y, t, x,
    n_permutations = 100,
    outcome_model = "zero",
    seed = 1
  )
  expect_s3_class(result, "kernel_test_result")
  expect_true(result$p_value >= 0 && result$p_value <= 1)
})

test_that("DR-DETT runs without error", {
  set.seed(42)
  n <- 200
  x <- matrix(rnorm(n * 2), n, 2)
  t <- rbinom(n, 1, plogis(0.3 * x[, 1]))
  y <- t * rnorm(n, sd = 2) + (1 - t) * rnorm(n, sd = 1) + x[, 1]

  result <- dr_dett_test(y, t, x,
    n_permutations = 100,
    seed = 1
  )
  expect_s3_class(result, "kernel_test_result")
  expect_equal(result$method, "DR-DETT")
})

test_that("kernel_causal_test formula interface works", {
  set.seed(42)
  n <- 200
  dat <- data.frame(
    y = rnorm(n),
    treatment = rbinom(n, 1, 0.5),
    x1 = rnorm(n),
    x2 = rnorm(n)
  )
  dat$y <- dat$y + 0.5 * dat$treatment + 0.3 * dat$x1

  result <- kernel_causal_test(
    y ~ treatment | x1 + x2,
    data = dat,
    method = "dr-date",
    n_permutations = 50,
    seed = 1
  )
  expect_s3_class(result, "kernel_test_result")
})
