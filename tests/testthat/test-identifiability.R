test_that("lhs_design returns correct shape and is reproducible", {
  bounds <- rbind(a = c(0, 1), b = c(-2, 2), c = c(10, 20))
  d1 <- lhs_design(50, bounds, seed = 1)
  d2 <- lhs_design(50, bounds, seed = 1)

  expect_equal(dim(d1), c(50L, 3L))
  expect_equal(colnames(d1), c("a", "b", "c"))
  expect_identical(d1, d2)
})

test_that("lhs_design respects bounds and stratifies coverage", {
  bounds <- rbind(a = c(0, 1), b = c(5, 10))
  n <- 40L
  d <- lhs_design(n, bounds, seed = 2)

  expect_true(all(d[, "a"] >= 0 & d[, "a"] <= 1))
  expect_true(all(d[, "b"] >= 5 & d[, "b"] <= 10))

  # Each of n equal bins should contain exactly one point per column.
  bin_a <- findInterval(d[, "a"], seq(0, 1, length.out = n + 1L),
                        rightmost.closed = TRUE)
  expect_equal(sort(bin_a), seq_len(n))
})

test_that("lhs_design validates inputs", {
  expect_error(lhs_design(1, rbind(c(0, 1))), "integer")
  expect_error(lhs_design(5, rbind(c(1, 0))), "upper bound")
  expect_error(lhs_design(5, rbind(c(0, Inf))), "finite")
  expect_error(lhs_design(5, matrix(0, 2, 3)), "two columns")
})

test_that("hsic_identifiability flags active vs inert parameters", {
  set.seed(7)
  n <- 80L
  bounds <- rbind(
    active1 = c(0, 1),
    active2 = c(-1, 1),
    inert   = c(0, 1)
  )
  theta <- lhs_design(n, bounds, seed = 7)
  y1 <- 2 * theta[, "active1"] + stats::rnorm(n, sd = 0.05)
  y2 <- theta[, "active2"]^2 + stats::rnorm(n, sd = 0.05)
  y  <- cbind(yield = y1, biomass = y2)

  fit <- hsic_identifiability(theta, y,
                              n_permutations = 199L, seed = 7)
  expect_s3_class(fit, "hsic_identifiability")
  expect_equal(dim(fit$statistic), c(3L, 2L))
  expect_equal(dim(fit$p_value), c(3L, 2L))
  expect_true(fit$identifiable["active1"])
  expect_true(fit$identifiable["active2"])
  expect_false(fit$identifiable["inert"])
})

test_that("hsic_identifiability is reproducible with seed", {
  set.seed(1)
  n <- 60L
  theta <- matrix(stats::runif(n * 2L), nrow = n,
                  dimnames = list(NULL, c("p1", "p2")))
  y <- cbind(o = theta[, 1] + stats::rnorm(n, sd = 0.1))
  f1 <- hsic_identifiability(theta, y, n_permutations = 99L, seed = 42L)
  f2 <- hsic_identifiability(theta, y, n_permutations = 99L, seed = 42L)
  expect_identical(f1$statistic, f2$statistic)
  expect_identical(f1$p_value,   f2$p_value)
})

test_that("hsic_identifiability accepts vector y and supports p_adjust='none'", {
  set.seed(2)
  n <- 60L
  theta <- matrix(stats::runif(n * 2L), nrow = n,
                  dimnames = list(NULL, c("p1", "p2")))
  y <- theta[, 1] + stats::rnorm(n, sd = 0.1)
  fit <- hsic_identifiability(theta, y, p_adjust = "none",
                              n_permutations = 99L, seed = 2L)
  expect_equal(fit$p_adjust, "none")
  expect_identical(fit$p_value, fit$p_value_adjusted)
  expect_equal(dim(fit$statistic), c(2L, 1L))
})

test_that("hsic_identifiability validates inputs", {
  expect_error(
    hsic_identifiability(matrix(0, 5L, 2L), matrix(0, 5L, 1L)),
    "At least 10"
  )
  expect_error(
    hsic_identifiability(matrix(0, 20L, 2L), matrix(0, 10L, 1L)),
    "same number of rows"
  )
  expect_error(
    hsic_identifiability(matrix(0, 20L, 2L), matrix(0, 20L, 1L), alpha = 1.5),
    "alpha"
  )
})

test_that("print, summary, plot and as.data.frame work", {
  set.seed(3)
  n <- 40L
  theta <- matrix(stats::runif(n * 2L), nrow = n,
                  dimnames = list(NULL, c("a", "b")))
  y <- cbind(o1 = theta[, 1] + stats::rnorm(n, sd = 0.2),
             o2 = stats::rnorm(n))
  fit <- hsic_identifiability(theta, y, n_permutations = 99L, seed = 3L)

  expect_output(print(fit), "Identifiability")
  expect_output(summary(fit), "Summary")
  expect_silent({
    grDevices::pdf(tempfile(fileext = ".pdf"))
    on.exit(grDevices::dev.off(), add = TRUE)
    plot(fit)
  })

  df <- as.data.frame(fit)
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 4L)
  expect_setequal(colnames(df),
                  c("parameter", "output", "statistic",
                    "p_value", "p_value_adjusted"))
})
