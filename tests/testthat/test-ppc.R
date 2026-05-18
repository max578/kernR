test_that("pesto_ensemble validates inputs and round-trips", {
  set.seed(1)
  post <- matrix(stats::rnorm(40L), ncol = 2L)
  obs  <- matrix(stats::rnorm(20L), ncol = 2L)
  ens  <- pesto_ensemble(post, obs, metadata = list(year = 2018L))
  expect_s3_class(ens, "pesto_ensemble")
  expect_equal(ens$posterior, post)
  expect_equal(ens$observed, obs)
  expect_equal(ens$metadata$year, 2018L)
  expect_output(print(ens), "PESTO ensemble")
})

test_that("pesto_ensemble rejects bad inputs", {
  expect_error(
    pesto_ensemble(matrix(1, 5L, 2L), matrix(1, 5L, 3L)),
    "same number of columns"
  )
  expect_error(
    pesto_ensemble(matrix(NA_real_, 5L, 2L)),
    "finite"
  )
  expect_error(
    pesto_ensemble(matrix(1, 5L, 2L), metadata = "not a list"),
    "list"
  )
})

test_that("mmd_ppc accepts a calibrated posterior (no rejection)", {
  set.seed(42)
  post <- matrix(stats::rnorm(400L), ncol = 2L)
  obs  <- matrix(stats::rnorm(40L),  ncol = 2L)
  fit <- mmd_ppc(post, obs, n_permutations = 199L, seed = 1L)
  expect_s3_class(fit, "mmd_ppc")
  expect_s3_class(fit, "kernel_test_result")
  expect_equal(fit$method, "MMD PPC")
  expect_equal(fit$n_posterior, 200L)
  expect_equal(fit$n_observed, 20L)
  expect_false(fit$reject)
  expect_gt(fit$p_value, 0.05)
})

test_that("mmd_ppc rejects a mean-shifted posterior", {
  set.seed(42)
  post <- matrix(stats::rnorm(400L), ncol = 2L)
  obs  <- matrix(stats::rnorm(40L, mean = 1.5), ncol = 2L)
  fit <- mmd_ppc(post, obs, n_permutations = 199L, seed = 1L)
  expect_true(fit$reject)
  expect_lt(fit$p_value, 0.05)
})

test_that("mmd_ppc surprise diagnostic is bounded and consistent", {
  set.seed(7)
  post <- matrix(stats::rnorm(200L), ncol = 1L)
  obs  <- matrix(stats::rnorm(20L),  ncol = 1L)
  B <- 99L
  fit <- mmd_ppc(post, obs, n_permutations = B, seed = 7L)
  expect_equal(fit$surprise_bits, -log2(fit$p_value))
  expect_lte(fit$surprise_bits, log2(B + 1L))
  expect_gte(fit$surprise_bits, 0)
})

test_that("mmd_ppc dispatches on pesto_ensemble", {
  set.seed(3)
  post <- matrix(stats::rnorm(400L), ncol = 2L)
  obs  <- matrix(stats::rnorm(40L),  ncol = 2L)
  ens <- pesto_ensemble(post, obs, metadata = list(holdout = 2018L))
  fit <- mmd_ppc(ens, n_permutations = 199L, seed = 3L)
  expect_s3_class(fit, "mmd_ppc")
  expect_equal(fit$pesto_metadata$holdout, 2018L)
  # Dispatch with explicit observed override
  fit2 <- mmd_ppc(ens, observed = obs, n_permutations = 199L, seed = 3L)
  expect_equal(fit$statistic, fit2$statistic)
})

test_that("mmd_ppc.pesto_ensemble errors when observed is missing entirely", {
  ens <- pesto_ensemble(matrix(stats::rnorm(40L), ncol = 2L))
  expect_error(mmd_ppc(ens, n_permutations = 99L, seed = 1L),
               "carries no")
})

test_that("mmd_ppc validates inputs", {
  post <- matrix(stats::rnorm(40L), ncol = 2L)
  expect_error(mmd_ppc(post),
               "observed")
  expect_error(mmd_ppc(post, observed = matrix(0, 20L, 3L)),
               "same number of columns")
  expect_error(mmd_ppc(post, observed = matrix(0, 3L, 2L)),
               "at least 5")
  expect_error(mmd_ppc(post, observed = matrix(0, 20L, 2L), alpha = 0),
               "alpha")
})

test_that("mmd_ppc is reproducible with seed", {
  set.seed(2)
  post <- matrix(stats::rnorm(200L), ncol = 1L)
  obs  <- matrix(stats::rnorm(20L),  ncol = 1L)
  f1 <- mmd_ppc(post, obs, n_permutations = 99L, seed = 99L)
  f2 <- mmd_ppc(post, obs, n_permutations = 99L, seed = 99L)
  expect_identical(f1$statistic, f2$statistic)
  expect_identical(f1$p_value,   f2$p_value)
})

test_that("mmd_ppc print method emits verdict line", {
  set.seed(8)
  post <- matrix(stats::rnorm(200L), ncol = 1L)
  obs  <- matrix(stats::rnorm(20L),  ncol = 1L)
  fit <- mmd_ppc(post, obs, n_permutations = 99L, seed = 8L)
  expect_output(print(fit), "PPC verdict")
  expect_output(print(fit), "Surprise")
})

# ----- mmd_ppc() against the PESTO 0.3.0 manifest contract --------------

.make_demo_manifest <- function(seed = 1L, nreal = 80L, nobs = 3L,
                                obs_target = NULL,
                                post_shift = 0.0) {
  set.seed(seed)
  post <- matrix(stats::rnorm(nreal * nobs), nreal, nobs) + post_shift
  cols <- paste0("o", seq_len(nobs))
  colnames(post) <- cols
  if (is.null(obs_target)) obs_target <- stats::rnorm(nobs)
  names(obs_target) <- cols
  PESTO::pesto_ensemble_manifest(
    run_id          = paste0("ppc_demo_", seed),
    params          = data.frame(real_name = paste0("r", seq_len(nreal)),
                                  p1 = stats::rnorm(nreal),
                                  check.names = FALSE),
    outputs         = data.frame(real_name = paste0("r", seq_len(nreal)),
                                  post, check.names = FALSE),
    weights         = stats::setNames(rep(1, nobs), cols),
    obs_target      = obs_target,
    data_hash       = paste0("sha256:demo_", seed),
    pesto_version   = as.character(utils::packageVersion("PESTO")),
    timestamp       = Sys.time(),
    method          = "ies_callback",
    noptmax         = 1L,
    lambda_schedule = 1
  )
}

test_that("mmd_ppc dispatches on pesto_ensemble_manifest", {
  testthat::skip_if_not_installed("PESTO", minimum_version = "0.3.0")
  m <- .make_demo_manifest(seed = 11L)
  observed <- matrix(stats::rnorm(30L), 10L, 3L,
                     dimnames = list(NULL, paste0("o", 1:3)))
  res <- mmd_ppc(m, observed = observed,
                 n_permutations = 99L, seed = 1L)
  expect_s3_class(res, "mmd_ppc")
  expect_s3_class(res, "kernel_test_result")
  expect_equal(res$pesto_metadata$run_id, "ppc_demo_11")
  expect_setequal(res$pesto_metadata$outputs_used, paste0("o", 1:3))
})

test_that("mmd_ppc rejects under a clear distributional shift", {
  testthat::skip_if_not_installed("PESTO", minimum_version = "0.3.0")
  m <- .make_demo_manifest(seed = 21L, nreal = 200L, nobs = 2L)
  observed <- matrix(stats::rnorm(200L, mean = 3.0),
                     ncol = 2L,
                     dimnames = list(NULL, paste0("o", 1:2)))
  res <- mmd_ppc(m, observed = observed,
                 n_permutations = 199L, seed = 1L)
  expect_lt(res$p_value, 0.05)
  expect_true(isTRUE(res$reject))
})

test_that("mmd_ppc on a manifest errors when `observed` is omitted", {
  testthat::skip_if_not_installed("PESTO", minimum_version = "0.3.0")
  m <- .make_demo_manifest(seed = 31L)
  expect_error(
    mmd_ppc(m, n_permutations = 50L),
    "must be supplied"
  )
})

test_that("mmd_ppc honours outputs subselection on a manifest", {
  testthat::skip_if_not_installed("PESTO", minimum_version = "0.3.0")
  m <- .make_demo_manifest(seed = 41L)
  observed <- matrix(stats::rnorm(20L), 10L, 2L,
                     dimnames = list(NULL, c("o1", "o3")))
  res <- mmd_ppc(m, observed = observed, outputs = c("o1", "o3"),
                 n_permutations = 99L, seed = 1L)
  expect_setequal(res$pesto_metadata$outputs_used, c("o1", "o3"))
})

test_that("mmd_ppc errors on unknown outputs column for manifest", {
  testthat::skip_if_not_installed("PESTO", minimum_version = "0.3.0")
  m <- .make_demo_manifest(seed = 51L)
  observed <- matrix(stats::rnorm(20L), 10L, 2L)
  expect_error(
    mmd_ppc(m, observed = observed,
            outputs = c("o1", "no_such_col"),
            n_permutations = 50L),
    "not found"
  )
})

test_that("mmd_ppc errors on observed/posterior column-count mismatch", {
  testthat::skip_if_not_installed("PESTO", minimum_version = "0.3.0")
  m <- .make_demo_manifest(seed = 61L)
  bad <- matrix(stats::rnorm(10L), 5L, 2L)        # only 2 cols, mismatch
  expect_error(
    mmd_ppc(m, observed = bad, n_permutations = 50L),
    "same number of columns"
  )
})
