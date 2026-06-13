# Tests for dr_date_scenario() — Year-1 §B2 of the UQ ag-stack roadmap.
# Requires PESTO (>= 0.3.0) — wired through Imports.

skip_unless_pesto <- function() {
  testthat::skip_if_not_installed("PESTO", minimum_version = "0.3.0")
}

# Build a pair of PESTO manifests sharing one calibrated posterior and
# differing only in the *forward-simulated* outputs under two scenarios.
# This is the canonical ag-systems "did intervention shift outputs?"
# setup: identical parameters (positivity satisfied by construction),
# different outputs encode the intervention's effect.
.make_scenario_pair <- function(intervention_shift = 0.0,
                                seed = 42L, nreal = 60L,
                                npar = 2L, nobs = 4L, sigma = 0.05,
                                par_names = c("p1", "p2"),
                                obs_names = paste0("o", seq_len(nobs))) {
  set.seed(seed)
  G   <- matrix(stats::rnorm(nobs * npar), nobs, npar)
  y0  <- as.numeric(G %*% c(1.0, -0.5)) +
         stats::rnorm(nobs, sd = sigma)
  names(y0) <- obs_names

  prior <- matrix(stats::rnorm(nreal * npar), nreal, npar,
                  dimnames = list(NULL, par_names))

  fit <- PESTO::pesto_ies_callback(
    forward_model  = function(t) t %*% t(G),
    prior_ensemble = prior, obs = y0, obs_sd = sigma,
    noptmax = 2L, verbose = FALSE
  )

  # One posterior; forward-simulate it twice (baseline + intervention)
  par_post <- as.matrix(fit$par_ensemble[, par_names, with = FALSE])
  out_base <- par_post %*% t(G)
  out_intv <- out_base + intervention_shift
  colnames(out_base) <- colnames(out_intv) <- obs_names

  real_names <- fit$par_ensemble$real_name

  params_df <- data.frame(real_name = real_names,
                          par_post, stringsAsFactors = FALSE,
                          check.names = FALSE)
  outputs_base_df <- data.frame(real_name = real_names,
                                out_base, stringsAsFactors = FALSE,
                                check.names = FALSE)
  outputs_intv_df <- data.frame(real_name = real_names,
                                out_intv, stringsAsFactors = FALSE,
                                check.names = FALSE)

  weights    <- stats::setNames(rep(1 / sigma, nobs), obs_names)
  obs_target <- stats::setNames(y0, obs_names)
  ts         <- Sys.time()

  m_base <- PESTO::pesto_ensemble_manifest(
    run_id          = "test_baseline",
    params          = params_df,
    outputs         = outputs_base_df,
    weights         = weights,
    obs_target      = obs_target,
    data_hash       = "sha256:test_only_baseline",
    pesto_version   = as.character(utils::packageVersion("PESTO")),
    timestamp       = ts,
    method          = "ies_callback",
    noptmax         = 1L,
    lambda_schedule = 1
  )
  m_intv <- PESTO::pesto_ensemble_manifest(
    run_id          = "test_intervention",
    params          = params_df,                    # SAME params
    outputs         = outputs_intv_df,
    weights         = weights,
    obs_target      = obs_target,
    data_hash       = "sha256:test_only_intervention",
    pesto_version   = as.character(utils::packageVersion("PESTO")),
    timestamp       = ts,
    method          = "ies_callback",
    noptmax         = 1L,
    lambda_schedule = 1
  )
  list(baseline = m_base, intervention = m_intv)
}

test_that("dr_date_scenario rejects under a real distributional shift", {
  skip_unless_pesto()
  pair <- .make_scenario_pair(intervention_shift = 0.6, seed = 11L)
  res <- dr_date_scenario(
    baseline       = pair$baseline,
    intervention   = pair$intervention,
    n_permutations = 199L,
    seed           = 1L,
    verbose        = FALSE
  )
  expect_s3_class(res, "dr_date_scenario")
  expect_s3_class(res, "kernel_test_result")
  expect_equal(res$method, "DR-DATE (scenario)")
  expect_equal(res$n_baseline, 60L)
  expect_equal(res$n_intervention, 60L)
  expect_equal(res$baseline_run_id, "test_baseline")
  expect_equal(res$intervention_run_id, "test_intervention")
  expect_setequal(res$outputs_tested, paste0("o", 1:4))
  expect_lt(res$p_value, 0.05)
  expect_true(isTRUE(res$reject))
})

test_that("dr_date_scenario fails to reject when scenarios are identical", {
  skip_unless_pesto()
  # Null case: shift = 0 — same observations, so PESTO calibrations
  # produce identical posteriors; outputs come from the same
  # distribution. Different seeds for each replicate.
  set.seed(99L)
  null_pvals <- vapply(seq_len(6L), function(rep) {
    pair <- .make_scenario_pair(intervention_shift = 0.0,
                                seed = 100L + rep)
    res  <- dr_date_scenario(pair$baseline, pair$intervention,
                              n_permutations = 99L,
                              seed = rep)
    res$p_value
  }, numeric(1L))
  # Under the null, p-values should be roughly uniform — mean ≈ 0.5,
  # well above the 0.05 rejection threshold on average.
  expect_gt(mean(null_pvals), 0.20)
})

test_that("dr_date_scenario validates manifest pair classes", {
  skip_unless_pesto()
  pair <- .make_scenario_pair(seed = 1L)
  expect_error(
    dr_date_scenario(baseline = list(),
                     intervention = pair$intervention,
                     n_permutations = 50L),
    "pesto_ensemble_manifest"
  )
  expect_error(
    dr_date_scenario(baseline = pair$baseline,
                     intervention = "not_a_manifest",
                     n_permutations = 50L),
    "pesto_ensemble_manifest"
  )
})

test_that("dr_date_scenario rejects parameter-schema mismatch", {
  skip_unless_pesto()
  pair <- .make_scenario_pair(seed = 2L)
  # Build a manifest with a different parameter name set
  set.seed(3L)
  G  <- matrix(stats::rnorm(8L), 4L, 2L)
  y  <- as.numeric(G %*% c(0.1, 0.2)) + stats::rnorm(4L, sd = 0.05)
  names(y) <- paste0("o", 1:4)
  prior_bad <- matrix(stats::rnorm(60L * 2L), 60L, 2L,
                      dimnames = list(NULL, c("Q1", "Q2")))   # different par names
  fit_bad <- PESTO::pesto_ies_callback(
    forward_model  = function(t) t %*% t(G),
    prior_ensemble = prior_bad, obs = y, obs_sd = 0.05,
    noptmax = 2L, verbose = FALSE
  )
  m_bad <- PESTO::as_manifest(fit_bad, run_id = "bad_params")
  expect_error(
    dr_date_scenario(pair$baseline, m_bad, n_permutations = 50L),
    "parameter schemas differ"
  )
})

test_that("dr_date_scenario rejects output-schema mismatch", {
  skip_unless_pesto()
  pair <- .make_scenario_pair(seed = 4L, nobs = 4L)
  # Build a manifest with different obs count
  set.seed(5L)
  G  <- matrix(stats::rnorm(6L), 3L, 2L)   # 3 obs (not 4)
  y  <- as.numeric(G %*% c(0.1, 0.2)) + stats::rnorm(3L, sd = 0.05)
  names(y) <- paste0("o", 1:3)
  prior <- matrix(stats::rnorm(60L * 2L), 60L, 2L,
                  dimnames = list(NULL, c("p1", "p2")))
  fit_bad <- PESTO::pesto_ies_callback(
    forward_model  = function(t) t %*% t(G),
    prior_ensemble = prior, obs = y, obs_sd = 0.05,
    noptmax = 2L, verbose = FALSE
  )
  m_bad <- PESTO::as_manifest(fit_bad, run_id = "bad_outputs")
  expect_error(
    dr_date_scenario(pair$baseline, m_bad, n_permutations = 50L),
    "observation schemas differ"
  )
})

test_that("dr_date_scenario honours output subselection", {
  skip_unless_pesto()
  pair <- .make_scenario_pair(intervention_shift = 0.6, seed = 17L)
  res  <- dr_date_scenario(
    pair$baseline, pair$intervention,
    output = c("o1", "o3"),
    n_permutations = 99L, seed = 1L
  )
  expect_setequal(res$outputs_tested, c("o1", "o3"))
})

test_that("dr_date_scenario errors on unknown output column", {
  skip_unless_pesto()
  pair <- .make_scenario_pair(seed = 7L)
  expect_error(
    dr_date_scenario(pair$baseline, pair$intervention,
                     output = c("o1", "no_such_column"),
                     n_permutations = 50L),
    "not found"
  )
})

test_that("print.dr_date_scenario emits expected verdict header", {
  skip_unless_pesto()
  pair <- .make_scenario_pair(intervention_shift = 0.6, seed = 19L)
  res  <- dr_date_scenario(pair$baseline, pair$intervention,
                            n_permutations = 99L, seed = 1L)
  out <- utils::capture.output(print(res))
  expect_true(any(grepl("DR-DATE", out)))
  expect_true(any(grepl("baseline", out)))
  expect_true(any(grepl("intervention", out)))
  expect_true(any(grepl("test_baseline", out)))
  expect_true(any(grepl("PESTO versions", out)))
})

# --- Independent Oracle Principle: grounding refusals (Phase-0 FX-2/FX-3) ----

test_that(".validate_manifest_pair refuses an incompatible APSIM major (FX-3)", {
  skip_unless_pesto()
  pair <- .make_scenario_pair(intervention_shift = 0.3, seed = 7L)
  b <- pair$baseline;     b@apsim_version <- "2024.1.7000.0"
  i <- pair$intervention; i@apsim_version <- "2026.5.8046.0"
  expect_error(kernR:::.validate_manifest_pair(b, i),
               "incompatible APSIM majors")
  # Same major passes the APSIM gate.
  b2 <- pair$baseline;     b2@apsim_version <- "2024.1.7000.0"
  i2 <- pair$intervention; i2@apsim_version <- "2024.9.9999.0"
  expect_silent(kernR:::.validate_manifest_pair(b2, i2))
})

test_that(".validate_manifest_pair refuses an obs_schema unit disagreement (FX-2)", {
  skip_unless_pesto()
  skip_if_not_installed("PESTO", minimum_version = "0.6.0.9000")
  pair <- .make_scenario_pair(intervention_shift = 0.3, seed = 8L)
  os_t  <- PESTO::pesto_obs_schema(
    outputs = data.frame(name = paste0("o", 1:4),
                         quantity = rep("grain_yield", 4),
                         unit = rep("t/ha", 4), stringsAsFactors = FALSE))
  os_kg <- PESTO::pesto_obs_schema(
    outputs = data.frame(name = paste0("o", 1:4),
                         quantity = rep("grain_yield", 4),
                         unit = c("kg/ha", rep("t/ha", 3)),  # 1000x wrong on o1
                         stringsAsFactors = FALSE))
  b <- pair$baseline;     b@obs_schema <- os_t
  i <- pair$intervention; i@obs_schema <- os_kg
  expect_error(kernR:::.validate_manifest_pair(b, i), "disagrees on output")
  # Matching schemas pass.
  i2 <- pair$intervention; i2@obs_schema <- os_t
  expect_silent(kernR:::.validate_manifest_pair(b, i2))
})

test_that(".validate_manifest_pair is a graceful no-op without obs_schema/apsim", {
  skip_unless_pesto()
  pair <- .make_scenario_pair(intervention_shift = 0.3, seed = 9L)
  expect_silent(kernR:::.validate_manifest_pair(pair$baseline,
                                                pair$intervention))
})
