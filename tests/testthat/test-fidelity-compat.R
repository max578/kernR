# Fidelity-provenance awareness when consuming PESTO manifests: the
# soft compatibility check in dr_date_scenario(), its escalation under
# strict_fidelity, and provenance threading into the dr_date_scenario
# and mmd_ppc results.

skip_if_not_installed("PESTO")

# Build a minimal pesto_ensemble_manifest with an optional fidelity
# provenance record. Schema is whatever the two-scenario test needs:
# shared parameter + output columns, equal row counts.
mk_manifest <- function(run_id, shift = 0, fidelity = NULL,
                        n = 24L, seed = 1L) {
  set.seed(seed)
  params <- data.frame(
    real_name = paste0("r", seq_len(n)),
    p1 = stats::rnorm(n), p2 = stats::rnorm(n)
  )
  outputs <- data.frame(
    real_name = paste0("r", seq_len(n)),
    o1 = stats::rnorm(n) + shift, o2 = stats::rnorm(n) + shift
  )
  PESTO::pesto_ensemble_manifest(
    run_id = run_id, params = params, outputs = outputs,
    weights = c(o1 = 1, o2 = 1), obs_target = c(o1 = 0, o2 = 0),
    data_hash = "sha256:test", pesto_version = "0.4.1",
    timestamp = Sys.time(), method = "ies_callback",
    noptmax = 3L, lambda_schedule = c(1, 1, 1), fidelity = fidelity
  )
}

# A multi-fidelity provenance record of a given shape.
mf_record <- function(final_level = 1L, n_levels = 2L) {
  list(type = "multifidelity",
       schedule = c(0L, final_level),
       final_level = final_level, n_levels = n_levels,
       costs = c(1, 9))
}

test_that("matched / both-single-fidelity provenance passes silently", {
  expect_silent(
    .validate_fidelity_pair(mk_manifest("b"), mk_manifest("i"))
  )
  expect_silent(
    .validate_fidelity_pair(mk_manifest("b", fidelity = mf_record()),
                            mk_manifest("i", fidelity = mf_record()))
  )
})

test_that("single vs multi-fidelity mismatch warns, errors under strict", {
  b <- mk_manifest("b")
  i <- mk_manifest("i", fidelity = mf_record())
  expect_warning(.validate_fidelity_pair(b, i), "single-fidelity")
  expect_error(
    .validate_fidelity_pair(b, i, strict = TRUE), "single-fidelity"
  )
})

test_that("differing multi-fidelity stack shape warns / errors", {
  b <- mk_manifest("b", fidelity = mf_record(final_level = 1L, n_levels = 2L))
  i <- mk_manifest("i", fidelity = mf_record(final_level = 2L, n_levels = 3L))
  expect_warning(.validate_fidelity_pair(b, i), "final_level|n_levels")
  expect_error(.validate_fidelity_pair(b, i, strict = TRUE))
})

test_that(".fidelity_mismatch_message returns NULL for compatible records", {
  expect_null(.fidelity_mismatch_message(NULL, NULL))
  expect_null(.fidelity_mismatch_message(mf_record(), mf_record()))
  expect_false(is.null(.fidelity_mismatch_message(NULL, mf_record())))
})

test_that("dr_date_scenario threads fidelity provenance into the result", {
  b <- mk_manifest("baseline",     fidelity = mf_record())
  i <- mk_manifest("intervention", shift = 0.6, fidelity = mf_record())
  res <- dr_date_scenario(b, i, n_permutations = 40L, seed = 1L)
  expect_false(is.null(res$fidelity))
  expect_identical(res$fidelity$baseline$type, "multifidelity")
  expect_identical(res$fidelity$intervention$final_level, 1L)
  # print() must not error and should mention fidelity.
  expect_output(print(res), "fidelity")
})

test_that("dr_date_scenario warns on mismatched fidelity, stops if strict", {
  b <- mk_manifest("baseline")                                   # single
  i <- mk_manifest("intervention", shift = 0.6, fidelity = mf_record())
  expect_warning(dr_date_scenario(b, i, n_permutations = 20L, seed = 1L),
                 "fidelity")
  expect_error(
    dr_date_scenario(b, i, n_permutations = 20L, seed = 1L,
                     strict_fidelity = TRUE),
    "fidelity"
  )
})

test_that("mmd_ppc records fidelity provenance in pesto_metadata", {
  m <- mk_manifest("m", fidelity = mf_record())
  observed <- matrix(stats::rnorm(12L), ncol = 2L)   # 6 rows, 2 outputs
  out <- mmd_ppc(m, observed = observed)
  expect_identical(out$pesto_metadata$fidelity$type, "multifidelity")
  # Single-fidelity manifest records NULL provenance.
  m0 <- mk_manifest("m0")
  out0 <- mmd_ppc(m0, observed = observed)
  expect_null(out0$pesto_metadata$fidelity)
})
