# test-manifest-contract.R -- conformance of kernR's emitted manifest to the
# federation's REFERENCE contract (K4).
#
# The oracle for every expectation below is the reference implementation
# `integration/orchestra_manifest.R` in the ORCHESTRA coordination workspace,
# which kernR did not author: the class identity it dispatches on, the payload
# hash recipe its `verify_manifest()` recomputes, and the property set its
# validator reads. Emitting a same-shaped object is not conformance -- a
# consumer that refuses the object is the failure mode these tests exist to
# catch, and it is invisible to a check that uses kernR's own
# `verify_manifest()` on kernR's own manifest.
#
# The final block runs the real reference file when the environment variable
# `KERNR_ORCHESTRA_CONTRACT` points at it (a coordination-workspace path that
# does not exist on a public runner); the blocks before it encode the same
# contract requirements so the gate still bites on CI.

# The reference class is declared outside any package namespace
# (`S7::new_class("orchestra_manifest", ...)`, reference contract lines
# 99-122), so its S7 identity is the bare name. A class declared under a
# package carries "<pkg>::<name>" instead and is a foreign look-alike that no
# orchestra consumer recognises.
test_that("the emitted manifest carries the federation's bare class identity", {
  set.seed(11L)
  x <- matrix(stats::rnorm(80L), ncol = 1L)
  y <- matrix(stats::rnorm(80L), ncol = 1L)
  m <- as_orchestra_manifest(hsic_test(x, y, n_permutations = 60L))

  expect_identical(class(m)[[1L]], "orchestra_manifest")
  expect_null(attr(S7::S7_class(m), "package"))
})

# The reference property set and its order, read off the reference contract's
# `properties = list(...)` (lines 101-122). A consumer reads these slots by
# name; a missing or renamed slot breaks it silently.
test_that("the emitted manifest carries the reference property set", {
  set.seed(12L)
  x <- matrix(stats::rnorm(80L), ncol = 1L)
  y <- matrix(stats::rnorm(80L), ncol = 1L)
  m <- as_orchestra_manifest(hsic_test(x, y, n_permutations = 60L))

  reference_properties <- c(
    "manifest_version", "emitter_package", "emitter_version",
    "inferential_target", "run_id", "method", "seed", "params", "outputs",
    "weights", "obs_target", "obs_schema", "summary", "consumed_manifests",
    "metadata", "timestamp", "data_hash"
  )
  expect_identical(S7::prop_names(m), reference_properties)
  expect_identical(m@manifest_version, "2.0.0-draft")
  expect_identical(m@emitter_package, "kernR")
  expect_identical(m@emitter_version,
                   as.character(utils::packageVersion("kernR")))
})

# The reference hash recipe (contract lines 59-76): sha256 over
# `serialize(list(params, outputs, weights, obs_target, seed, summary),
# version = 2L)` with the fixed 14-byte header dropped, prefixed "sha256:".
# Recomputed here from that description alone -- not by calling kernR's own
# hasher, which would only prove self-consistency.
test_that("data_hash follows the reference payload-hash recipe", {
  reference_hash <- function(m) {
    obj <- list(m@params, m@outputs, m@weights, m@obs_target, m@seed,
                m@summary)
    raw <- serialize(obj, connection = NULL, version = 2L)
    raw <- raw[-seq_len(14L)]
    paste0("sha256:", digest::digest(raw, algo = "sha256", serialize = FALSE))
  }

  set.seed(13L)
  x <- matrix(stats::rnorm(80L), ncol = 1L)
  y <- matrix(stats::rnorm(80L), ncol = 1L)
  m_kernel <- as_orchestra_manifest(hsic_test(x, y, n_permutations = 60L))
  expect_identical(m_kernel@data_hash, reference_hash(m_kernel))

  n <- 60L
  X <- matrix(stats::rnorm(n), ncol = 1L)
  tt <- stats::rbinom(n, 1L, 0.5)
  mech <- function(theta, X, t) theta[1] + theta[2] * t
  yy <- mech(c(1, 2), X, tt) + stats::rnorm(n, sd = 0.3)
  theta <- cbind(stats::rnorm(200L, 1, 0.05), stats::rnorm(200L, 2, 0.3))
  m_taci <- as_orchestra_manifest(
    taci_test(posterior = theta, mechanism = mech, X = X, treatment = tt,
              outcome = yy, treatment_type = "binary", n_perm = 100L,
              seed = 13L))
  expect_identical(m_taci@data_hash, reference_hash(m_taci))
})

# The end-to-end oracle: the reference file itself, sourced and applied.
test_that("the reference contract accepts a kernR manifest (KERNR_ORCHESTRA_CONTRACT)", {
  ref_path <- Sys.getenv("KERNR_ORCHESTRA_CONTRACT", unset = "")
  skip_if(!nzchar(ref_path) || !file.exists(ref_path),
          "reference orchestra contract not available")

  ref <- new.env(parent = globalenv())
  sys.source(ref_path, envir = ref)

  n <- 60L
  set.seed(14L)
  X <- matrix(stats::rnorm(n), ncol = 1L)
  tt <- stats::rbinom(n, 1L, 0.5)
  mech <- function(theta, X, t) theta[1] + theta[2] * t
  yy <- mech(c(1, 2), X, tt) + stats::rnorm(n, sd = 0.3)

  theta_ok <- cbind(stats::rnorm(200L, 1, 0.05), stats::rnorm(200L, 2, 0.3))
  m_ok <- as_orchestra_manifest(
    taci_test(posterior = theta_ok, mechanism = mech, X = X, treatment = tt,
              outcome = yy, treatment_type = "binary", n_perm = 100L,
              seed = 14L))

  expect_true(S7::S7_inherits(m_ok, ref$orchestra_manifest))
  expect_true(ref$verify_manifest(m_ok)$ok)
  expect_no_error(ref$consume_manifest(m_ok))

  # A verdict whose posterior-adequacy gate failed must reach the consumer as
  # an abstention, not as a clean decision.
  theta_bad <- cbind(stats::rnorm(200L, 1, 1e-6), stats::rnorm(200L, 2, 1e-6))
  m_abs <- suppressWarnings(as_orchestra_manifest(
    taci_test(posterior = theta_bad, mechanism = mech, X = X, treatment = tt,
              outcome = yy, treatment_type = "binary", n_perm = 50L,
              seed = 14L)))
  expect_no_error(ref$consume_manifest(m_abs))
  expect_true(m_abs@summary$abstained)
  expect_false(is.na(m_abs@summary$abstain_reason))
})
