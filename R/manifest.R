# manifest.R -- the orchestra contract-emitting result object.
#
# kernR consumes the C2 manifest contract (`mmd_ppc()` / `coverage_test()` /
# `joint_coverage_test()` on a `pesto_ensemble_manifest`) but, until now,
# emitted nothing: a TACI mechanism-consistency verdict or a kernel hypothesis
# test had no typed edge into `decideR`. This file carries a stand-in
# `orchestra_manifest` class whose property names and integrity hash match the
# federation's reference implementation
# (`ORCHESTRA_dev/integration/orchestra_manifest.R`) field for field, so a
# consumer reads a kernR manifest with no special casing. The contract is
# defined here (rather than depended on) because the reference implementation
# lives in the composition layer, not an installable package -- the same
# choice `gpfield` and `optimix` made.
#
# kernR's verdicts are causal / distributional-treatment-effect statements
# (TACI's mechanism-consistency decision; a kernel test's reject/accept call
# on an effect), so the emit target is `"treatment_effects"` throughout. The
# reliability story (K2/K4 elsewhere in this audit pass) rides in the typed
# `summary`: `abstained` is set from exactly the same reliability flags a
# caller would otherwise have to read off the raw result by hand
# (`ess_warning`, `density_ratio_warning`, `posterior_adequacy$ok`).

# The manifest schema version this member emits. Tracks the reference
# implementation; bump only with an additive (x.y) or breaking (x.0) change
# there.
MANIFEST_VERSION <- "2.0.0-draft"

# The inferential-target enum, kept identical to the reference contract so a
# consumer's dispatch never sees an unknown token.
.INFERENTIAL_TARGETS <- c("parameters", "predictions",
                          "treatment_effects", "decisions",
                          "breeding_values", "marker_associations",
                          "structure")

#' Null-coalescing operator
#'
#' @param a,b The candidate and the fallback.
#'
#' @returns `a` unless it is `NULL`, in which case `b`.
#' @noRd
#' @keywords internal
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Payload integrity hash matching the orchestra contract
#'
#' SHA-256 over the load-bearing data slots only (the schema, version and
#' timestamp are metadata and are not hashed -- spec section 8.4). The
#' `"sha256:"` prefix and the slot ordering match the reference implementation
#' so a manifest hashed by kernR verifies under the federation's
#' `verify_manifest()`. Uses `tools::sha256sum()` (base R) so the emitter adds
#' no new hard dependency.
#'
#' @param params,outputs,weights,obs_target,seed,summary The load-bearing
#'   slots.
#'
#' @returns A single `"sha256:"`-prefixed hex string.
#' @noRd
#' @keywords internal
.manifest_hash_payload <- function(params, outputs, weights, obs_target,
                                   seed, summary = NULL) {
  obj <- list(params, outputs, weights, obs_target, seed, summary)
  raw <- serialize(obj, connection = NULL)
  paste0("sha256:", as.character(tools::sha256sum(bytes = raw)))
}

#' A typed home for a kernR verdict (contract v1.1)
#'
#' Mirrors the reference contract's `manifest_summary()`: the result's
#' headline decision, an abstain flag (set when a reliability gate --
#' insufficient ESS, a poor density-ratio fit, an inadequate TACI posterior --
#' failed), and a list of key metrics. Integrity-hashed with the rest of the
#' payload.
#'
#' @param headline A single non-empty status string.
#' @param abstained A single logical; `TRUE` when the verdict is not
#'   trustworthy (a reliability gate failed).
#' @param metrics A named list of scalar metrics.
#' @param abstain_reason A single string, or `NA`.
#'
#' @returns A classed `manifest_summary` list.
#' @noRd
#' @keywords internal
manifest_summary <- function(headline, abstained = FALSE, metrics = list(),
                             abstain_reason = NA_character_) {
  stopifnot(is.character(headline), length(headline) == 1L, nzchar(headline),
            is.logical(abstained), length(abstained) == 1L, is.list(metrics))
  structure(
    list(headline = headline, abstained = isTRUE(abstained),
         abstain_reason = abstain_reason, metrics = metrics),
    class = "manifest_summary"
  )
}

# The contract class -----------------------------------------------------

#' The orchestra ensemble-manifest contract (kernR-side implementation)
#'
#' A versioned, hashed, provenance-complete S7 result object,
#' property-compatible with the federation's reference `orchestra_manifest`.
#' A `taci_result` or a `kernel_test_result` is adapted into one through
#' [as_orchestra_manifest()]: the verdict rides in the typed `summary`
#' (`inferential_target = "treatment_effects"`), the test provenance rides in
#' `metadata`, and the payload is integrity-hashed so a tampered verdict is
#' detected downstream.
#'
#' @usage NULL
#'
#' @returns An S7 object of class `orchestra_manifest`.
#'
#' @seealso [as_orchestra_manifest()], [verify_manifest()]
#' @family orchestra manifest
#' @export
orchestra_manifest <- S7::new_class(
  "orchestra_manifest",
  package = "kernR",
  properties = list(
    manifest_version   = S7::new_property(S7::class_character,
                                          default = MANIFEST_VERSION),
    emitter_package    = S7::class_character,
    emitter_version    = S7::class_character,
    inferential_target = S7::class_character,
    run_id             = S7::class_character,
    method             = S7::class_character,
    seed               = S7::new_property(S7::class_integer,
                                          default = NA_integer_),
    params             = S7::new_property(S7::class_data.frame,
                                          default = data.frame()),
    outputs            = S7::new_property(S7::class_any, default = NULL),
    weights            = S7::new_property(S7::class_any, default = NULL),
    obs_target         = S7::new_property(S7::class_any, default = NULL),
    obs_schema         = S7::new_property(S7::class_any, default = NULL),
    summary            = S7::new_property(S7::class_any, default = NULL),
    consumed_manifests = S7::new_property(S7::class_list, default = list()),
    metadata           = S7::new_property(S7::class_list, default = list()),
    timestamp          = S7::class_POSIXct,
    data_hash          = S7::class_character
  ),
  validator = function(self) {
    errs <- character(0)
    if (length(self@run_id) != 1L || !nzchar(self@run_id)) {
      errs <- c(errs, "`run_id` must be a single non-empty string")
    }
    if (length(self@inferential_target) != 1L ||
        !self@inferential_target %in% .INFERENTIAL_TARGETS) {
      errs <- c(errs, sprintf("`inferential_target` must be one of %s",
                              paste(.INFERENTIAL_TARGETS, collapse = ", ")))
    }
    if (length(self@data_hash) != 1L) {
      errs <- c(errs, "`data_hash` must be a single string")
    }
    if (!is.null(self@summary) &&
        !inherits(self@summary, "manifest_summary")) {
      errs <- c(errs, "`summary`, when set, must be a manifest_summary() object")
    }
    if (length(errs) == 0L) NULL else paste(errs, collapse = "; ")
  }
)

#' Verify an orchestra manifest's payload integrity
#'
#' Recomputes the payload hash from the object's own load-bearing slots and
#' compares it to the stored `data_hash`. A mismatch means the payload was
#' modified after emission. Matches the reference contract's
#' `verify_manifest()`.
#'
#' @param m An `orchestra_manifest` object.
#'
#' @returns A list with logical `ok` and a human-readable `message`.
#'
#' @examples
#' mechanism <- function(theta, X, t) theta[3] + theta[1] * (1 - exp(-theta[2] * t))
#' set.seed(1L)
#' n <- 60L
#' nrate <- runif(n, 0, 200)
#' yield <- 1.1 + 4.2 * (1 - exp(-0.018 * nrate)) + rnorm(n, 0, 0.25)
#' post <- cbind(ymax = rnorm(100L, 4.2, 0.30),
#'              rate = rnorm(100L, 0.018, 0.004),
#'              y0   = rnorm(100L, 1.1, 0.15))
#' verdict <- taci_test(post, mechanism, X = matrix(1, n, 1),
#'                      treatment = nrate, outcome = yield,
#'                      n_perm = 99L, seed = 1L)
#' m <- as_orchestra_manifest(verdict)
#' verify_manifest(m)$ok  # TRUE -- payload untampered
#'
#' @seealso [as_orchestra_manifest()]
#' @family orchestra manifest
#' @export
verify_manifest <- function(m) {
  if (!S7::S7_inherits(m, orchestra_manifest)) {
    return(list(ok = FALSE, message = "not an orchestra_manifest"))
  }
  recomputed <- .manifest_hash_payload(m@params, m@outputs, m@weights,
                                       m@obs_target, m@seed, m@summary)
  ok <- identical(recomputed, m@data_hash)
  list(ok = ok,
       message = if (ok) "data_hash verified"
                 else "data_hash MISMATCH -- manifest payload was modified")
}

# The emit generic --------------------------------------------------------

#' Emit a kernR verdict as an orchestra manifest
#'
#' The federation's emit / adapt generic. Both methods set
#' `inferential_target = "treatment_effects"`: kernR's verdicts are
#' hypothesis-test or mechanism-consistency calls about a treatment effect,
#' never a point estimate or a prediction. `abstained` in the typed `summary`
#' is set from the result's own reliability flags -- `posterior_adequacy$ok`
#' for a TACI verdict, `ess_warning` / `density_ratio_warning` for a kernel
#' test -- so a downstream consumer reading only the manifest sees the same
#' reliability signal the raw result already carried, without having to know
#' which flag belongs to which verb.
#'
#' @param x A `taci_result` or `kernel_test_result` object.
#' @param ... Method-specific arguments.
#' @param run_id Optional character run identifier; a content hash is derived
#'   when omitted, so two identical verdicts get the same identifier.
#'
#' @returns An `orchestra_manifest` S7 object.
#'
#' @examples
#' mechanism <- function(theta, X, t) theta[3] + theta[1] * (1 - exp(-theta[2] * t))
#' set.seed(1L)
#' n <- 60L
#' nrate <- runif(n, 0, 200)
#' yield <- 1.1 + 4.2 * (1 - exp(-0.018 * nrate)) + rnorm(n, 0, 0.25)
#' post <- cbind(ymax = rnorm(100L, 4.2, 0.30),
#'              rate = rnorm(100L, 0.018, 0.004),
#'              y0   = rnorm(100L, 1.1, 0.15))
#' verdict <- taci_test(post, mechanism, X = matrix(1, n, 1),
#'                      treatment = nrate, outcome = yield,
#'                      n_perm = 99L, seed = 1L)
#' m <- as_orchestra_manifest(verdict)
#' m@summary$headline
#' m@summary$abstained  # FALSE: the posterior-adequacy gate passed
#'
#' @seealso [verify_manifest()], [taci_test()]
#' @family orchestra manifest
#' @export
as_orchestra_manifest <- function(x, ...) {
  UseMethod("as_orchestra_manifest")
}

#' @rdname as_orchestra_manifest
#' @export
as_orchestra_manifest.taci_result <- function(x, ..., run_id = NULL) {
  reliable <- isTRUE(x$posterior_adequacy$ok)

  summ <- manifest_summary(
    headline = x$verdict,
    abstained = !reliable,
    abstain_reason = if (reliable) NA_character_ else x$posterior_adequacy$reason,
    metrics = list(
      decision = x$decision,
      grounding = x$grounding,
      observed_statistic = x$observed_statistic,
      p_h0 = x$p_h0,
      h1_consistent = x$h1_consistent,
      h1_percentile = x$h1_percentile,
      borderline = x$borderline,
      dr_fit_quality = x$dr_fit_quality
    )
  )

  meta <- list(
    derived = FALSE,
    treatment_type = x$treatment_type,
    adjusted = x$adjusted,
    h0_mode = x$h0_mode,
    density_ratio = x$density_ratio,
    n = x$n,
    n_posterior = x$n_posterior,
    alpha = x$alpha,
    mechanism_provenance = x$mechanism_provenance,
    posterior_provenance = x$posterior_provenance
  )

  ev <- as.character(utils::packageVersion("kernR"))
  dh <- .manifest_hash_payload(data.frame(), NULL, NULL, NULL, NA_integer_,
                               summ)
  rid <- run_id %||% paste0("kernR-taci-",
                            substr(sub("^sha256:", "", dh), 1L, 12L))

  orchestra_manifest(
    manifest_version   = MANIFEST_VERSION,
    emitter_package    = "kernR",
    emitter_version    = ev,
    inferential_target = "treatment_effects",
    run_id             = rid,
    method             = "kernR:taci",
    seed               = NA_integer_,
    params             = data.frame(),
    summary            = summ,
    consumed_manifests = list(),
    metadata           = meta,
    timestamp          = Sys.time(),
    data_hash          = dh
  )
}

#' @rdname as_orchestra_manifest
#' @export
as_orchestra_manifest.kernel_test_result <- function(x, ..., run_id = NULL) {
  ess_warn <- isTRUE(x$ess_warning)
  dr_warn  <- isTRUE(x$density_ratio_warning)
  reliable <- !ess_warn && !dr_warn
  reason <- if (reliable) {
    NA_character_
  } else if (ess_warn && dr_warn) {
    "effective sample size below floor and a poor density-ratio fit"
  } else if (ess_warn) {
    "effective sample size below floor"
  } else {
    "poor density-ratio fit"
  }

  summ <- manifest_summary(
    headline = sprintf("%s: statistic = %.6g, p = %.4g", x$method,
                       x$statistic, x$p_value),
    abstained = !reliable,
    abstain_reason = reason,
    metrics = list(
      statistic = x$statistic,
      p_value = x$p_value,
      method = x$method,
      n = x$n,
      ess = x$ess,
      ess_warning = ess_warn,
      density_ratio_warning = dr_warn
    )
  )

  meta <- list(
    derived = FALSE,
    method = x$method,
    n_permutations = x$n_permutations,
    kernel_x = x$kernel_x,
    kernel_y = x$kernel_y
  )

  ev <- as.character(utils::packageVersion("kernR"))
  dh <- .manifest_hash_payload(data.frame(), NULL, x$weights, NULL,
                               NA_integer_, summ)
  rid <- run_id %||% paste0("kernR-", tolower(x$method), "-",
                            substr(sub("^sha256:", "", dh), 1L, 12L))

  orchestra_manifest(
    manifest_version   = MANIFEST_VERSION,
    emitter_package    = "kernR",
    emitter_version    = ev,
    inferential_target = "treatment_effects",
    run_id             = rid,
    method             = paste0("kernR:", x$method),
    seed               = NA_integer_,
    params             = data.frame(),
    summary            = summ,
    consumed_manifests = list(),
    metadata           = meta,
    timestamp          = Sys.time(),
    data_hash          = dh
  )
}
