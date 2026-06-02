# coverage.R -- graded coverage / calibration diagnostic for ensembles.
#
# Where mmd_ppc() and ksd_test() return a binary reject/accept verdict on
# whether an ensemble matches observations, coverage_test() reports *how* and
# *which way* an ensemble is mis-calibrated. It computes the probability
# integral transform (PIT) of each held-out observation within the predictive
# ensemble, then summarises: empirical coverage at nominal interval levels, a
# signed dispersion ratio (Var(PIT) vs the uniform 1/12 -- above one is
# under-dispersion / over-confidence, below one is over-dispersion), a bias
# indicator (mean PIT vs 0.5), and a rank-histogram uniformity test. This is
# the interpretable companion to the kernel tests, built for the validation
# role: it turns "the ensemble is mis-calibrated" into "your 90% intervals
# cover 62% -- under-dispersed by a factor of 1.6".

#' Coverage / Calibration Diagnostic for a Predictive Ensemble
#'
#' Quantifies *how* a predictive ensemble is calibrated against held-out
#' observations, rather than only testing whether it differs from them. Each
#' observation is mapped to its probability integral transform (PIT) within the
#' ensemble; calibrated draws give uniform PITs. The result reports empirical
#' coverage at nominal interval levels, a signed dispersion ratio, a bias
#' indicator, and a rank-histogram uniformity test, and classifies the ensemble
#' as calibrated, under-dispersed (over-confident), over-dispersed, or biased.
#'
#' This is the graded complement to the binary kernel verdicts [mmd_ppc()] and
#' [ksd_test()]. Its motivating use is ensemble over-confidence: an
#' under-dispersed ensemble (for example a collapsed iterative-ensemble-smoother
#' posterior) gives a U-shaped rank histogram, a dispersion ratio above one, and
#' empirical coverage below nominal -- all of which this function names and
#' measures.
#'
#' Calibration is assessed **per output dimension and pooled across dimensions**
#' (marginal calibration); it does not test the joint dependence structure, for
#' which the two-sample [mmd_ppc()] is the right tool. With few held-out
#' observations the uniformity test has low power; read the coverage and
#' dispersion summaries (which are informative at any sample size) alongside it.
#'
#' @param x Numeric matrix `n_draws x d` of predictive draws, a
#'   `pesto_ensemble` (see [pesto_ensemble()]), or a `pesto_ensemble_manifest`.
#' @param observed Numeric matrix `n_obs x d` of held-out observations. When `x`
#'   is a `pesto_ensemble` carrying an `observed` slot, may be `NULL`.
#' @param levels Numeric vector in `(0, 1)`. Central predictive-interval levels
#'   at which to report empirical coverage. Default `c(0.5, 0.8, 0.9)`.
#' @param n_bins Integer. Number of equal-width bins for the rank-histogram
#'   uniformity test. Default `10`.
#' @param alpha Numeric in `(0, 1)`. Significance level for the calibration
#'   verdict. Default `0.05`.
#'
#' @returns An object of class `"coverage_test"` with components:
#'   \describe{
#'     \item{coverage}{Data frame of `nominal` vs pooled `empirical` coverage at
#'       each requested level.}
#'     \item{coverage_by_dim}{Matrix of empirical coverage per dimension x
#'       level.}
#'     \item{dispersion_ratio}{`Var(PIT) / (1/12)`: above one under-dispersed,
#'       below one over-dispersed.}
#'     \item{mean_pit}{Pooled mean PIT (0.5 under calibration; a bias
#'       indicator).}
#'     \item{calibration}{List with the rank-histogram chi-squared `statistic`,
#'       `p_value`, and `n_bins`.}
#'     \item{verdict}{Character: the calibration classification.}
#'     \item{reject}{Logical: `calibration$p_value <= alpha`.}
#'     \item{pit}{The `n_obs x d` matrix of PIT values.}
#'     \item{n_draws, n_obs, dimension, levels, alpha}{Inputs / sizes.}
#'     \item{pesto_metadata}{Provenance carried from a `pesto_ensemble_manifest`
#'       (including `fidelity`), or `NULL`.}
#'   }
#'
#' @references
#' Gneiting, T., Balabdaoui, F., & Raftery, A. E. (2007). Probabilistic
#' forecasts, calibration and sharpness. *Journal of the Royal Statistical
#' Society B*, 69(2), 243-268.
#'
#' Hamill, T. M. (2001). Interpretation of rank histograms for verifying
#' ensemble forecasts. *Monthly Weather Review*, 129(3), 550-560.
#'
#' @examples
#' set.seed(1)
#' obs <- matrix(stats::rnorm(120L), ncol = 2L)
#'
#' # Calibrated: predictive ensemble from the same law
#' ens_ok <- matrix(stats::rnorm(2000L), ncol = 2L)
#' coverage_test(ens_ok, obs)
#'
#' # Under-dispersed (over-confident): ensemble too narrow
#' ens_tight <- matrix(stats::rnorm(2000L, sd = 0.4), ncol = 2L)
#' coverage_test(ens_tight, obs)
#'
#' @seealso [mmd_ppc()], [ksd_test()], [pesto_ensemble()]
#' @family goodness-of-fit tests
#' @author Max Moldovan, \email{max.moldovan@@adelaide.edu.au}
#' @export
coverage_test <- function(x, ...) UseMethod("coverage_test")

#' @rdname coverage_test
#' @export
coverage_test.default <- function(x, observed,
                                  levels = c(0.5, 0.8, 0.9),
                                  n_bins = 10L,
                                  alpha = 0.05,
                                  ...) {
  cl <- match.call()
  ensemble <- as.matrix(x)
  if (missing(observed) || is.null(observed)) {
    stop("`observed` must be supplied when `x` is a matrix.", call. = FALSE)
  }
  observed <- as.matrix(observed)
  .check_coverage_input(ensemble, observed, levels, n_bins, alpha)

  n_draws <- nrow(ensemble)
  n_obs   <- nrow(observed)
  d       <- ncol(ensemble)

  # Probability integral transform (mid-rank, tie-safe) ------------------
  pit <- .ensemble_pit(ensemble, observed)

  # Coverage at each nominal central-interval level ----------------------
  cov_dim <- vapply(levels, function(l) {
    lo <- (1 - l) / 2
    hi <- (1 + l) / 2
    colMeans(pit >= lo & pit <= hi)
  }, numeric(d))
  cov_dim <- matrix(cov_dim, nrow = d,
                    dimnames = list(paste0("dim", seq_len(d)),
                                    paste0(round(100 * levels), "%")))
  empirical <- vapply(levels, function(l) {
    lo <- (1 - l) / 2
    hi <- (1 + l) / 2
    mean(pit >= lo & pit <= hi)
  }, numeric(1L))
  coverage <- data.frame(nominal = levels, empirical = empirical)

  # Dispersion + bias + rank-histogram uniformity ------------------------
  pit_vec <- as.numeric(pit)
  dispersion_ratio <- stats::var(pit_vec) / (1 / 12)
  mean_pit <- mean(pit_vec)
  calib <- .pit_uniformity(pit_vec, n_bins)

  verdict <- .coverage_verdict(calib$p_value, alpha, dispersion_ratio,
                               mean_pit)

  structure(
    list(
      coverage         = coverage,
      coverage_by_dim  = t(cov_dim),
      dispersion_ratio = dispersion_ratio,
      mean_pit         = mean_pit,
      calibration      = calib,
      verdict          = verdict,
      reject           = calib$p_value <= alpha,
      pit              = pit,
      n_draws          = n_draws,
      n_obs            = n_obs,
      dimension        = d,
      levels           = levels,
      alpha            = alpha,
      pesto_metadata   = NULL,
      call             = cl
    ),
    class = "coverage_test"
  )
}

#' @rdname coverage_test
#' @export
coverage_test.pesto_ensemble <- function(x, observed = NULL, ...) {
  if (is.null(observed)) observed <- x$observed
  if (is.null(observed)) {
    stop("`observed` must be supplied: ensemble carries no `observed` slot.",
         call. = FALSE)
  }
  out <- coverage_test.default(x$posterior, observed = observed, ...)
  out$pesto_metadata <- x$metadata
  out$call <- match.call()
  out
}

#' @rdname coverage_test
#' @param outputs Optional character vector of manifest output columns to test
#'   (the `real_name` column is always excluded). Defaults to all numeric output
#'   columns. Used only for the `pesto_ensemble_manifest` method.
#' @export
coverage_test.pesto_ensemble_manifest <- function(x, observed,
                                                  outputs = NULL, ...) {
  cl <- match.call()
  if (missing(observed) || is.null(observed)) {
    stop("`observed` must be supplied for the manifest method (a held-out ",
         "matrix with one column per tested output).", call. = FALSE)
  }
  outputs_df <- as.data.frame(x@outputs)
  out_cols   <- setdiff(names(outputs_df), "real_name")
  if (is.null(outputs)) {
    outputs <- out_cols
  } else {
    outputs <- as.character(outputs)
    missing_out <- setdiff(outputs, out_cols)
    if (length(missing_out)) {
      stop("`outputs` columns not found in manifest outputs: ",
           paste(missing_out, collapse = ", "), call. = FALSE)
    }
  }
  ensemble <- as.matrix(outputs_df[, outputs, drop = FALSE])

  out <- coverage_test.default(ensemble, observed = observed, ...)
  out$pesto_metadata <- list(
    run_id        = x@run_id,
    pesto_version = x@pesto_version,
    method        = x@method,
    outputs_used  = outputs,
    fidelity      = x@fidelity
  )
  out$call <- cl
  out
}

#' @export
print.coverage_test <- function(x, ...) {
  cat("\n  Coverage / calibration diagnostic\n\n")
  cat("Ensemble:   ", x$n_draws, " draws x ", x$dimension, " dims; ",
      x$n_obs, " held-out obs\n", sep = "")
  cov_str <- paste(
    sprintf("%g%%->%.0f%%", 100 * x$coverage$nominal,
            100 * x$coverage$empirical),
    collapse = "   "
  )
  cat("Coverage:   ", cov_str, "  (nominal -> empirical)\n", sep = "")
  cat("Dispersion: ratio = ",
      formatC(x$dispersion_ratio, digits = 3L, format = "f"),
      " (Var(PIT)/(1/12); >1 under-, <1 over-dispersed)\n", sep = "")
  cat("Mean PIT:   ", formatC(x$mean_pit, digits = 3L, format = "f"),
      " (0.5 = unbiased)\n", sep = "")
  cat("Calibration: chi2 = ",
      formatC(x$calibration$statistic, digits = 3L, format = "g"),
      ", p = ", formatC(x$calibration$p_value, digits = 4L, format = "f"),
      " (rank histogram, ", x$calibration$n_bins, " bins)\n", sep = "")
  cat("Verdict:    ", x$verdict, "\n", sep = "")
  if (length(x$pesto_metadata)) {
    cat("Metadata:   ", paste(names(x$pesto_metadata), collapse = ", "),
        "\n", sep = "")
  }
  cat("\n")
  invisible(x)
}

# Internal helpers -----------------------------------------------------------

#' Mid-rank probability integral transform of observations within an ensemble
#'
#' @param ensemble Numeric `n_draws x d` matrix.
#' @param observed Numeric `n_obs x d` matrix.
#' @returns `n_obs x d` matrix of PIT values, tie-safe (mid-distribution).
#' @noRd
#' @keywords internal
.ensemble_pit <- function(ensemble, observed) {
  d <- ncol(ensemble)
  n_obs <- nrow(observed)
  pit <- matrix(0, nrow = n_obs, ncol = d)
  for (j in seq_len(d)) {
    ej <- ensemble[, j]
    pit[, j] <- vapply(observed[, j], function(o) {
      (mean(ej < o) + mean(ej <= o)) / 2
    }, numeric(1L))
  }
  pit
}

#' Rank-histogram chi-squared test of PIT uniformity
#'
#' @param pit_vec Numeric vector of PIT values in `[0, 1]`.
#' @param n_bins Integer bin count.
#' @returns List with `statistic`, `p_value`, `n_bins`.
#' @noRd
#' @keywords internal
.pit_uniformity <- function(pit_vec, n_bins) {
  breaks <- seq(0, 1, length.out = n_bins + 1L)
  counts <- tabulate(
    findInterval(pit_vec, breaks, rightmost.closed = TRUE, all.inside = TRUE),
    nbins = n_bins
  )
  expected <- length(pit_vec) / n_bins
  stat <- sum((counts - expected)^2 / expected)
  p_value <- stats::pchisq(stat, df = n_bins - 1L, lower.tail = FALSE)
  list(statistic = stat, p_value = p_value, n_bins = as.integer(n_bins))
}

#' Classify a calibration result
#'
#' @param p_value Rank-histogram uniformity p-value.
#' @param alpha Verdict level.
#' @param ratio Dispersion ratio.
#' @param mean_pit Pooled mean PIT.
#' @returns A single human-readable verdict string.
#' @noRd
#' @keywords internal
.coverage_verdict <- function(p_value, alpha, ratio, mean_pit) {
  if (p_value > alpha) {
    return("calibrated (PIT consistent with uniform)")
  }
  issues <- character(0)
  if (abs(mean_pit - 0.5) > 0.05) {
    issues <- c(issues,
                if (mean_pit < 0.5) "biased (ensemble too high)"
                else "biased (ensemble too low)")
  }
  if (ratio > 1.1) {
    issues <- c(issues, "under-dispersed (over-confident)")
  } else if (ratio < 0.9) {
    issues <- c(issues, "over-dispersed")
  }
  if (length(issues) == 0L) {
    issues <- "mis-calibrated (non-uniform PIT; shape not mean/variance)"
  }
  paste0("REJECT -- ", paste(issues, collapse = "; "))
}

#' Validate coverage_test() inputs
#'
#' @param ensemble Coerced ensemble matrix.
#' @param observed Coerced observation matrix.
#' @param levels Requested interval levels.
#' @param n_bins Bin count.
#' @param alpha Verdict level.
#' @returns `invisible(NULL)`; called for its error side effects.
#' @noRd
#' @keywords internal
.check_coverage_input <- function(ensemble, observed, levels, n_bins, alpha) {
  if (!is.numeric(ensemble) || any(!is.finite(ensemble))) {
    stop("`x` must be numeric and contain only finite values.", call. = FALSE)
  }
  if (!is.numeric(observed) || any(!is.finite(observed))) {
    stop("`observed` must be numeric and contain only finite values.",
         call. = FALSE)
  }
  if (ncol(ensemble) != ncol(observed)) {
    stop("`x` and `observed` must have the same number of columns.",
         call. = FALSE)
  }
  if (nrow(ensemble) < 10L) {
    stop("`x` (the ensemble) must have at least 10 draws.", call. = FALSE)
  }
  if (nrow(observed) < 1L) {
    stop("`observed` must have at least one row.", call. = FALSE)
  }
  if (!is.numeric(levels) || any(levels <= 0) || any(levels >= 1)) {
    stop("`levels` must be numbers in (0, 1).", call. = FALSE)
  }
  if (!is.numeric(n_bins) || length(n_bins) != 1L || n_bins < 2) {
    stop("`n_bins` must be a single integer >= 2.", call. = FALSE)
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a single number in (0, 1).", call. = FALSE)
  }
  invisible(NULL)
}
