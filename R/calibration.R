# calibration.R -- calibration of a point prediction against a held-out outcome.
#
# The sibling of coverage.R. `coverage_test()` asks whether a predictive
# *distribution* is the right width; this asks whether a predictive *level* is
# right, which is the quantity a decision layer prices. The four statistics are
# the continuous-outcome analogues of the clinical prediction-model suite
# (calibration slope, calibration-in-the-large, intercept, and the integrated
# calibration index), where the distinction between a model that ranks well and
# a model that is right about magnitude has been standard practice for years.
#
# The verdict is power-aware, which is the part that is not standard. A
# calibration test on a handful of observations cannot resolve a slope of 0.8
# from a slope of 1.0, and a test that reports "calibrated" in that situation is
# reporting the absence of power, not the presence of calibration. This suite
# refuses instead: when the interval is wider than the deviation the caller says
# they care about, the verdict is `"undetermined"` and the result carries the
# typed `kernR_abstention` marker the rest of the federation routes on.

#' Calibration suite for a continuous prediction
#'
#' Quantifies whether predictions are right in *level*, not merely in order.
#' A model can rank cases correctly while being systematically wrong about
#' magnitude, and it is magnitude that a loss-optimal decision prices.
#'
#' Four statistics, each with a bootstrap percentile interval:
#' \describe{
#'   \item{slope}{Coefficient of `observed ~ predicted`. Ideal 1. Below 1 means
#'     predictions are too spread out -- the classic over-fitting signature, and
#'     the quantity that collapses first when a model moves to a new site.}
#'   \item{intercept}{Intercept of `observed ~ offset(predicted)`. Ideal 0.}
#'   \item{citl}{Calibration-in-the-large, `mean(observed) - mean(predicted)`.
#'     Ideal 0. The pure level shift, free of any slope effect.}
#'   \item{ici}{Integrated calibration index: the mean absolute distance between
#'     a smooth of the observed outcome on the prediction and the prediction
#'     itself. In the units of the outcome.}
#' }
#'
#' @section The verdict is power-aware:
#' `verdict` is `"miscalibrated"` when the slope interval excludes 1;
#' `"calibrated"` when it contains 1 *and* is narrow enough to have resolved a
#' deviation of `tolerance`; and `"undetermined"` otherwise -- the test could
#' not have detected the miscalibration the caller cares about, so it declines
#' to certify calibration. An `"undetermined"` result is stamped
#' `kernR_abstention`, so a caller gating on calibration cannot mistake an
#' under-powered pass for a real one. `detectable` reports the smallest
#' deviation from a slope of 1 that this sample could resolve.
#'
#' @param x Numeric vector of predictions, or an object with a method.
#' @param observed Numeric vector of realised outcomes, the same length as `x`.
#' @param n_boot Bootstrap replicates for the percentile intervals.
#' @param tolerance The deviation from a slope of 1 the caller would want to
#'   detect. Drives the power-aware verdict; it does not affect the estimates.
#' @param conf_level Interval level for every reported statistic.
#' @param seed Optional integer seed for the bootstrap, for reproducibility.
#' @param ... Passed to methods.
#'
#' @returns An object of class `calibration_suite` (additionally
#'   `kernR_abstention` when the verdict is `"undetermined"`), a list with
#'   `slope`, `intercept`, `citl` and `ici` (each a named numeric of
#'   `estimate`, `lower`, `upper`), `curve` (a data.frame for plotting),
#'   `detectable`, `verdict`, `abstained`, and the inputs `n`, `n_boot`,
#'   `tolerance`, `conf_level`.
#'
#' @references
#' Van Calster, B., McLernon, D. J., van Smeden, M., Wynants, L., &
#' Steyerberg, E. W. (2019). Calibration: the Achilles heel of predictive
#' analytics. *BMC Medicine*, 17, 230.
#'
#' Austin, P. C., & Steyerberg, E. W. (2019). The Integrated Calibration Index
#' (ICI) and related metrics for quantifying the calibration of logistic
#' regression models. *Statistics in Medicine*, 38(21), 4051-4065.
#'
#' @examples
#' set.seed(1)
#' # A calibrated prediction is one the outcome varies around: generate the
#' # prediction first, then the outcome from it. Note that `truth + noise` is
#' # NOT calibrated as a prediction of `truth` -- it carries noise the outcome
#' # does not, so its slope attenuates below 1 by construction.
#' pred <- stats::rnorm(300L, mean = 5, sd = 2)
#' obs <- pred + stats::rnorm(300L, sd = 0.5)
#' calibration_suite(pred, obs, n_boot = 200L)
#'
#' # Over-spread predictions: slope below 1
#' calibration_suite(5 + 1.6 * (pred - 5), pred, n_boot = 200L)
#'
#' # A noisy handful cannot certify anything: verdict is "undetermined"
#' calibration_suite(pred[1:6], pred[1:6] + stats::rnorm(6L, sd = 3),
#'                   n_boot = 200L)
#'
#' @seealso [coverage_test()] for the width of a predictive distribution.
#' @family goodness-of-fit tests
#' @author Max Moldovan, \email{max.moldovan@@adelaide.edu.au}
#' @export
calibration_suite <- function(x, ...) UseMethod("calibration_suite")

#' @rdname calibration_suite
#' @export
calibration_suite.default <- function(x, observed,
                                      n_boot = 1000L,
                                      tolerance = 0.1,
                                      conf_level = 0.95,
                                      seed = NULL,
                                      ...) {
  cl <- match.call()
  predicted <- as.numeric(x)
  if (missing(observed) || is.null(observed)) {
    stop("`observed` must be supplied.", call. = FALSE)
  }
  observed <- as.numeric(observed)
  .check_calibration_input(predicted, observed, n_boot, tolerance, conf_level)

  n <- length(predicted)
  if (!is.null(seed)) set.seed(as.integer(seed))

  est <- .calibration_point(predicted, observed)

  # Bootstrap the pairs. Resamples that degenerate (a constant prediction
  # vector has no slope) are dropped and counted rather than silently
  # contributing an NA to a quantile.
  boot <- matrix(NA_real_, nrow = n_boot, ncol = 4L,
                 dimnames = list(NULL, c("slope", "intercept", "citl", "ici")))
  for (i1 in seq_len(n_boot)) {
    idx <- sample.int(n, n, replace = TRUE)
    b <- .calibration_point(predicted[idx], observed[idx])
    boot[i1, ] <- c(b$slope, b$intercept, b$citl, b$ici)
  }
  n_ok <- sum(stats::complete.cases(boot[, c("slope", "intercept", "citl")]))

  ci <- function(nm) {
    v <- boot[, nm]
    v <- v[is.finite(v)]
    if (!length(v)) {
      return(c(estimate = est[[nm]], lower = NA_real_, upper = NA_real_))
    }
    a <- (1 - conf_level) / 2
    q <- stats::quantile(v, probs = c(a, 1 - a), names = FALSE, na.rm = TRUE)
    c(estimate = est[[nm]], lower = q[1L], upper = q[2L])
  }
  slope <- ci("slope")
  intercept <- ci("intercept")
  citl <- ci("citl")
  ici <- ci("ici")

  # Power-aware verdict. `detectable` is the half-width of the slope interval:
  # the smallest departure from 1 this sample could have separated from 1.
  detectable <- unname((slope[["upper"]] - slope[["lower"]]) / 2)
  excludes_one <- is.finite(slope[["lower"]]) && is.finite(slope[["upper"]]) &&
    (slope[["lower"]] > 1 || slope[["upper"]] < 1)
  verdict <- if (!is.finite(detectable)) {
    "undetermined"
  } else if (excludes_one) {
    "miscalibrated"
  } else if (detectable > tolerance) {
    "undetermined"
  } else {
    "calibrated"
  }

  out <- structure(
    list(
      slope       = slope,
      intercept   = intercept,
      citl        = citl,
      ici         = ici,
      curve       = est$curve,
      detectable  = detectable,
      verdict     = verdict,
      abstained   = identical(verdict, "undetermined"),
      n           = n,
      n_boot      = n_boot,
      n_boot_used = n_ok,
      tolerance   = tolerance,
      conf_level  = conf_level,
      call        = cl
    ),
    class = "calibration_suite"
  )
  if (out$abstained) out <- .mark_kernR_abstention(out)
  out
}

#' @export
print.calibration_suite <- function(x, ...) {
  cat("\n  Calibration suite (continuous outcome)\n\n")
  cat(sprintf("  n = %d observations; %d bootstrap replicates (%d usable)\n",
              x$n, x$n_boot, x$n_boot_used))
  fmt <- function(v, ideal, nm) {
    sprintf("  %-12s %8.3f  [%7.3f, %7.3f]   (ideal %g)\n",
            nm, v[["estimate"]], v[["lower"]], v[["upper"]], ideal)
  }
  cat(fmt(x$slope, 1, "slope"))
  cat(fmt(x$intercept, 0, "intercept"))
  cat(fmt(x$citl, 0, "CITL"))
  cat(sprintf("  %-12s %8.3f  [%7.3f, %7.3f]   (ideal 0, outcome units)\n",
              "ICI", x$ici[["estimate"]], x$ici[["lower"]], x$ici[["upper"]]))
  cat(sprintf("\n  Detectable:   %.3f  (smallest slope departure this n resolves)\n",
              x$detectable))
  cat(sprintf("  Tolerance:    %.3f  (what the caller asked to detect)\n",
              x$tolerance))
  cat(sprintf("  Verdict:      %s\n", x$verdict))
  if (isTRUE(x$abstained)) {
    cat("                the interval is wider than the tolerance, so this\n")
    cat("                sample cannot certify calibration -- not a pass\n")
  }
  cat("\n")
  invisible(x)
}

# --- internals ---------------------------------------------------------------

.check_calibration_input <- function(predicted, observed, n_boot, tolerance,
                                     conf_level) {
  if (length(predicted) != length(observed)) {
    stop("`x` and `observed` must have the same length.", call. = FALSE)
  }
  if (!length(predicted)) {
    stop("`x` and `observed` must be non-empty.", call. = FALSE)
  }
  if (any(!is.finite(predicted)) || any(!is.finite(observed))) {
    stop("`x` and `observed` must contain only finite values.", call. = FALSE)
  }
  if (!is.numeric(n_boot) || length(n_boot) != 1L || n_boot < 1) {
    stop("`n_boot` must be a single positive integer.", call. = FALSE)
  }
  if (!is.numeric(tolerance) || length(tolerance) != 1L || tolerance <= 0) {
    stop("`tolerance` must be a single positive number.", call. = FALSE)
  }
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      conf_level <= 0 || conf_level >= 1) {
    stop("`conf_level` must be a single number in (0, 1).", call. = FALSE)
  }
  invisible(TRUE)
}

# The four point estimates plus the smoothed curve. A constant prediction
# vector has no slope, so the linear pieces return NA rather than a fitted
# artefact; the bootstrap drops those resamples.
.calibration_point <- function(predicted, observed) {
  n <- length(predicted)
  degenerate <- stats::var(predicted) == 0 || n < 3L
  if (degenerate) {
    slope <- NA_real_
    intercept <- NA_real_
  } else {
    slope <- unname(stats::coef(stats::lm(observed ~ predicted))[2L])
    intercept <- unname(
      stats::coef(stats::lm(observed ~ 1, offset = predicted))[1L])
  }
  citl <- mean(observed) - mean(predicted)

  # ICI: |smooth(observed | predicted) - predicted|, averaged. loess needs
  # enough distinct support to be meaningful; below that the index is NA and
  # the verdict falls through to "undetermined" on the interval instead.
  ici <- NA_real_
  curve <- NULL
  if (!degenerate && length(unique(predicted)) >= 10L) {
    sm <- try(stats::loess(observed ~ predicted, degree = 1L,
                           surface = "direct"), silent = TRUE)
    if (!inherits(sm, "try-error")) {
      fitted_obs <- stats::fitted(sm)
      ici <- mean(abs(fitted_obs - predicted))
      ord <- order(predicted)
      curve <- data.frame(predicted = predicted[ord],
                          observed_smooth = fitted_obs[ord])
    }
  }
  list(slope = slope, intercept = intercept, citl = citl, ici = ici,
       curve = curve)
}
