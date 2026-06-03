# joint-coverage.R -- joint (multivariate) calibration diagnostic for ensembles.
#
# coverage_test() assesses marginal calibration (per output dimension, pooled);
# it is blind to the dependence structure -- an ensemble with correct margins
# but a wrong correlation passes. joint_coverage_test() closes that gap with a
# multivariate rank histogram: each held-out observation is reduced, together
# with the ensemble, to a single multivariate rank via a pre-rank function, and
# the histogram of those ranks is tested for uniformity. The band-depth
# pre-rank (default) is sensitive to dependence miscalibration that marginal
# methods miss; the average-rank pre-rank reproduces the familiar U-shape =
# under-dispersed reading. This is the joint complement to coverage_test() in
# the validation surface, and the right tool for asking whether an ensemble's
# covariance -- not just its margins -- is calibrated.

#' Joint (Multivariate) Calibration Diagnostic for a Predictive Ensemble
#'
#' Tests whether a predictive ensemble is calibrated *jointly*, including its
#' dependence structure, against held-out observations. Where [coverage_test()]
#' assesses each output dimension separately (marginal calibration), this
#' function builds a multivariate rank histogram: each observation is reduced,
#' together with the ensemble, to one multivariate rank through a pre-rank
#' function, and the histogram of those ranks is tested for uniformity. An
#' ensemble with correct margins but a mis-specified correlation -- which
#' [coverage_test()] passes -- is caught here.
#'
#' The pre-rank function determines what miscalibration is visible:
#' \describe{
#'   \item{`"band_depth"` (default)}{Ranks points by multivariate centrality
#'     (band depth; Thorarinsdottir et al. 2016). Sensitive to dependence /
#'     correlation miscalibration that marginal and average-rank methods miss.
#'     The reading is a *slope*: a mean rank below `0.5` means observations fall
#'     outside the ensemble cloud (jointly under-dispersed); above `0.5` means
#'     they sit too centrally (over-dispersed); a non-uniform histogram with a
#'     central mean signals a dependence error that is not a pure dispersion
#'     shift.}
#'   \item{`"average"`}{Ranks points by the sum of their per-dimension ranks
#'     (Gneiting et al. 2008). Gives the familiar rank-histogram reading -- a
#'     U-shape and a dispersion ratio above one mean under-dispersion, an
#'     inverted-U and a ratio below one mean over-dispersion -- but is weaker
#'     against correlation-only errors.}
#' }
#'
#' Under calibration the multivariate rank is uniform on `{1, ..., n_draws + 1}`
#' for either pre-rank, so the chi-squared uniformity test, the dispersion
#' ratio, and the coverage table are all referenced against the uniform. With
#' few held-out observations the uniformity test has low power; read the
#' dispersion and mean-rank summaries alongside it. A joint test needs at least
#' two output dimensions; for a single output use [coverage_test()].
#'
#' @inheritParams coverage_test
#' @param prerank Character. Pre-rank function: `"band_depth"` (default,
#'   dependence-sensitive) or `"average"` (familiar dispersion reading).
#' @param seed Integer or `NULL`. Random seed: multivariate ranks break
#'   pre-rank ties at random, so a non-`NULL` seed makes the result
#'   reproducible.
#'
#' @returns An object of class `"joint_coverage_test"` with components:
#'   \describe{
#'     \item{prerank}{The pre-rank function used.}
#'     \item{coverage}{Data frame of `nominal` vs `empirical` central-rank-band
#'       coverage at each requested level.}
#'     \item{dispersion_ratio}{`Var(u) / (1/12)` of the normalised multivariate
#'       ranks `u`; above one under-dispersed, below one over-dispersed (the
#'       primary signal for `prerank = "average"`).}
#'     \item{mean_rank}{Mean normalised rank (`0.5` under calibration; the
#'       primary direction signal for `prerank = "band_depth"`).}
#'     \item{calibration}{List with the rank-histogram chi-squared `statistic`,
#'       `p_value`, and `n_bins`.}
#'     \item{verdict}{Character: the joint-calibration classification.}
#'     \item{reject}{Logical: `calibration$p_value <= alpha`.}
#'     \item{ranks}{Integer multivariate ranks, one per observation.}
#'     \item{n_draws, n_obs, dimension, levels, alpha}{Inputs / sizes.}
#'     \item{pesto_metadata}{Provenance carried from a `pesto_ensemble_manifest`
#'       (including `fidelity`), or `NULL`.}
#'   }
#'
#' @references
#' Gneiting, T., Stanberry, L. I., Grimit, E. P., Held, L., & Johnson, N. A.
#' (2008). Assessing probabilistic forecasts of multivariate quantities, with
#' an application to ensemble predictions of surface winds. *TEST*, 17(2),
#' 211-235.
#'
#' Thorarinsdottir, T. L., Scheuerer, M., & Heinz, C. (2016). Assessing the
#' calibration of high-dimensional ensemble forecasts using rank histograms.
#' *Journal of Computational and Graphical Statistics*, 25(1), 105-122.
#'
#' @examples
#' set.seed(1)
#' chol2 <- chol(matrix(c(1, 0.9, 0.9, 1), 2L))
#'
#' # Correctly correlated ensemble: calibrated jointly
#' obs <- matrix(stats::rnorm(120L), ncol = 2L) %*% chol2
#' ens_ok <- matrix(stats::rnorm(4000L), ncol = 2L) %*% chol2
#' joint_coverage_test(ens_ok, obs, seed = 1L)
#'
#' # Right margins, wrong dependence: ensemble independent, obs correlated.
#' # coverage_test() passes; joint_coverage_test() catches it.
#' ens_indep <- matrix(stats::rnorm(4000L), ncol = 2L)
#' joint_coverage_test(ens_indep, obs, seed = 1L)
#'
#' @seealso [coverage_test()], [mmd_ppc()], [concordance_test()]
#' @family goodness-of-fit tests
#' @author Max Moldovan, \email{max.moldovan@@adelaide.edu.au}
#' @export
joint_coverage_test <- function(x, ...) UseMethod("joint_coverage_test")

#' @rdname joint_coverage_test
#' @export
joint_coverage_test.default <- function(x, observed,
                                        prerank = c("band_depth", "average"),
                                        levels = c(0.5, 0.8, 0.9),
                                        n_bins = 10L,
                                        alpha = 0.05,
                                        seed = NULL,
                                        ...) {
  cl <- match.call()
  prerank <- match.arg(prerank)
  ensemble <- as.matrix(x)
  if (missing(observed) || is.null(observed)) {
    stop("`observed` must be supplied when `x` is a matrix.", call. = FALSE)
  }
  observed <- as.matrix(observed)
  .check_coverage_input(ensemble, observed, levels, n_bins, alpha)
  if (ncol(ensemble) < 2L) {
    stop("`joint_coverage_test()` needs at least 2 output dimensions; ",
         "use `coverage_test()` for a single output.", call. = FALSE)
  }

  if (!is.null(seed)) set.seed(seed)

  n_draws <- nrow(ensemble)
  n_obs   <- nrow(observed)
  d       <- ncol(ensemble)
  m_aug   <- n_draws + 1L

  # Multivariate ranks via the chosen pre-rank, normalised to (0, 1) ------
  ranks <- .mv_prerank_ranks(ensemble, observed, prerank)
  u <- (ranks - 0.5) / m_aug

  # Central-rank-band coverage at each nominal level ---------------------
  empirical <- vapply(levels, function(l) {
    lo <- (1 - l) / 2
    hi <- (1 + l) / 2
    mean(u >= lo & u <= hi)
  }, numeric(1L))
  coverage <- data.frame(nominal = levels, empirical = empirical)

  # Dispersion, mean rank, and rank-histogram uniformity -----------------
  dispersion_ratio <- stats::var(u) / (1 / 12)
  mean_rank <- mean(u)
  calib <- .pit_uniformity(u, n_bins)
  verdict <- .joint_coverage_verdict(prerank, calib$p_value, alpha,
                                     dispersion_ratio, mean_rank)

  structure(
    list(
      prerank          = prerank,
      coverage         = coverage,
      dispersion_ratio = dispersion_ratio,
      mean_rank        = mean_rank,
      calibration      = calib,
      verdict          = verdict,
      reject           = calib$p_value <= alpha,
      ranks            = ranks,
      n_draws          = n_draws,
      n_obs            = n_obs,
      dimension        = d,
      levels           = levels,
      alpha            = alpha,
      pesto_metadata   = NULL,
      call             = cl
    ),
    class = "joint_coverage_test"
  )
}

#' @rdname joint_coverage_test
#' @export
joint_coverage_test.pesto_ensemble <- function(x, observed = NULL, ...) {
  if (is.null(observed)) observed <- x$observed
  if (is.null(observed)) {
    stop("`observed` must be supplied: ensemble carries no `observed` slot.",
         call. = FALSE)
  }
  out <- joint_coverage_test.default(x$posterior, observed = observed, ...)
  out$pesto_metadata <- x$metadata
  out$call <- match.call()
  out
}

#' @rdname joint_coverage_test
#' @export
joint_coverage_test.pesto_ensemble_manifest <- function(x, observed,
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

  out <- joint_coverage_test.default(ensemble, observed = observed, ...)
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
print.joint_coverage_test <- function(x, ...) {
  cat("\n  Joint (multivariate) calibration diagnostic\n\n")
  cat("Ensemble:    ", x$n_draws, " draws x ", x$dimension, " dims; ",
      x$n_obs, " held-out obs\n", sep = "")
  cat("Pre-rank:    ", x$prerank,
      if (identical(x$prerank, "band_depth")) {
        " (dependence-sensitive)"
      } else {
        " (dispersion reading)"
      }, "\n", sep = "")
  cov_str <- paste(
    sprintf("%g%%->%.0f%%", 100 * x$coverage$nominal,
            100 * x$coverage$empirical),
    collapse = "   "
  )
  cat("Coverage:    ", cov_str, "  (nominal -> empirical rank band)\n",
      sep = "")
  cat("Mean rank:   ", formatC(x$mean_rank, digits = 3L, format = "f"),
      " (0.5 = calibrated; <0.5 obs outlying, >0.5 obs central)\n", sep = "")
  cat("Dispersion:  ratio = ",
      formatC(x$dispersion_ratio, digits = 3L, format = "f"),
      " (Var(rank)/(1/12))\n", sep = "")
  cat("Calibration: chi2 = ",
      formatC(x$calibration$statistic, digits = 3L, format = "g"),
      ", p = ", formatC(x$calibration$p_value, digits = 4L, format = "f"),
      " (rank histogram, ", x$calibration$n_bins, " bins)\n", sep = "")
  cat("Verdict:     ", x$verdict, "\n", sep = "")
  if (length(x$pesto_metadata)) {
    cat("Metadata:    ", paste(names(x$pesto_metadata), collapse = ", "),
        "\n", sep = "")
  }
  cat("\n")
  invisible(x)
}

# Internal helpers -----------------------------------------------------------

#' Multivariate ranks of observations within an ensemble via a pre-rank
#'
#' For each observation, the observation is augmented into the ensemble
#' (`M = n_draws + 1` points), every point is given its per-dimension rank in
#' the augmented set, those ranks are collapsed to a scalar pre-rank, and the
#' observation's rank among the `M` pre-ranks is returned. Under calibration
#' this rank is uniform on `{1, ..., M}`. The average pre-rank sums the
#' per-dimension ranks (Gneiting et al. 2008); the band-depth pre-rank uses
#' `sum_l (r_l - 1)(M - r_l)`, monotone in the modified band depth
#' (Thorarinsdottir et al. 2016), so a higher value is more central. Pre-rank
#' ties are broken at random, so callers honour `set.seed()`.
#'
#' @param ensemble Numeric `n_draws x d` matrix.
#' @param observed Numeric `n_obs x d` matrix.
#' @param prerank `"band_depth"` or `"average"`.
#' @returns Integer vector of `n_obs` multivariate ranks in `{1, ..., M}`.
#' @noRd
#' @keywords internal
.mv_prerank_ranks <- function(ensemble, observed, prerank) {
  m_aug <- nrow(ensemble) + 1L
  n_obs <- nrow(observed)
  out <- integer(n_obs)
  for (i in seq_len(n_obs)) {
    aug <- rbind(observed[i, , drop = FALSE], ensemble)   # row 1 = observation
    r <- apply(aug, 2L, rank, ties.method = "average")    # M x d per-dim ranks
    pr <- if (identical(prerank, "average")) {
      rowSums(r)
    } else {
      rowSums((r - 1) * (m_aug - r))                       # band-depth pre-rank
    }
    out[i] <- rank(pr, ties.method = "random")[1L]         # obs rank among M
  }
  out
}

#' Classify a joint-calibration result
#'
#' @param prerank The pre-rank used (verdict logic differs by pre-rank).
#' @param p_value Rank-histogram uniformity p-value.
#' @param alpha Verdict level.
#' @param ratio Dispersion ratio `Var(u)/(1/12)`.
#' @param mean_rank Mean normalised multivariate rank.
#' @returns A single human-readable verdict string.
#' @noRd
#' @keywords internal
.joint_coverage_verdict <- function(prerank, p_value, alpha, ratio,
                                    mean_rank) {
  if (p_value > alpha) {
    return("calibrated (multivariate ranks consistent with uniform)")
  }
  if (identical(prerank, "average")) {
    issues <- character(0)
    if (abs(mean_rank - 0.5) > 0.05) {
      issues <- c(issues,
                  if (mean_rank < 0.5) "biased (ensemble too high)"
                  else "biased (ensemble too low)")
    }
    if (ratio > 1.1) {
      issues <- c(issues, "jointly under-dispersed (over-confident)")
    } else if (ratio < 0.9) {
      issues <- c(issues, "jointly over-dispersed")
    }
    if (length(issues) == 0L) {
      issues <- "dependence mis-specified (non-uniform multivariate ranks)"
    }
    return(paste0("REJECT -- ", paste(issues, collapse = "; ")))
  }
  # band_depth: a slope reading (mean rank), with a central-but-non-uniform
  # histogram signalling a dependence error rather than a dispersion shift.
  if (mean_rank < 0.45) {
    return(paste0("REJECT -- jointly under-dispersed ",
                  "(observations outlying vs ensemble)"))
  }
  if (mean_rank > 0.55) {
    return("REJECT -- jointly over-dispersed (observations too central)")
  }
  paste0("REJECT -- dependence mis-specified ",
         "(multivariate ranks non-uniform; not a pure dispersion shift)")
}
