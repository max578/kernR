#' PESTO Ensemble Manifest (Constructor)
#'
#' Lightweight constructor for a posterior-predictive ensemble produced by
#' PESTO (or any compatible UQ engine). Bundles a posterior-predictive
#' sample matrix with optional held-out observations and free-form
#' metadata, providing a stable interface for [mmd_ppc()] and other
#' downstream verdict-layer tools.
#'
#' This is the kernR-side schema for the cross-package contract; until
#' PESTO ships its native manifest emitter, callers can construct the
#' object directly from in-memory matrices.
#'
#' @param posterior Numeric matrix `M x d`: `M` posterior-predictive
#'   draws over `d` output dimensions (e.g. yield, biomass).
#' @param observed Optional numeric matrix `n_obs x d` of held-out
#'   observations. May be `NULL` when observations are supplied later
#'   to [mmd_ppc()].
#' @param metadata Optional named list of free-form metadata
#'   (run id, ensemble seed, holdout year, etc.).
#'
#' @return An object of class `"pesto_ensemble"`.
#' @seealso [mmd_ppc()]
#'
#' @examples
#' set.seed(1)
#' post <- matrix(stats::rnorm(200L), ncol = 2L)
#' obs  <- matrix(stats::rnorm(20L),  ncol = 2L)
#' ens  <- pesto_ensemble(post, obs, metadata = list(holdout_year = 2018))
#' ens
#'
#' @export
pesto_ensemble <- function(posterior, observed = NULL, metadata = list()) {
  posterior <- as.matrix(posterior)
  if (!is.numeric(posterior) || any(!is.finite(posterior))) {
    stop("`posterior` must be a numeric matrix of finite values.",
         call. = FALSE)
  }
  if (!is.null(observed)) {
    observed <- as.matrix(observed)
    if (ncol(observed) != ncol(posterior)) {
      stop("`observed` must have the same number of columns as `posterior`.",
           call. = FALSE)
    }
    if (any(!is.finite(observed))) {
      stop("`observed` must contain only finite numeric values.",
           call. = FALSE)
    }
  }
  if (!is.list(metadata)) {
    stop("`metadata` must be a list.", call. = FALSE)
  }
  structure(
    list(posterior = posterior, observed = observed, metadata = metadata),
    class = "pesto_ensemble"
  )
}

#' @export
print.pesto_ensemble <- function(x, ...) {
  cat("\n  PESTO ensemble manifest\n\n")
  cat("Posterior: ", nrow(x$posterior), " draws x ",
      ncol(x$posterior), " dims\n", sep = "")
  if (!is.null(x$observed)) {
    cat("Observed:  ", nrow(x$observed), " obs x ",
        ncol(x$observed), " dims\n", sep = "")
  } else {
    cat("Observed:  <none attached>\n")
  }
  if (length(x$metadata)) {
    cat("Metadata:  ", paste(names(x$metadata), collapse = ", "), "\n",
        sep = "")
  }
  cat("\n")
  invisible(x)
}


#' MMD Posterior-Predictive Check
#'
#' Model-free verdict on whether a posterior-predictive ensemble is
#' consistent with held-out observations, via the Maximum Mean Discrepancy
#' (MMD) two-sample test. Wraps [mmd_test()] with a posterior-predictive
#' framing and adds a Shannon-information *surprise* diagnostic
#' (`-log2(p)`) for intuitive interpretation: 0 bits = no surprise
#' (`p = 1`); ~4.32 bits = `p = 0.05`; the maximum achievable surprise
#' at `n_permutations = B` is `log2(B + 1)`.
#'
#' Use after an ensemble-smoother run (PESTO IES, EnKF, etc.) to ask:
#' *does the calibrated model produce predictive draws that match the
#' held-out year / paddock / season at the distributional level?* The
#' MMD test is sensitive to mean, variance, and tail differences --
#' strictly more informative than RMSE on the posterior predictive mean.
#'
#' @param x Either a numeric matrix `M x d` of posterior-predictive
#'   draws, or a `pesto_ensemble` object (see [pesto_ensemble()]).
#' @param observed Numeric matrix `n_obs x d` of held-out observations.
#'   When `x` is a `pesto_ensemble` carrying its own `observed` slot,
#'   may be left `NULL` to use the bundled observations.
#' @param kernel Kernel specification. Default is RBF with median
#'   heuristic over the pooled posterior + observed sample.
#' @param n_permutations Integer. Permutations for the null. Default 500.
#' @param alpha Numeric in `(0, 1)`. Significance level used for the
#'   verdict. Default `0.05`.
#' @param seed Integer or `NULL`. Random seed for reproducibility.
#' @param ... Additional arguments (currently unused; reserved for
#'   future Nystrom acceleration).
#'
#' @return An object of class `c("mmd_ppc", "kernel_test_result")` with
#'   the standard `kernel_test_result` fields plus:
#'   \describe{
#'     \item{n_posterior}{Number of posterior-predictive draws.}
#'     \item{n_observed}{Number of held-out observations.}
#'     \item{surprise_bits}{Shannon-information surprise `-log2(p_value)`.}
#'     \item{alpha}{Verdict significance level.}
#'     \item{reject}{Logical: `p_value <= alpha`.}
#'     \item{pesto_metadata}{Carried through from `pesto_ensemble` input,
#'       when provided; otherwise `NULL`.}
#'   }
#'
#' @references
#' Gretton, A., Borgwardt, K. M., Rasch, M. J., Scholkopf, B., & Smola,
#' A. (2012). A kernel two-sample test. *JMLR*, 13, 723-773.
#'
#' Gelman, A., Meng, X.-L., & Stern, H. (1996). Posterior predictive
#' assessment of model fitness via realized discrepancies. *Statistica
#' Sinica*, 6(4), 733-760.
#'
#' @examples
#' set.seed(1)
#' # Calibrated model: posterior matches truth
#' post <- matrix(stats::rnorm(400L), ncol = 2L)
#' obs  <- matrix(stats::rnorm(40L),  ncol = 2L)
#' fit_ok <- mmd_ppc(post, obs, n_permutations = 199L, seed = 1L)
#' fit_ok
#'
#' # Miscalibrated model: posterior is mean-shifted
#' obs_shift <- obs + 1.5
#' fit_bad <- mmd_ppc(post, obs_shift, n_permutations = 199L, seed = 1L)
#' fit_bad
#'
#' @seealso [mmd_test()], [pesto_ensemble()]
#' @export
mmd_ppc <- function(x, ...) UseMethod("mmd_ppc")

#' @rdname mmd_ppc
#' @export
mmd_ppc.default <- function(x, observed,
                            kernel = kernel_spec(),
                            n_permutations = 500L,
                            alpha = 0.05,
                            seed = NULL,
                            ...) {
  cl <- match.call()
  posterior <- as.matrix(x)
  if (missing(observed) || is.null(observed)) {
    stop("`observed` must be supplied when `x` is a matrix.",
         call. = FALSE)
  }
  observed <- as.matrix(observed)
  if (ncol(posterior) != ncol(observed)) {
    stop("`x` (posterior) and `observed` must have the same number of columns.",
         call. = FALSE)
  }
  if (nrow(posterior) < 5L || nrow(observed) < 5L) {
    stop("Each of posterior and observed must have at least 5 rows.",
         call. = FALSE)
  }
  if (!is.numeric(alpha) || length(alpha) != 1L ||
      alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a single number in (0, 1).", call. = FALSE)
  }

  res <- mmd_test(
    x              = posterior,
    y              = observed,
    kernel         = kernel,
    n_permutations = n_permutations,
    alpha          = alpha,
    seed           = seed
  )

  res$method        <- "MMD PPC"
  res$n_posterior   <- nrow(posterior)
  res$n_observed    <- nrow(observed)
  res$surprise_bits <- -log2(res$p_value)
  res$alpha         <- alpha
  res$reject        <- res$p_value <= alpha
  res$pesto_metadata <- NULL
  res$call          <- cl
  class(res)        <- c("mmd_ppc", class(res))
  res
}

#' @rdname mmd_ppc
#' @export
mmd_ppc.pesto_ensemble <- function(x, observed = NULL, ...) {
  if (is.null(observed)) observed <- x$observed
  if (is.null(observed)) {
    stop("`observed` must be supplied: ensemble carries no `observed` slot.",
         call. = FALSE)
  }
  out <- mmd_ppc.default(x$posterior, observed = observed, ...)
  out$pesto_metadata <- x$metadata
  out$call <- match.call()
  out
}

#' @rdname mmd_ppc
#' @param outputs Optional character vector of output column names from
#'   the manifest to test against. Defaults to all numeric output columns
#'   (the `real_name` column is excluded). Used only for the
#'   `pesto_ensemble_manifest` method.
#' @export
mmd_ppc.pesto_ensemble_manifest <- function(x, observed,
                                            outputs = NULL, ...) {
  cl <- match.call()
  if (missing(observed) || is.null(observed)) {
    stop("`observed` must be supplied for the manifest method. The ",
         "manifest's `obs_target` slot is a single nobs-dim point ",
         "(not a sample) and is unsuitable for two-sample PPC. ",
         "Pass a held-out matrix with at least 5 rows.",
         call. = FALSE)
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

  posterior <- as.matrix(outputs_df[, outputs, drop = FALSE])
  observed  <- as.matrix(observed)
  if (ncol(observed) != ncol(posterior)) {
    stop("`observed` must have the same number of columns as the ",
         "selected manifest outputs (got ", ncol(observed),
         "; expected ", ncol(posterior), ").", call. = FALSE)
  }

  out <- mmd_ppc.default(posterior, observed = observed, ...)
  out$pesto_metadata <- list(
    run_id        = x@run_id,
    pesto_version = x@pesto_version,
    method        = x@method,
    outputs_used  = outputs
  )
  out$call <- cl
  out
}

#' @export
print.mmd_ppc <- function(x, ...) {
  NextMethod()
  cat("PPC verdict\n")
  cat("  Posterior:  ", x$n_posterior, " draws\n", sep = "")
  cat("  Observed:   ", x$n_observed, " obs\n", sep = "")
  cat("  Surprise:   ",
      formatC(x$surprise_bits, digits = 3L, format = "f"), " bits\n",
      sep = "")
  cat("  Verdict:    ",
      if (isTRUE(x$reject)) "REJECT (posterior inconsistent with observations)"
      else "consistent with observations",
      "\n", sep = "")
  if (length(x$pesto_metadata)) {
    cat("  Metadata:   ", paste(names(x$pesto_metadata), collapse = ", "),
        "\n", sep = "")
  }
  cat("\n")
  invisible(x)
}
