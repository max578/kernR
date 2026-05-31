# DR-DATE for APSIM scenarios -- Year-1 §B2 of the UQ ag-stack roadmap.
# Wraps the existing dr_date_test() (Fawkes-Hu-Evans-Sejdinovic 2024) for
# the two-scenario interpretation: baseline vs intervention APSIM
# ensembles encoded as PESTO 0.3.0 pesto_ensemble_manifest objects.

#' DR-DATE for Two PESTO Ensemble Scenarios
#'
#' Tests whether the *distribution* of simulated outputs differs between
#' two APSIM (or other) scenario ensembles -- e.g. baseline management vs
#' an intervention such as stubble retention -- by running the doubly
#' robust DR-DATE statistic of Fawkes, Hu, Evans & Sejdinovic (2024)
#' over the pooled ensembles.
#'
#' This is a thin scenario-facing wrapper around [dr_date_test()]:
#' parameters are treated as covariates (so the test adjusts for any
#' systematic difference in the parameter posteriors that came from the
#' two PESTO runs), outputs are the outcome, and the scenario label is
#' the binary treatment. Sensitive to distributional differences
#' (variance, shape, tails), not just mean shifts.
#'
#' The PESTO 0.3.0 [`PESTO::pesto_ensemble_manifest`] S7 contract is the
#' supported input shape; the per-realisation file-I/O for ingestion is
#' handled by `PESTO::read_manifest()` upstream of this call.
#'
#' @param baseline A `pesto_ensemble_manifest` (S7) -- the reference
#'   scenario.
#' @param intervention A `pesto_ensemble_manifest` (S7) -- the
#'   alternative scenario. Must share `pesto_version` (major.minor) plus
#'   parameter and observation schemas with `baseline`.
#' @param output Optional character vector of observation column names
#'   to test against. Defaults to all numeric output columns shared by
#'   the two manifests (the `real_name` column is excluded). Pass a
#'   subset to focus the test on specific outputs (e.g. end-of-season
#'   yield only).
#' @param propensity_model Forwarded to [dr_date_test()]. Default
#'   `"logistic"`. In the scenario context the true propensity is 50/50
#'   by design; logistic recovers that and absorbs any sampling
#'   imbalance.
#' @param outcome_model Forwarded to [dr_date_test()]. Default `"krr"`.
#' @param n_permutations Forwarded to [dr_date_test()]. Default 500.
#' @param n_bins Forwarded to [dr_date_test()]. Default 10.
#' @param regularisation Forwarded to [dr_date_test()]. Default `"cv"`.
#' @param alpha Forwarded to [dr_date_test()]. Default 0.05.
#' @param seed Integer or `NULL`. Random seed.
#' @param verbose Logical. Default `FALSE`.
#' @param ... Reserved.
#'
#' @return An object of class
#'   `c("dr_date_scenario", "kernel_test_result")` with the standard
#'   `kernel_test_result` fields plus:
#'   \describe{
#'     \item{baseline_run_id}{Run id from the baseline manifest.}
#'     \item{intervention_run_id}{Run id from the intervention manifest.}
#'     \item{n_baseline}{Realisations in baseline ensemble.}
#'     \item{n_intervention}{Realisations in intervention ensemble.}
#'     \item{outputs_tested}{Character vector of output columns used.}
#'     \item{pesto_versions}{Named character -- baseline / intervention.}
#'   }
#'
#' @references
#' Fawkes, J., Hu, R., Evans, R. J., & Sejdinovic, D. (2024). Doubly
#' robust kernel statistics for testing distributional treatment effects.
#' *Transactions on Machine Learning Research*.
#'
#' @seealso [dr_date_test()] for the underlying observational-causal
#'   test; [`PESTO::pesto_ensemble_manifest`] for the input contract;
#'   [`PESTO::pesto_ies_callback()`] for producing the ensembles upstream.
#'
#' @examples
#' \donttest{
#' # Requires PESTO (>= 0.4.1) -- wired through Imports.
#' library(PESTO)
#' npar <- 2L; nobs <- 4L; nreal <- 60L
#' G  <- matrix(stats::rnorm(nobs * npar), nobs, npar)
#' y0 <- as.numeric(G %*% c(1.0, -0.5)) + stats::rnorm(nobs, sd = 0.05)
#' y1 <- y0 + c(0.6, 0.6, 0.6, 0.6)   # intervention shifts outputs
#' names(y0) <- names(y1) <- paste0("o", seq_len(nobs))
#'
#' prior <- matrix(stats::rnorm(nreal * npar), nreal, npar,
#'                 dimnames = list(NULL, c("p1", "p2")))
#' fit0 <- pesto_ies_callback(function(t) t %*% t(G), prior, y0, 0.05,
#'                            noptmax = 3, verbose = FALSE)
#' fit1 <- pesto_ies_callback(function(t) t %*% t(G), prior, y1, 0.05,
#'                            noptmax = 3, verbose = FALSE)
#' m_base <- as_manifest(fit0, run_id = "baseline")
#' m_intv <- as_manifest(fit1, run_id = "intervention")
#' res <- dr_date_scenario(m_base, m_intv,
#'                          n_permutations = 200L, seed = 1L)
#' print(res)
#' }
#' @family distributional treatment effects
#' @export
dr_date_scenario <- function(baseline, intervention,
                              output            = NULL,
                              propensity_model  = c("logistic", "ranger",
                                                    "xgboost"),
                              outcome_model     = c("krr", "zero"),
                              n_permutations    = 500L,
                              n_bins            = 10L,
                              regularisation    = "cv",
                              alpha             = 0.05,
                              seed              = NULL,
                              verbose           = FALSE,
                              ...) {

  cl <- match.call()
  propensity_model <- match.arg(propensity_model)
  outcome_model    <- match.arg(outcome_model)

  .validate_manifest_pair(baseline, intervention)

  params_b  <- as.data.frame(baseline@params)
  params_i  <- as.data.frame(intervention@params)
  outputs_b <- as.data.frame(baseline@outputs)
  outputs_i <- as.data.frame(intervention@outputs)

  par_cols   <- setdiff(names(params_b),  "real_name")
  out_cols   <- setdiff(names(outputs_b), "real_name")

  if (is.null(output)) {
    output <- out_cols
  } else {
    output <- as.character(output)
    missing_out <- setdiff(output, out_cols)
    if (length(missing_out)) {
      stop("`output` columns not found in manifest outputs: ",
           paste(missing_out, collapse = ", "), call. = FALSE)
    }
  }

  covariates <- rbind(
    as.matrix(params_b[,  par_cols, drop = FALSE]),
    as.matrix(params_i[,  par_cols, drop = FALSE])
  )
  y <- rbind(
    as.matrix(outputs_b[, output, drop = FALSE]),
    as.matrix(outputs_i[, output, drop = FALSE])
  )
  treatment <- c(rep(0L, nrow(params_b)), rep(1L, nrow(params_i)))

  res <- dr_date_test(
    y                = y,
    treatment        = treatment,
    covariates       = covariates,
    propensity_model = propensity_model,
    outcome_model    = outcome_model,
    n_permutations   = n_permutations,
    n_bins           = n_bins,
    regularisation   = regularisation,
    alpha            = alpha,
    seed             = seed,
    verbose          = verbose
  )

  res$method              <- "DR-DATE (scenario)"
  res$baseline_run_id     <- baseline@run_id
  res$intervention_run_id <- intervention@run_id
  res$n_baseline          <- nrow(params_b)
  res$n_intervention      <- nrow(params_i)
  res$outputs_tested      <- output
  res$pesto_versions      <- c(baseline     = baseline@pesto_version,
                               intervention = intervention@pesto_version)
  res$alpha               <- alpha
  res$reject              <- res$p_value <= alpha
  res$call                <- cl
  class(res) <- c("dr_date_scenario", class(res))
  res
}

# Validate that the two manifests can be compared. Hard-stop on
# schema mismatch -- silent comparison of incompatible scenarios would
# be a worse error mode than a noisy failure here.
.validate_manifest_pair <- function(baseline, intervention) {
  if (!.is_pesto_manifest(baseline)) {
    stop("`baseline` must be a PESTO::pesto_ensemble_manifest (S7) ",
         "object. Got: ", paste(class(baseline), collapse = "/"),
         call. = FALSE)
  }
  if (!.is_pesto_manifest(intervention)) {
    stop("`intervention` must be a PESTO::pesto_ensemble_manifest (S7) ",
         "object. Got: ", paste(class(intervention), collapse = "/"),
         call. = FALSE)
  }

  par_b <- setdiff(names(baseline@params),     "real_name")
  par_i <- setdiff(names(intervention@params), "real_name")
  if (!identical(par_b, par_i)) {
    stop("Manifest parameter schemas differ. Baseline: [",
         paste(par_b, collapse = ", "),
         "]; intervention: [", paste(par_i, collapse = ", "), "].",
         call. = FALSE)
  }

  out_b <- setdiff(names(baseline@outputs),     "real_name")
  out_i <- setdiff(names(intervention@outputs), "real_name")
  if (!identical(out_b, out_i)) {
    stop("Manifest observation schemas differ. Baseline outputs: [",
         paste(out_b, collapse = ", "),
         "]; intervention outputs: [",
         paste(out_i, collapse = ", "), "].",
         call. = FALSE)
  }

  vb <- baseline@pesto_version
  vi <- intervention@pesto_version
  if (.major_minor(vb) != .major_minor(vi)) {
    stop("Manifests come from incompatible PESTO versions: baseline=",
         vb, "; intervention=", vi,
         ". Compare only within the same major.minor.", call. = FALSE)
  }

  invisible(NULL)
}

# S7 classes carry a package-qualified S3 class string
# ("PESTO::pesto_ensemble_manifest"), not the bare class name. Accept
# both qualified and bare forms so this works whether or not S7 has
# been re-loaded with a different package-prefix convention.
.is_pesto_manifest <- function(x) {
  inherits(x, c("PESTO::pesto_ensemble_manifest",
                "pesto_ensemble_manifest"))
}

.major_minor <- function(v) {
  parts <- strsplit(as.character(v), "[.\\-]")[[1L]]
  if (length(parts) < 2L) return(as.character(v))
  paste(parts[1L:2L], collapse = ".")
}

#' @export
print.dr_date_scenario <- function(x, ...) {
  NextMethod()
  cat("Scenario contrast\n")
  cat("  baseline      : ", x$baseline_run_id, " (n=", x$n_baseline,
      ")\n", sep = "")
  cat("  intervention  : ", x$intervention_run_id, " (n=",
      x$n_intervention, ")\n", sep = "")
  cat("  outputs tested: ", paste(x$outputs_tested, collapse = ", "),
      "\n", sep = "")
  cat("  PESTO versions: ", paste(names(x$pesto_versions),
                                  unname(x$pesto_versions),
                                  sep = "=", collapse = ", "),
      "\n", sep = "")
  cat("  Verdict:        ",
      if (isTRUE(x$reject)) {
        "REJECT (distributions differ; intervention has effect)"
      } else {
        "fail to reject (no distributional difference detected)"
      },
      "\n", sep = "")
  cat("\n")
  invisible(x)
}
