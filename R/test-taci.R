# test-taci.R -- Theory-Anchored Causal Inference (TACI) mechanism-consistency test.
#
# TACI asks a question the model-free tests cannot: does an observed treatment
# effect agree with the effect a calibrated mechanistic model predicts? The
# reference distribution is built from the model's own posterior-predictive
# draws, so the null and the alternative are both mechanism-implied. The test
# statistic is kernR's weighted bd-HSIC, assembled here against that custom
# reference via the exported `weighted_hsic_stat()` / `resolve_bandwidth()`
# primitives -- TACI's statistic IS kernR's bd-HSIC, scored against a model band
# rather than a permutation band. Backdoor adjustment reuses kernR's
# density-ratio machinery so an observational, confounded treatment is handled
# correctly.

# --- internal helpers --------------------------------------------------------

#' Gaussian gram matrix at the median-heuristic bandwidth
#'
#' @param v Numeric vector or matrix of observations.
#' @returns A square kernel matrix, the RBF gram at the data-resolved bandwidth.
#' @noRd
#' @keywords internal
.taci_kgram <- function(v) {
  v <- as.matrix(v)
  kernel_matrix(v, kernel = resolve_bandwidth(kernel_spec(), v))
}

#' Backdoor density-ratio weights w(t, z) = p*(t) / p(t | z)
#'
#' Breaks the T <- Z -> Y confounding path by reweighting the observational
#' sample toward the interventional treatment marginal. Returns unit weights
#' (plain HSIC) when there are no confounders. Carries the proxymix C4
#' fit-quality flag through as an attribute when that backend is used.
#'
#' @param treatment Numeric treatment vector.
#' @param z Confounder matrix, or `NULL` for the unadjusted statistic.
#' @param method Character density-ratio backend passed to [fit_density_ratio()].
#' @returns A numeric weight vector with a `fit_quality` attribute.
#' @noRd
#' @keywords internal
.taci_backdoor_weights <- function(treatment, z, method) {
  if (is.null(z)) {
    return(rep(1, length(treatment)))
  }
  fit <- fit_density_ratio(x = as.matrix(treatment), z = as.matrix(z),
                           method = method)
  w <- as.numeric(predict_density_ratio(fit, as.matrix(treatment),
                                         as.matrix(z), type = "weight"))
  attr(w, "fit_quality") <- fit$fit_quality
  w
}

#' Coerce a posterior argument to a numeric draw matrix
#'
#' @param posterior A numeric matrix or data.frame of posterior parameter draws,
#'   one row per draw.
#' @returns A numeric matrix of draws.
#' @noRd
#' @keywords internal
.taci_posterior <- function(posterior) {
  if (!is.matrix(posterior) && !is.data.frame(posterior)) {
    stop("`posterior` must be a matrix or data.frame of parameter draws ",
         "(one row per draw).", call. = FALSE)
  }
  m <- as.matrix(posterior)
  if (!is.numeric(m)) {
    stop("`posterior` must contain only numeric parameter columns; drop any ",
         "label column before calling `taci_test()`.", call. = FALSE)
  }
  m
}

# --- the TACI test -----------------------------------------------------------

#' Theory-anchored causal inference (TACI) mechanism-consistency test
#'
#' Tests whether an observed treatment effect is consistent with the effect a
#' calibrated mechanistic model predicts. Unlike a model-free causal test, which
#' asks only whether treatment and outcome are associated, TACI anchors both the
#' null and the alternative in a process model: the reference distribution is
#' built from the model's own posterior-predictive draws. The verdict is
#' three-way -- the data are consistent with the model-implied effect, the data
#' show an effect the model does not predict, or there is no detectable effect.
#'
#' @param posterior A numeric matrix or data.frame of posterior parameter draws,
#'   one row per draw, columns in the order the `mechanism` expects. Typically
#'   the parameter ensemble of a fitted simulator (for example a PESTO IES
#'   posterior); any numeric draw matrix is accepted.
#' @param mechanism A function `mechanism(theta, X, t)` that, given one posterior
#'   draw `theta` (a numeric vector), a covariate matrix `X` (`n` by `p`), and a
#'   treatment vector `t` (length `n` or scalar), returns the model-implied mean
#'   outcome `E[Y]` as a length-`n` numeric vector.
#' @param X Numeric matrix of covariates, `n` rows. Use a single constant column
#'   when the mechanism has no covariate dependence.
#' @param treatment Numeric treatment vector of length `n`. Binary or continuous;
#'   the construction adapts via `treatment_type`.
#' @param outcome Numeric outcome vector of length `n`.
#' @param confounders Numeric matrix of backdoor confounders, or `NULL` (default)
#'   for the unadjusted statistic. When supplied, density-ratio weights break the
#'   backdoor path so an observational treatment is handled correctly.
#' @param density_ratio Character density-ratio backend used for backdoor
#'   adjustment, passed to [fit_density_ratio()]. Default is `"logistic"`.
#' @param h0_mode Character null construction. `"permute_within_model"` (default)
#'   permutes the model-implied outcome to break the treatment association while
#'   keeping the model's marginal; `"model_without_treatment"` simulates the
#'   outcome with treatment held at the control baseline so the model itself
#'   asserts no effect.
#' @param treatment_type Character. `"auto"` (default) calls a treatment with two
#'   or fewer distinct levels binary, otherwise continuous. `"binary"` recovers
#'   the `t = 0` to `t = 1` contrast exactly; `"continuous"` contrasts across the
#'   observed dose range.
#' @param mechanism_provenance Optional. A record of where the `mechanism`'s
#'   calibration came from (e.g. a PESTO manifest `run_id` + `apsim_version`, a
#'   citation, a fitted-model handle). TACI builds its entire reference band
#'   from the `mechanism` and cannot itself verify the calibration corresponds
#'   to reality; supplying this declares the mechanism grounded. When `NULL`
#'   (the default), the result's `grounding` is `"[unverified]"` and the
#'   human-facing `verdict` string is suffixed `[unverified]` so a verdict built
#'   on an un-grounded mechanism is never presented as unconditionally
#'   confident (Independent Oracle Principle). The `decision` enum is unchanged.
#' @param posterior_provenance Optional. The analogous record for the
#'   `posterior` draws; carried through to the result for completeness.
#' @param baseline Numeric control level at which the mechanism is switched off
#'   for the null. Defaults to `0` for a binary treatment and the treatment mean
#'   for a continuous one.
#' @param noise_sd Numeric observation-noise standard deviation for the
#'   model-implied draws, or `NULL` (default) to estimate it from the model
#'   residual at the posterior mean. The model residual is shape-correct for a
#'   saturating dose-response, where a treatment-only detrend would leave
#'   curvature and over-noise the reference band.
#' @param n_perm Integer number of model-implied reference replicates. Default is
#'   `300L`.
#' @param alpha Numeric significance level. Default is `0.05`.
#' @param seed Integer random seed, or `NULL`. Set it for a reproducible
#'   reference distribution.
#'
#' @returns An object of class `"taci_result"`: a list carrying the observed
#'   bd-HSIC statistic, the H0 tail p-value and in-tail flag, the H1 central
#'   interval, the H1 consistency flag and percentile, a `borderline` flag, the
#'   three-way `decision`, a `posterior_adequacy` diagnostic, and the reference
#'   draws. A `print()` method is provided.
#'
#' @details
#' The statistic is a weighted bd-HSIC between treatment and outcome, identical
#' to the engine [bd_hsic_test()] uses. TACI differs in the reference: rather
#' than a permutation null, it draws `n_perm` model-implied replicates, sampling
#' a posterior draw with replacement each time so the band integrates over the
#' posterior. The alternative band H1 is the model's prediction at the observed
#' treatment; the null band H0 is one of the two `h0_mode` constructions. The
#' observed statistic is read against both: a small H0 tail p-value means an
#' effect is present, and consistency with the H1 band means that effect matches
#' the model's prediction.
#'
#' A posterior-adequacy guard protects against a degenerate H1 band. When the
#' posterior pins the model-implied *effect* too precisely (effect coefficient
#' of variation below `0.02`), "consistency with H1" is not meaningful; the guard
#' warns and flags `posterior_adequacy$ok = FALSE`. The remedy is to widen the
#' posterior at its source, for example with ensemble inflation in the simulator
#' calibration. The guard reads the spread of the effect itself, not per-
#' parameter spread, so a well-identified nuisance covariate does not trip it.
#'
#' @references
#' Hu, R., Sejdinovic, D., & Evans, R. J. (2024). A kernel test for causal
#' association via noise contrastive backdoor adjustment. *JMLR*, 25(160), 1-56.
#'
#' @examples
#' set.seed(1)
#' n <- 80
#' nrate <- runif(n, 0, 200)
#' yield <- 1.1 + 4.2 * (1 - exp(-0.018 * nrate)) + rnorm(n, 0, 0.25)
#' post <- cbind(ymax = rnorm(200, 4.2, 0.30),
#'               rate = rnorm(200, 0.018, 0.004),
#'               y0   = rnorm(200, 1.1, 0.15))
#' mitscherlich <- function(theta, X, t) {
#'   theta[3] + theta[1] * (1 - exp(-theta[2] * t))
#' }
#' res <- taci_test(post, mitscherlich, X = matrix(1, n, 1),
#'                  treatment = nrate, outcome = yield, n_perm = 100, seed = 1)
#' print(res)
#' res$grounding            # "[unverified]" -- no mechanism provenance declared
#'
#' # Declaring mechanism provenance grounds the verdict (Independent Oracle
#' # Principle). TACI builds its whole reference band from `mechanism`, which it
#' # cannot itself check against reality; naming where the calibration came from
#' # -- a citation, a fitted-model handle, or a PESTO manifest's run id and the
#' # simulator version it was calibrated with -- moves `grounding` from
#' # "[unverified]" to "grounded".
#' res_grounded <- taci_test(
#'   post, mitscherlich, X = matrix(1, n, 1),
#'   treatment = nrate, outcome = yield, n_perm = 100, seed = 1,
#'   mechanism_provenance = list(
#'     run_id = "ies-2026-06-12-0042",
#'     simulator = "APSIM NG 2024.6.7579",
#'     reference = "Mitscherlich (1909) N-response form"
#'   )
#' )
#' res_grounded$grounding   # "grounded"
#'
#' @family causal association tests
#' @export
taci_test <- function(posterior, mechanism, X, treatment, outcome,
                      confounders = NULL,
                      density_ratio = "logistic",
                      h0_mode = c("permute_within_model",
                                  "model_without_treatment"),
                      treatment_type = c("auto", "binary", "continuous"),
                      baseline = NULL,
                      noise_sd = NULL,
                      n_perm = 300L,
                      alpha = 0.05,
                      mechanism_provenance = NULL,
                      posterior_provenance = NULL,
                      seed = NULL) {
  h0_mode <- match.arg(h0_mode)
  treatment_type <- match.arg(treatment_type)
  if (!is.function(mechanism)) {
    stop("`mechanism` must be a function(theta, X, t) returning E[Y].",
         call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)

  Theta <- .taci_posterior(posterior)
  X  <- as.matrix(X)
  Tt <- as.numeric(treatment)
  Y  <- as.numeric(outcome)
  n  <- length(Y)
  S  <- nrow(Theta)

  # Resolve the treatment contrast ----------------------------------------
  # "auto" calls anything with two or fewer distinct levels binary. The
  # continuous path generalises the dose endpoints of the model-implied effect
  # contrast and the baseline at which the mechanism is switched off for the
  # null; binary recovers the t = 0 -> t = 1 construction exactly.
  if (treatment_type == "auto") {
    treatment_type <- if (length(unique(Tt)) <= 2L) "binary" else "continuous"
  }
  if (treatment_type == "binary") {
    t_lo <- 0
    t_hi <- 1
    if (is.null(baseline)) baseline <- 0L
  } else {
    t_lo <- min(Tt)
    t_hi <- max(Tt)
    if (is.null(baseline)) baseline <- mean(Tt)
  }

  # Noise scale for the model-implied draws -------------------------------
  # TACI is conditional on the calibrated mechanism, so the natural noise scale
  # is the model residual at the posterior mean -- Y minus the model's own
  # prediction at the observed treatment. A treatment-only detrend of a
  # saturating dose-response leaves curvature in the residual and over-noises
  # the reference draws; the model residual is shape-correct. Falls back to a
  # treatment detrend only when the mechanism prediction is non-finite.
  if (is.null(noise_sd)) {
    mhat <- tryCatch(mechanism(colMeans(Theta), X, Tt),
                     error = function(e) rep(NA_real_, n))
    noise_sd <- if (all(is.finite(mhat))) {
      stats::sd(Y - mhat)
    } else {
      detrend <- if (treatment_type == "binary") {
        stats::lm(Y ~ factor(Tt))
      } else {
        stats::lm(Y ~ stats::poly(Tt, min(2L, length(unique(Tt)) - 1L)))
      }
      stats::sd(stats::residuals(detrend))
    }
  }

  # Posterior-adequacy guard ----------------------------------------------
  # The concern is a model-implied treatment EFFECT the posterior pins too
  # precisely: then the H1 band is degenerate and "consistency with H1" is
  # meaningless. The right quantity is the posterior spread of the effect itself
  # -- mean[ M(theta; X, t_hi) - M(theta; X, t_lo) ] across draws -- not per-
  # parameter spread, so a well-identified nuisance covariate does not trip it.
  # Computed over ALL draws (deterministic), so the guard consumes no RNG and
  # never perturbs the reference-loop stream.
  .eff_of <- function(th) {
    mean(mechanism(th, X, rep(t_hi, n)) - mechanism(th, X, rep(t_lo, n)))
  }
  eff    <- apply(Theta, 1L, .eff_of)
  eff_cv <- stats::sd(eff) / (abs(mean(eff)) + 1e-12)
  post_ok <- eff_cv >= 0.02
  adequacy <- list(
    ok = post_ok, effect_cv = eff_cv, effect_mean = mean(eff),
    reason = if (post_ok) NA_character_ else paste0(
      "model-implied treatment effect is over-determined (effect coeff. of ",
      "variation ", sprintf("%.3g", eff_cv), " < 0.02): the posterior pins ",
      "the effect too precisely, so the H1 band is degenerate and the ",
      "mechanism-consistency verdict is unreliable -- widen the posterior ",
      "(for example with ensemble inflation) at its source first"))
  if (!post_ok) {
    warning("taci_test(): ", adequacy$reason, call. = FALSE)
  }

  # Backdoor-adjusted weighted bd-HSIC ------------------------------------
  # Density-ratio weights break the T <- Z -> Y backdoor path so an
  # observational, confounded treatment is handled correctly. The treatment
  # kernel and weights are fixed across replicates (only Y varies), so the H0
  # null permutes the outcome with the weighting held -- matching
  # bd_hsic_test()'s permute-Y null.
  adjusted <- !is.null(confounders)
  w   <- .taci_backdoor_weights(Tt, confounders, density_ratio)
  Kt  <- .taci_kgram(Tt)
  whsic <- function(yv) weighted_hsic_stat(Kt, .taci_kgram(yv), w)
  obs_stat <- whsic(Y)

  # Model-implied reference draws -----------------------------------------
  # Each replicate samples a posterior draw with replacement, integrating over
  # p(theta | D_obs). H1 is the model's prediction at the OBSERVED treatment. H0
  # is one of the two constructions: permute the model outcome to break the
  # association, or simulate with treatment held at control so the model itself
  # asserts no effect.
  h1 <- numeric(n_perm)
  h0 <- numeric(n_perm)
  t_control <- rep(baseline, n)
  for (b in seq_len(n_perm)) {
    theta <- Theta[sample.int(S, 1L), ]
    Ysim  <- mechanism(theta, X, Tt) + stats::rnorm(n, sd = noise_sd)
    h1[b] <- whsic(Ysim)
    h0[b] <- if (h0_mode == "model_without_treatment") {
      whsic(mechanism(theta, X, t_control) + stats::rnorm(n, sd = noise_sd))
    } else {
      whsic(sample(Ysim))
    }
  }

  # Verdict ----------------------------------------------------------------
  p_h0 <- (1 + sum(h0 >= obs_stat)) / (1 + length(h0))
  in_tail <- p_h0 < alpha
  q <- stats::quantile(h1, c(alpha / 2, 1 - alpha / 2), names = FALSE)
  h1_consistent <- obs_stat >= q[1] && obs_stat <= q[2]
  # Where obs sits within H1 (0 = below the model's effect, 1 = above). Values
  # near the interval edges are BORDERLINE -- the label is then fragile.
  h1_percentile <- mean(h1 <= obs_stat)
  borderline <- abs(h1_percentile - (1 - alpha / 2)) < 0.05 ||
                abs(h1_percentile - (alpha / 2)) < 0.05

  decision <- if (!in_tail) {
    "no_effect"
  } else if (h1_consistent) {
    "mechanism_consistent_effect"
  } else {
    "mechanism_inconsistent_effect"
  }

  # FX-9 (Independent Oracle Principle): TACI's null + alternative bands are
  # built entirely from the caller-supplied `mechanism` (+ `posterior`), and
  # nothing here verifies that the mechanism's calibration corresponds to
  # reality. If the mechanism encodes a hallucinated APSIM fact, the verdict is
  # confidently wrong. We therefore LABEL the verdict's grounding rather than
  # emit it as unconditionally confident: `grounding` is "grounded" only when
  # the caller declares where the mechanism came from (`mechanism_provenance`),
  # else "[unverified]" (orchestra provenance vocabulary). The `decision` enum
  # is left intact for downstream consumers; the honesty rides in `grounding`
  # and the human-facing `verdict` string. Pure metadata -- nil statistical risk.
  grounding <- if (is.null(mechanism_provenance)) "[unverified]" else "grounded"
  verdict <- if (identical(grounding, "[unverified]") &&
                 decision %in% c("mechanism_consistent_effect",
                                 "mechanism_inconsistent_effect")) {
    paste0(decision, " [unverified]")
  } else {
    decision
  }

  structure(
    list(observed_statistic = obs_stat,
         p_h0 = p_h0, in_tail = in_tail,
         h1_interval = q, h1_consistent = h1_consistent,
         h1_percentile = h1_percentile, borderline = borderline,
         decision = decision,
         grounding = grounding, verdict = verdict,
         mechanism_provenance = if (is.null(mechanism_provenance)) NA
                                else mechanism_provenance,
         posterior_provenance = if (is.null(posterior_provenance)) NA
                                else posterior_provenance,
         posterior_adequacy = adequacy,
         treatment_type = treatment_type, baseline = baseline,
         adjusted = adjusted, h0_mode = h0_mode,
         density_ratio = if (adjusted) density_ratio else NA_character_,
         dr_fit_quality = attr(w, "fit_quality"),
         h0_ref = h0, h1_ref = h1,
         n = n, n_posterior = S, alpha = alpha, noise_sd = noise_sd),
    class = "taci_result"
  )
}
