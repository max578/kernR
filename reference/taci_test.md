# Theory-anchored causal inference (TACI) mechanism-consistency test

Tests whether an observed treatment effect is consistent with the effect
a calibrated mechanistic model predicts. Unlike a model-free causal
test, which asks only whether treatment and outcome are associated, TACI
anchors both the null and the alternative in a process model: the
reference distribution is built from the model's own
posterior-predictive draws. The verdict is three-way – the data are
consistent with the model-implied effect, the data show an effect the
model does not predict, or there is no detectable effect.

## Usage

``` r
taci_test(
  posterior,
  mechanism,
  X,
  treatment,
  outcome,
  confounders = NULL,
  density_ratio = "logistic",
  h0_mode = c("permute_within_model", "model_without_treatment"),
  treatment_type = c("auto", "binary", "continuous"),
  baseline = NULL,
  noise_sd = NULL,
  n_perm = 300L,
  alpha = 0.05,
  mechanism_provenance = NULL,
  posterior_provenance = NULL,
  seed = NULL
)
```

## Arguments

- posterior:

  A numeric matrix or data.frame of posterior parameter draws, one row
  per draw, columns in the order the `mechanism` expects. Typically the
  parameter ensemble of a fitted simulator (for example a PESTO IES
  posterior); any numeric draw matrix is accepted.

- mechanism:

  A function `mechanism(theta, X, t)` that, given one posterior draw
  `theta` (a numeric vector), a covariate matrix `X` (`n` by `p`), and a
  treatment vector `t` (length `n` or scalar), returns the model-implied
  mean outcome `E[Y]` as a length-`n` numeric vector.

- X:

  Numeric matrix of covariates, `n` rows. Use a single constant column
  when the mechanism has no covariate dependence.

- treatment:

  Numeric treatment vector of length `n`. Binary or continuous; the
  construction adapts via `treatment_type`.

- outcome:

  Numeric outcome vector of length `n`.

- confounders:

  Numeric matrix of backdoor confounders, or `NULL` (default) for the
  unadjusted statistic. When supplied, density-ratio weights break the
  backdoor path so an observational treatment is handled correctly.

- density_ratio:

  Character density-ratio backend used for backdoor adjustment, passed
  to
  [`fit_density_ratio()`](https://max578.github.io/kernR/reference/fit_density_ratio.md).
  Default is `"logistic"`.

- h0_mode:

  Character null construction. `"permute_within_model"` (default)
  permutes the model-implied outcome to break the treatment association
  while keeping the model's marginal; `"model_without_treatment"`
  simulates the outcome with treatment held at the control baseline so
  the model itself asserts no effect.

- treatment_type:

  Character. `"auto"` (default) calls a treatment with two or fewer
  distinct levels binary, otherwise continuous. `"binary"` recovers the
  `t = 0` to `t = 1` contrast exactly; `"continuous"` contrasts across
  the observed dose range.

- baseline:

  Numeric control level at which the mechanism is switched off for the
  null. Defaults to `0` for a binary treatment and the treatment mean
  for a continuous one.

- noise_sd:

  Numeric observation-noise standard deviation for the model-implied
  draws, or `NULL` (default) to estimate it from the model residual at
  the posterior mean. The model residual is shape-correct for a
  saturating dose-response, where a treatment-only detrend would leave
  curvature and over-noise the reference band.

- n_perm:

  Integer number of model-implied reference replicates. Default is
  `300L`.

- alpha:

  Numeric significance level. Default is `0.05`.

- mechanism_provenance:

  Optional. A record of where the `mechanism`'s calibration came from
  (e.g. a PESTO manifest `run_id` + `apsim_version`, a citation, a
  fitted-model handle). TACI builds its entire reference band from the
  `mechanism` and cannot itself verify the calibration corresponds to
  reality; supplying this declares the mechanism grounded. When `NULL`
  (the default), the result's `grounding` is `"[unverified]"` and the
  human-facing `verdict` string is suffixed `[unverified]` so a verdict
  built on an un-grounded mechanism is never presented as
  unconditionally confident (Independent Oracle Principle). The
  `decision` enum is unchanged.

- posterior_provenance:

  Optional. The analogous record for the `posterior` draws; carried
  through to the result for completeness.

- seed:

  Integer random seed, or `NULL`. Set it for a reproducible reference
  distribution.

## Value

An object of class `"taci_result"`: a list carrying the observed bd-HSIC
statistic, the H0 tail p-value and in-tail flag, the H1 central
interval, the H1 consistency flag and percentile, a `borderline` flag,
the three-way `decision`, a `posterior_adequacy` diagnostic, and the
reference draws. A [`print()`](https://rdrr.io/r/base/print.html) method
is provided.

## Details

The statistic is a weighted bd-HSIC between treatment and outcome,
identical to the engine
[`bd_hsic_test()`](https://max578.github.io/kernR/reference/bd_hsic_test.md)
uses. TACI differs in the reference: rather than a permutation null, it
draws `n_perm` model-implied replicates, sampling a posterior draw with
replacement each time so the band integrates over the posterior. The
alternative band H1 is the model's prediction at the observed treatment;
the null band H0 is one of the two `h0_mode` constructions. The observed
statistic is read against both: a small H0 tail p-value means an effect
is present, and consistency with the H1 band means that effect matches
the model's prediction.

A posterior-adequacy guard protects against a degenerate H1 band. When
the posterior pins the model-implied *effect* too precisely (effect
coefficient of variation below `0.02`), "consistency with H1" is not
meaningful; the guard warns and flags `posterior_adequacy$ok = FALSE`.
The remedy is to widen the posterior at its source, for example with
ensemble inflation in the simulator calibration. The guard reads the
spread of the effect itself, not per- parameter spread, so a
well-identified nuisance covariate does not trip it.

## References

Hu, R., Sejdinovic, D., & Evans, R. J. (2024). A kernel test for causal
association via noise contrastive backdoor adjustment. *JMLR*, 25(160),
1-56.

## See also

Other causal association tests:
[`bd_hsic_test()`](https://max578.github.io/kernR/reference/bd_hsic_test.md),
[`hierarchical_test()`](https://max578.github.io/kernR/reference/hierarchical_test.md),
[`kernel_causal_test()`](https://max578.github.io/kernR/reference/kernel_causal_test.md)

## Examples

``` r
set.seed(1)
n <- 80
nrate <- runif(n, 0, 200)
yield <- 1.1 + 4.2 * (1 - exp(-0.018 * nrate)) + rnorm(n, 0, 0.25)
post <- cbind(ymax = rnorm(200, 4.2, 0.30),
              rate = rnorm(200, 0.018, 0.004),
              y0   = rnorm(200, 1.1, 0.15))
mitscherlich <- function(theta, X, t) {
  theta[3] + theta[1] * (1 - exp(-theta[2] * t))
}
res <- taci_test(post, mitscherlich, X = matrix(1, n, 1),
                 treatment = nrate, outcome = yield, n_perm = 100, seed = 1)
print(res)
#> TACI mechanism-consistency test
#>   treatment: continuous (H0 baseline = 105.1)
#>   statistic: unadjusted bd-HSIC
#>   observed bd-HSIC: 0.06608
#>   H0 tail p-value:  0.010  (in tail: TRUE)
#>   H1 central [0.05187, 0.07807]  obs at H1 pctile 0.35  consistent: TRUE
#>   DECISION: MECHANISM_CONSISTENT_EFFECT
#>   GROUNDING: [unverified] (mechanism provenance not declared -- verdict not grounded)
res$grounding            # "[unverified]" -- no mechanism provenance declared
#> [1] "[unverified]"

# Declaring mechanism provenance grounds the verdict (Independent Oracle
# Principle). TACI builds its whole reference band from `mechanism`, which it
# cannot itself check against reality; naming where the calibration came from
# -- a citation, a fitted-model handle, or a PESTO manifest's run id and the
# simulator version it was calibrated with -- moves `grounding` from
# "[unverified]" to "grounded".
res_grounded <- taci_test(
  post, mitscherlich, X = matrix(1, n, 1),
  treatment = nrate, outcome = yield, n_perm = 100, seed = 1,
  mechanism_provenance = list(
    run_id = "ies-2026-06-12-0042",
    simulator = "APSIM NG 2024.6.7579",
    reference = "Mitscherlich (1909) N-response form"
  )
)
res_grounded$grounding   # "grounded"
#> [1] "grounded"
```
