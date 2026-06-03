# DR-DATE for Two PESTO Ensemble Scenarios

Tests whether the *distribution* of simulated outputs differs between
two APSIM (or other) scenario ensembles – e.g. baseline management vs an
intervention such as stubble retention – by running the doubly robust
DR-DATE statistic of Fawkes, Hu, Evans & Sejdinovic (2024) over the
pooled ensembles.

## Usage

``` r
dr_date_scenario(
  baseline,
  intervention,
  output = NULL,
  propensity_model = c("logistic", "ranger", "xgboost"),
  outcome_model = c("krr", "zero"),
  n_permutations = 500L,
  n_bins = 10L,
  regularisation = "cv",
  alpha = 0.05,
  seed = NULL,
  verbose = FALSE,
  strict_fidelity = FALSE,
  ...
)
```

## Arguments

- baseline:

  A `pesto_ensemble_manifest` (S7) – the reference scenario.

- intervention:

  A `pesto_ensemble_manifest` (S7) – the alternative scenario. Must
  share `pesto_version` (major.minor) plus parameter and observation
  schemas with `baseline`.

- output:

  Optional character vector of observation column names to test against.
  Defaults to all numeric output columns shared by the two manifests
  (the `real_name` column is excluded). Pass a subset to focus the test
  on specific outputs (e.g. end-of-season yield only).

- propensity_model:

  Forwarded to
  [`dr_date_test()`](https://max578.github.io/kernR/reference/dr_date_test.md).
  Default `"logistic"`. In the scenario context the true propensity is
  50/50 by design; logistic recovers that and absorbs any sampling
  imbalance.

- outcome_model:

  Forwarded to
  [`dr_date_test()`](https://max578.github.io/kernR/reference/dr_date_test.md).
  Default `"krr"`.

- n_permutations:

  Forwarded to
  [`dr_date_test()`](https://max578.github.io/kernR/reference/dr_date_test.md).
  Default 500.

- n_bins:

  Forwarded to
  [`dr_date_test()`](https://max578.github.io/kernR/reference/dr_date_test.md).
  Default 10.

- regularisation:

  Forwarded to
  [`dr_date_test()`](https://max578.github.io/kernR/reference/dr_date_test.md).
  Default `"cv"`.

- alpha:

  Forwarded to
  [`dr_date_test()`](https://max578.github.io/kernR/reference/dr_date_test.md).
  Default 0.05.

- seed:

  Integer or `NULL`. Random seed.

- verbose:

  Logical. Default `FALSE`.

- strict_fidelity:

  Logical. If `FALSE` (default) a mismatch in the two manifests'
  multi-fidelity provenance raises a `warning`; if `TRUE` it raises an
  error. See the *Fidelity provenance* section.

- ...:

  Reserved.

## Value

An object of class `c("dr_date_scenario", "kernel_test_result")` with
the standard `kernel_test_result` fields plus:

- baseline_run_id:

  Run id from the baseline manifest.

- intervention_run_id:

  Run id from the intervention manifest.

- n_baseline:

  Realisations in baseline ensemble.

- n_intervention:

  Realisations in intervention ensemble.

- outputs_tested:

  Character vector of output columns used.

- pesto_versions:

  Named character – baseline / intervention.

- fidelity:

  List with `baseline` / `intervention` fidelity provenance from the
  PESTO manifests (a
  `list(type, schedule, final_level, n_levels, costs)` for a
  multi-fidelity run, or `NULL` for a single-fidelity run).

## Details

This is a thin scenario-facing wrapper around
[`dr_date_test()`](https://max578.github.io/kernR/reference/dr_date_test.md):
parameters are treated as covariates (so the test adjusts for any
systematic difference in the parameter posteriors that came from the two
PESTO runs), outputs are the outcome, and the scenario label is the
binary treatment. Sensitive to distributional differences (variance,
shape, tails), not just mean shifts.

The PESTO 0.3.0
[`PESTO::pesto_ensemble_manifest`](https://rdrr.io/pkg/PESTO/man/pesto_ensemble_manifest.html)
S7 contract is the supported input shape; the per-realisation file-I/O
for ingestion is handled by
[`PESTO::read_manifest()`](https://rdrr.io/pkg/PESTO/man/read_manifest.html)
upstream of this call.

## Fidelity provenance

A PESTO multi-fidelity run records, in the manifest `fidelity` slot,
which fidelity levels produced the ensemble. Comparing a baseline and an
intervention that were calibrated at different fidelities risks
confounding the distributional contrast with fidelity bias. This wrapper
therefore surfaces a provenance mismatch – one scenario single-fidelity
and the other multi-fidelity, or two multi-fidelity runs with different
stack shapes / final levels – as a `warning` (default) or, with
`strict_fidelity = TRUE`, a hard error. Matched or both-single-fidelity
provenance passes silently. The check is forward-compatible: manifests
from PESTO versions that do not populate the slot read as `NULL` and
pass.

## References

Fawkes, J., Hu, R., Evans, R. J., & Sejdinovic, D. (2024). Doubly robust
kernel statistics for testing distributional treatment effects.
*Transactions on Machine Learning Research*.

## See also

[`dr_date_test()`](https://max578.github.io/kernR/reference/dr_date_test.md)
for the underlying observational-causal test;
[`PESTO::pesto_ensemble_manifest`](https://rdrr.io/pkg/PESTO/man/pesto_ensemble_manifest.html)
for the input contract;
[`PESTO::pesto_ies_callback()`](https://rdrr.io/pkg/PESTO/man/pesto_ies_callback.html)
for producing the ensembles upstream.

Other distributional treatment effects:
[`dr_date_test()`](https://max578.github.io/kernR/reference/dr_date_test.md),
[`dr_dett_test()`](https://max578.github.io/kernR/reference/dr_dett_test.md)

## Examples

``` r
# \donttest{
# Requires PESTO (>= 0.4.1) -- wired through Imports.
library(PESTO)
npar <- 2L; nobs <- 4L; nreal <- 60L
G  <- matrix(stats::rnorm(nobs * npar), nobs, npar)
y0 <- as.numeric(G %*% c(1.0, -0.5)) + stats::rnorm(nobs, sd = 0.05)
y1 <- y0 + c(0.6, 0.6, 0.6, 0.6)   # intervention shifts outputs
names(y0) <- names(y1) <- paste0("o", seq_len(nobs))

prior <- matrix(stats::rnorm(nreal * npar), nreal, npar,
                dimnames = list(NULL, c("p1", "p2")))
fit0 <- pesto_ies_callback(function(t) t %*% t(G), prior, y0, 0.05,
                           noptmax = 3, verbose = FALSE)
fit1 <- pesto_ies_callback(function(t) t %*% t(G), prior, y1, 0.05,
                           noptmax = 3, verbose = FALSE)
m_base <- as_manifest(fit0, run_id = "baseline")
m_intv <- as_manifest(fit1, run_id = "intervention")
res <- dr_date_scenario(m_base, m_intv,
                         n_permutations = 200L, seed = 1L)
print(res)
#> 
#>    DR-DATE (scenario) Test
#> 
#> Statistic: 0.203779 
#> P-value:   1.0000 
#> N:         120 
#> Perms:     200 
#> Kernel Y:  rbf (bw = 0.8709)
#> ESS:       60.0 
#> 
#> Scenario contrast
#>   baseline      : baseline (n=60)
#>   intervention  : intervention (n=60)
#>   outputs tested: o1, o2, o3, o4
#>   PESTO versions: baseline=0.5.0, intervention=0.5.0
#>   fidelity      : baseline=single, intervention=single
#>   Verdict:        fail to reject (no distributional difference detected)
#> 
# }
```
