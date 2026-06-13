# MMD Posterior-Predictive Check

Model-free verdict on whether a posterior-predictive ensemble is
consistent with held-out observations, via the Maximum Mean Discrepancy
(MMD) two-sample test. Wraps
[`mmd_test()`](https://max578.github.io/kernR/reference/mmd_test.md)
with a posterior-predictive framing and adds a Shannon-information
*surprise* diagnostic (`-log2(p)`) for intuitive interpretation: 0 bits
= no surprise (`p = 1`); ~4.32 bits = `p = 0.05`; the maximum achievable
surprise at `n_permutations = B` is `log2(B + 1)`.

## Usage

``` r
mmd_ppc(x, ...)

# Default S3 method
mmd_ppc(
  x,
  observed,
  kernel = kernel_spec(),
  n_permutations = 500L,
  alpha = 0.05,
  seed = NULL,
  ...
)

# S3 method for class 'pesto_ensemble'
mmd_ppc(x, observed = NULL, ...)

# S3 method for class 'pesto_ensemble_manifest'
mmd_ppc(x, observed, outputs = NULL, ...)
```

## Arguments

- x:

  Either a numeric matrix `M x d` of posterior-predictive draws, or a
  `pesto_ensemble` object (see
  [`pesto_ensemble()`](https://max578.github.io/kernR/reference/pesto_ensemble.md)).

- ...:

  Additional arguments (currently unused; reserved for future Nystrom
  acceleration).

- observed:

  Numeric matrix `n_obs x d` of held-out observations. When `x` is a
  `pesto_ensemble` carrying its own `observed` slot, may be left `NULL`
  to use the bundled observations.

- kernel:

  Kernel specification. Default is RBF with median heuristic over the
  pooled posterior + observed sample.

- n_permutations:

  Integer. Permutations for the null. Default 500.

- alpha:

  Numeric in `(0, 1)`. Significance level used for the verdict. Default
  `0.05`.

- seed:

  Integer or `NULL`. Random seed for reproducibility.

- outputs:

  Optional character vector of output column names from the manifest to
  test against. Defaults to all numeric output columns (the `real_name`
  column is excluded). Used only for the `pesto_ensemble_manifest`
  method.

## Value

An object of class `c("mmd_ppc", "kernel_test_result")` with the
standard `kernel_test_result` fields plus:

- n_posterior:

  Number of posterior-predictive draws.

- n_observed:

  Number of held-out observations.

- surprise_bits:

  Shannon-information surprise `-log2(p_value)`.

- alpha:

  Verdict significance level.

- reject:

  Logical: `p_value <= alpha`.

- pesto_metadata:

  Carried through from `pesto_ensemble` input, when provided; otherwise
  `NULL`.

The p-value lives on `result$p_value` (with an underscore), **not**
`p.value`. For a flat one-row summary with the `broom`-canonical
`p.value` column – plus `surprise_bits` and `reject` – call
[`generics::tidy()`](https://generics.r-lib.org/reference/tidy.html) on
the result.

## Details

Use after an ensemble-smoother run (PESTO IES, EnKF, etc.) to ask: *does
the calibrated model produce predictive draws that match the held-out
year / paddock / season at the distributional level?* The MMD test is
sensitive to mean, variance, and tail differences – strictly more
informative than RMSE on the posterior predictive mean.

The `pesto_ensemble_manifest` method records provenance from the input
manifest in `result$pesto_metadata`: `run_id`, `pesto_version`,
`method`, `outputs_used`, and `fidelity` (the multi-fidelity provenance
record, or `NULL` for a single-fidelity run), so the PPC verdict traces
back to the producing ensemble and the fidelity it was calibrated at.

## References

Gretton, A., Borgwardt, K. M., Rasch, M. J., Scholkopf, B., & Smola, A.
(2012). A kernel two-sample test. *JMLR*, 13, 723-773.

Gelman, A., Meng, X.-L., & Stern, H. (1996). Posterior predictive
assessment of model fitness via realized discrepancies. *Statistica
Sinica*, 6(4), 733-760.

## See also

[`mmd_test()`](https://max578.github.io/kernR/reference/mmd_test.md),
[`pesto_ensemble()`](https://max578.github.io/kernR/reference/pesto_ensemble.md)

Other posterior predictive checks:
[`pesto_ensemble()`](https://max578.github.io/kernR/reference/pesto_ensemble.md)

## Examples

``` r
set.seed(1)
# Calibrated model: posterior matches truth
post <- matrix(stats::rnorm(400L), ncol = 2L)
obs  <- matrix(stats::rnorm(40L),  ncol = 2L)
fit_ok <- mmd_ppc(post, obs, n_permutations = 199L, seed = 1L)
fit_ok
#> 
#>    MMD PPC Test
#> 
#> Statistic: -0.0113161 
#> P-value:   0.7600 
#> N:         220 
#> Perms:     199 
#> Kernel X:  rbf (bw = 1.619)
#> 
#> PPC verdict
#>   Posterior:  200 draws
#>   Observed:   20 obs
#>   Surprise:   0.396 bits
#>   Verdict:    consistent with observations
#> 

# Miscalibrated model: posterior is mean-shifted
obs_shift <- obs + 1.5
fit_bad <- mmd_ppc(post, obs_shift, n_permutations = 199L, seed = 1L)
fit_bad
#> 
#>    MMD PPC Test
#> 
#> Statistic: 0.373352 
#> P-value:   0.0050 
#> N:         220 
#> Perms:     199 
#> Kernel X:  rbf (bw =  1.72)
#> 
#> PPC verdict
#>   Posterior:  200 draws
#>   Observed:   20 obs
#>   Surprise:   7.644 bits
#>   Verdict:    REJECT (posterior inconsistent with observations)
#> 
```
