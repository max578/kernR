# Coverage / Calibration Diagnostic for a Predictive Ensemble

Quantifies *how* a predictive ensemble is calibrated against held-out
observations, rather than only testing whether it differs from them.
Each observation is mapped to its probability integral transform (PIT)
within the ensemble; calibrated draws give uniform PITs. The result
reports empirical coverage at nominal interval levels, a signed
dispersion ratio, a bias indicator, and a rank-histogram uniformity
test, and classifies the ensemble as calibrated, under-dispersed
(over-confident), over-dispersed, or biased.

## Usage

``` r
coverage_test(x, ...)

# Default S3 method
coverage_test(
  x,
  observed,
  levels = c(0.5, 0.8, 0.9),
  n_bins = 10L,
  alpha = 0.05,
  ...
)

# S3 method for class 'pesto_ensemble'
coverage_test(x, observed = NULL, ...)

# S3 method for class 'pesto_ensemble_manifest'
coverage_test(x, observed, outputs = NULL, ...)
```

## Arguments

- x:

  Numeric matrix `n_draws x d` of predictive draws, a `pesto_ensemble`
  (see
  [`pesto_ensemble()`](https://max578.github.io/kernR/reference/pesto_ensemble.md)),
  or a `pesto_ensemble_manifest`.

- ...:

  Additional arguments passed between methods (currently unused).

- observed:

  Numeric matrix `n_obs x d` of held-out observations. When `x` is a
  `pesto_ensemble` carrying an `observed` slot, may be `NULL`.

- levels:

  Numeric vector in `(0, 1)`. Central predictive-interval levels at
  which to report empirical coverage. Default `c(0.5, 0.8, 0.9)`.

- n_bins:

  Integer. Number of equal-width bins for the rank-histogram uniformity
  test. Default `10`.

- alpha:

  Numeric in `(0, 1)`. Significance level for the calibration verdict.
  Default `0.05`.

- outputs:

  Optional character vector of manifest output columns to test (the
  `real_name` column is always excluded). Defaults to all numeric output
  columns. Used only for the `pesto_ensemble_manifest` method.

## Value

An object of class `"coverage_test"` with components:

- coverage:

  Data frame of `nominal` vs pooled `empirical` coverage at each
  requested level.

- coverage_by_dim:

  Matrix of empirical coverage per dimension x level.

- dispersion_ratio:

  `Var(PIT) / (1/12)`: above one under-dispersed, below one
  over-dispersed.

- mean_pit:

  Pooled mean PIT (0.5 under calibration; a bias indicator).

- calibration:

  List with the rank-histogram chi-squared `statistic`, `p_value`, and
  `n_bins`.

- verdict:

  Character: the calibration classification.

- reject:

  Logical: `calibration$p_value <= alpha`.

- pit:

  The `n_obs x d` matrix of PIT values.

- n_draws, n_obs, dimension, levels, alpha:

  Inputs / sizes.

- pesto_metadata:

  Provenance carried from a `pesto_ensemble_manifest` (including
  `fidelity`), or `NULL`.

## Details

This is the graded complement to the binary kernel verdicts
[`mmd_ppc()`](https://max578.github.io/kernR/reference/mmd_ppc.md) and
[`ksd_test()`](https://max578.github.io/kernR/reference/ksd_test.md).
Its motivating use is ensemble over-confidence: an under-dispersed
ensemble (for example a collapsed iterative-ensemble-smoother posterior)
gives a U-shaped rank histogram, a dispersion ratio above one, and
empirical coverage below nominal – all of which this function names and
measures.

Calibration is assessed **per output dimension and pooled across
dimensions** (marginal calibration); it does not test the joint
dependence structure, for which the two-sample
[`mmd_ppc()`](https://max578.github.io/kernR/reference/mmd_ppc.md) is
the right tool. With few held-out observations the uniformity test has
low power; read the coverage and dispersion summaries (which are
informative at any sample size) alongside it.

## References

Gneiting, T., Balabdaoui, F., & Raftery, A. E. (2007). Probabilistic
forecasts, calibration and sharpness. *Journal of the Royal Statistical
Society B*, 69(2), 243-268.

Hamill, T. M. (2001). Interpretation of rank histograms for verifying
ensemble forecasts. *Monthly Weather Review*, 129(3), 550-560.

## See also

[`mmd_ppc()`](https://max578.github.io/kernR/reference/mmd_ppc.md),
[`ksd_test()`](https://max578.github.io/kernR/reference/ksd_test.md),
[`pesto_ensemble()`](https://max578.github.io/kernR/reference/pesto_ensemble.md)

Other goodness-of-fit tests:
[`concordance_test()`](https://max578.github.io/kernR/reference/concordance_test.md),
[`concordance_test_nystrom()`](https://max578.github.io/kernR/reference/concordance_test_nystrom.md),
[`gaussian_score()`](https://max578.github.io/kernR/reference/gaussian_score.md),
[`joint_coverage_test()`](https://max578.github.io/kernR/reference/joint_coverage_test.md),
[`ksd_test()`](https://max578.github.io/kernR/reference/ksd_test.md),
[`ksd_test_nystrom()`](https://max578.github.io/kernR/reference/ksd_test_nystrom.md),
[`numeric_score()`](https://max578.github.io/kernR/reference/numeric_score.md)

## Author

Max Moldovan, <max.moldovan@adelaide.edu.au>

## Examples

``` r
set.seed(1)
obs <- matrix(stats::rnorm(120L), ncol = 2L)

# Calibrated: predictive ensemble from the same law
ens_ok <- matrix(stats::rnorm(2000L), ncol = 2L)
coverage_test(ens_ok, obs)
#> 
#>   Coverage / calibration diagnostic
#> 
#> Ensemble:   1000 draws x 2 dims; 60 held-out obs
#> Coverage:   50%->62%   80%->85%   90%->94%  (nominal -> empirical)
#> Dispersion: ratio = 0.806 (Var(PIT)/(1/12); >1 under-, <1 over-dispersed)
#> Mean PIT:   0.539 (0.5 = unbiased)
#> Calibration: chi2 = 11.3, p = 0.2536 (rank histogram, 10 bins)
#> Verdict:    calibrated (PIT consistent with uniform)
#> 

# Under-dispersed (over-confident): ensemble too narrow
ens_tight <- matrix(stats::rnorm(2000L, sd = 0.4), ncol = 2L)
coverage_test(ens_tight, obs)
#> 
#>   Coverage / calibration diagnostic
#> 
#> Ensemble:   1000 draws x 2 dims; 60 held-out obs
#> Coverage:   50%->20%   80%->43%   90%->59%  (nominal -> empirical)
#> Dispersion: ratio = 1.851 (Var(PIT)/(1/12); >1 under-, <1 over-dispersed)
#> Mean PIT:   0.535 (0.5 = unbiased)
#> Calibration: chi2 =  111, p = 0.0000 (rank histogram, 10 bins)
#> Verdict:    REJECT -- under-dispersed (over-confident)
#> 
```
