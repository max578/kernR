# Joint (Multivariate) Calibration Diagnostic for a Predictive Ensemble

Tests whether a predictive ensemble is calibrated *jointly*, including
its dependence structure, against held-out observations. Where
[`coverage_test()`](https://max578.github.io/kernR/reference/coverage_test.md)
assesses each output dimension separately (marginal calibration), this
function builds a multivariate rank histogram: each observation is
reduced, together with the ensemble, to one multivariate rank through a
pre-rank function, and the histogram of those ranks is tested for
uniformity. An ensemble with correct margins but a mis-specified
correlation – which
[`coverage_test()`](https://max578.github.io/kernR/reference/coverage_test.md)
passes – is caught here.

## Usage

``` r
joint_coverage_test(x, ...)

# Default S3 method
joint_coverage_test(
  x,
  observed,
  prerank = c("band_depth", "average"),
  levels = c(0.5, 0.8, 0.9),
  n_bins = 10L,
  alpha = 0.05,
  seed = NULL,
  ...
)

# S3 method for class 'pesto_ensemble'
joint_coverage_test(x, observed = NULL, ...)

# S3 method for class 'pesto_ensemble_manifest'
joint_coverage_test(x, observed, outputs = NULL, ...)
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

- prerank:

  Character. Pre-rank function: `"band_depth"` (default,
  dependence-sensitive) or `"average"` (familiar dispersion reading).

- levels:

  Numeric vector in `(0, 1)`. Central predictive-interval levels at
  which to report empirical coverage. Default `c(0.5, 0.8, 0.9)`.

- n_bins:

  Integer. Number of equal-width bins for the rank-histogram uniformity
  test. Default `10`.

- alpha:

  Numeric in `(0, 1)`. Significance level for the calibration verdict.
  Default `0.05`.

- seed:

  Integer or `NULL`. Random seed: multivariate ranks break pre-rank ties
  at random, so a non-`NULL` seed makes the result reproducible.

- outputs:

  Optional character vector of manifest output columns to test (the
  `real_name` column is always excluded). Defaults to all numeric output
  columns. Used only for the `pesto_ensemble_manifest` method.

## Value

An object of class `"joint_coverage_test"` with components:

- prerank:

  The pre-rank function used.

- coverage:

  Data frame of `nominal` vs `empirical` central-rank-band coverage at
  each requested level.

- dispersion_ratio:

  `Var(u) / (1/12)` of the normalised multivariate ranks `u`; above one
  under-dispersed, below one over-dispersed (the primary signal for
  `prerank = "average"`).

- mean_rank:

  Mean normalised rank (`0.5` under calibration; the primary direction
  signal for `prerank = "band_depth"`).

- calibration:

  List with the rank-histogram chi-squared `statistic`, `p_value`, and
  `n_bins`.

- verdict:

  Character: the joint-calibration classification.

- reject:

  Logical: `calibration$p_value <= alpha`.

- ranks:

  Integer multivariate ranks, one per observation.

- n_draws, n_obs, dimension, levels, alpha:

  Inputs / sizes.

- pesto_metadata:

  Provenance carried from a `pesto_ensemble_manifest` (including
  `fidelity`), or `NULL`.

## Details

The pre-rank function determines what miscalibration is visible:

- `"band_depth"` (default):

  Ranks points by multivariate centrality (band depth; Thorarinsdottir
  et al. 2016). Sensitive to dependence / correlation miscalibration
  that marginal and average-rank methods miss. The reading is a *slope*:
  a mean rank below `0.5` means observations fall outside the ensemble
  cloud (jointly under-dispersed); above `0.5` means they sit too
  centrally (over-dispersed); a non-uniform histogram with a central
  mean signals a dependence error that is not a pure dispersion shift.

- `"average"`:

  Ranks points by the sum of their per-dimension ranks (Gneiting et al.
  2008). Gives the familiar rank-histogram reading – a U-shape and a
  dispersion ratio above one mean under-dispersion, an inverted-U and a
  ratio below one mean over-dispersion – but is weaker against
  correlation-only errors.

Under calibration the multivariate rank is uniform on
`{1, ..., n_draws + 1}` for either pre-rank, so the chi-squared
uniformity test, the dispersion ratio, and the coverage table are all
referenced against the uniform. With few held-out observations the
uniformity test has low power; read the dispersion and mean-rank
summaries alongside it. A joint test needs at least two output
dimensions; for a single output use
[`coverage_test()`](https://max578.github.io/kernR/reference/coverage_test.md).

## References

Gneiting, T., Stanberry, L. I., Grimit, E. P., Held, L., & Johnson, N.
A. (2008). Assessing probabilistic forecasts of multivariate quantities,
with an application to ensemble predictions of surface winds. *TEST*,
17(2), 211-235.

Thorarinsdottir, T. L., Scheuerer, M., & Heinz, C. (2016). Assessing the
calibration of high-dimensional ensemble forecasts using rank
histograms. *Journal of Computational and Graphical Statistics*, 25(1),
105-122.

## See also

[`coverage_test()`](https://max578.github.io/kernR/reference/coverage_test.md),
[`mmd_ppc()`](https://max578.github.io/kernR/reference/mmd_ppc.md),
[`concordance_test()`](https://max578.github.io/kernR/reference/concordance_test.md)

Other goodness-of-fit tests:
[`concordance_test()`](https://max578.github.io/kernR/reference/concordance_test.md),
[`concordance_test_nystrom()`](https://max578.github.io/kernR/reference/concordance_test_nystrom.md),
[`coverage_test()`](https://max578.github.io/kernR/reference/coverage_test.md),
[`gaussian_score()`](https://max578.github.io/kernR/reference/gaussian_score.md),
[`ksd_test()`](https://max578.github.io/kernR/reference/ksd_test.md),
[`ksd_test_nystrom()`](https://max578.github.io/kernR/reference/ksd_test_nystrom.md),
[`numeric_score()`](https://max578.github.io/kernR/reference/numeric_score.md)

## Author

Max Moldovan, <max.moldovan@adelaide.edu.au>

## Examples

``` r
set.seed(1)
chol2 <- chol(matrix(c(1, 0.9, 0.9, 1), 2L))

# Correctly correlated ensemble: calibrated jointly
obs <- matrix(stats::rnorm(120L), ncol = 2L) %*% chol2
ens_ok <- matrix(stats::rnorm(4000L), ncol = 2L) %*% chol2
joint_coverage_test(ens_ok, obs, seed = 1L)
#> 
#>   Joint (multivariate) calibration diagnostic
#> 
#> Ensemble:    2000 draws x 2 dims; 60 held-out obs
#> Pre-rank:    band_depth (dependence-sensitive)
#> Coverage:    50%->53%   80%->80%   90%->88%  (nominal -> empirical rank band)
#> Mean rank:   0.581 (0.5 = calibrated; <0.5 obs outlying, >0.5 obs central)
#> Dispersion:  ratio = 0.883 (Var(rank)/(1/12))
#> Calibration: chi2 =    8, p = 0.5341 (rank histogram, 10 bins)
#> Verdict:     calibrated (multivariate ranks consistent with uniform)
#> 

# Right margins, wrong dependence: ensemble independent, obs correlated.
# coverage_test() passes; joint_coverage_test() catches it.
ens_indep <- matrix(stats::rnorm(4000L), ncol = 2L)
joint_coverage_test(ens_indep, obs, seed = 1L)
#> 
#>   Joint (multivariate) calibration diagnostic
#> 
#> Ensemble:    2000 draws x 2 dims; 60 held-out obs
#> Pre-rank:    band_depth (dependence-sensitive)
#> Coverage:    50%->37%   80%->62%   90%->77%  (nominal -> empirical rank band)
#> Mean rank:   0.622 (0.5 = calibrated; <0.5 obs outlying, >0.5 obs central)
#> Dispersion:  ratio = 1.225 (Var(rank)/(1/12))
#> Calibration: chi2 =   25, p = 0.0030 (rank histogram, 10 bins)
#> Verdict:     REJECT -- jointly over-dispersed (observations too central)
#> 
```
