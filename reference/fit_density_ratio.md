# Fit a Density-Ratio Model

Trains a density-ratio estimator for the do-null reweighting
`w(x, z) = p*(x) / p(x | z)` used by
[`bd_hsic_test()`](https://max578.github.io/kernR/reference/bd_hsic_test.md).
The fitted model is decoupled from evaluation:
[`predict_density_ratio()`](https://max578.github.io/kernR/reference/predict_density_ratio.md)
applies it to held-out rows, so train/test splits are honoured cleanly.

## Usage

``` r
fit_density_ratio(
  x,
  z,
  method = c("logistic", "ranger", "xgboost", "proxymix"),
  n_noise = 1L,
  proxymix_components = 2L,
  seed = NULL
)
```

## Arguments

- x:

  Numeric vector or matrix. Treatment variable (training).

- z:

  Numeric matrix or data.frame. Confounders (training).

- method:

  Character. Backend: `"logistic"` (default), `"ranger"`, `"xgboost"`,
  or `"proxymix"`.

- n_noise:

  Integer. Noise samples per real sample for classifier backends.
  Default `1L`.

- proxymix_components:

  Integer. Mixture components per density when `method = "proxymix"`.
  Default `2L`.

- seed:

  Integer or `NULL`. Random seed.

## Value

An object of class `density_ratio_fit` (plus
`density_ratio_fit_<method>` as the dispatch class). Carries: `method`,
the backend-specific fit (`model` for classifiers; `fit_joint` +
`fit_marg` for proxymix), `diagnostics`, `n_train`, `ncol_x`, `ncol_z`,
`seed`.

## Details

Four backends are supported (the `method` argument):

- `"logistic"` (default), `"ranger"`, `"xgboost"` – classifier-based
  noise-contrastive estimation. The classifier is trained to distinguish
  joint samples `(x, z)` from product-of-marginals samples
  `(x_perm, z)`; the density ratio is recovered from the calibrated
  class probabilities. Log-ratios are stored internally for numerical
  stability.

- `"proxymix"` – Gaussian-mixture density-ratio. Fits one GMM to the
  joint sample cloud `(x, z)` and one to a permuted product-of-marginals
  cloud via `proxymix::fit_proxymix(regime = "sample")`; ratios are
  evaluated in log-space from
  [`proxymix::dgmm()`](https://rdrr.io/pkg/proxymix/man/dgmm.html).
  Per-GMM convergence diagnostics (BIC, AIC, final log-likelihood,
  iteration count) are surfaced on the returned fit; query them via
  `fit$diagnostics`.

Introduced in kernR 0.0.0.9014 to close the documented-but-
unimplemented sample-split gap in
[`bd_hsic_test()`](https://max578.github.io/kernR/reference/bd_hsic_test.md)
(see NEWS).
[`estimate_density_ratio()`](https://max578.github.io/kernR/reference/estimate_density_ratio.md)
is now a thin backwards-compatible wrapper that fits and predicts on the
same data.

## See also

[`predict_density_ratio()`](https://max578.github.io/kernR/reference/predict_density_ratio.md),
[`estimate_density_ratio()`](https://max578.github.io/kernR/reference/estimate_density_ratio.md),
[`bd_hsic_test()`](https://max578.github.io/kernR/reference/bd_hsic_test.md).

Other density ratio and propensity:
[`assess_overlap()`](https://max578.github.io/kernR/reference/assess_overlap.md),
[`effective_sample_size()`](https://max578.github.io/kernR/reference/effective_sample_size.md),
[`estimate_density_ratio()`](https://max578.github.io/kernR/reference/estimate_density_ratio.md),
[`estimate_propensity()`](https://max578.github.io/kernR/reference/estimate_propensity.md),
[`plot_weights()`](https://max578.github.io/kernR/reference/plot_weights.md),
[`predict_density_ratio()`](https://max578.github.io/kernR/reference/predict_density_ratio.md)

## Examples

``` r
set.seed(1L)
n <- 200L
z <- matrix(rnorm(n * 2L), n, 2L)
x <- z[, 1L] + rnorm(n)
fit <- fit_density_ratio(x, z, method = "logistic", seed = 1L)
fit$diagnostics
#> $method
#> [1] "logistic"
#> 
#> $n_train
#> [1] 200
#> 
#> $n_noise
#> [1] 1
#> 
```
