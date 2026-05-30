# Estimate Density Ratios (backwards-compatible wrapper)

Wraps
[`fit_density_ratio()`](https://max578.github.io/kernR/reference/fit_density_ratio.md) +
[`predict_density_ratio()`](https://max578.github.io/kernR/reference/predict_density_ratio.md)
on the same data. Preserved for backwards compatibility with kernR
0.0.0.901x callers; new code should prefer the explicit fit/predict pair
so train/test splits are honoured.

## Usage

``` r
estimate_density_ratio(
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

A list of class `density_ratio_fit_estimate` with components `weights`,
`ratios`, `ess`, `method`, `n`, and `fit` (the underlying
`density_ratio_fit` for callers that want diagnostics).

## Details

The return shape (`weights`, `ratios`, `ess`, `method`, `n`) is
unchanged from previous versions. Internally, ratios are now computed in
log-space (which fixes pathological tail behaviour that the
classifier-based 0.0.0.9012 implementation occasionally showed under
extreme imbalance).

## See also

[`fit_density_ratio()`](https://max578.github.io/kernR/reference/fit_density_ratio.md)
for the fit/predict surface.

## Examples

``` r
set.seed(42)
n <- 200
z <- matrix(rnorm(n * 2), n, 2)
x <- z[, 1] + rnorm(n)
dr <- estimate_density_ratio(x, z)
dr$ess
#> [1] 198.4936
```
