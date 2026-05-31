# Compute Effective Sample Size

Computes ESS from importance weights: ESS = (sum(w))^2 / sum(w^2).

## Usage

``` r
effective_sample_size(w)
```

## Arguments

- w:

  Numeric vector of weights.

## Value

Scalar effective sample size.

## See also

Other density ratio and propensity:
[`assess_overlap()`](https://max578.github.io/kernR/reference/assess_overlap.md),
[`estimate_density_ratio()`](https://max578.github.io/kernR/reference/estimate_density_ratio.md),
[`estimate_propensity()`](https://max578.github.io/kernR/reference/estimate_propensity.md),
[`fit_density_ratio()`](https://max578.github.io/kernR/reference/fit_density_ratio.md),
[`plot_weights()`](https://max578.github.io/kernR/reference/plot_weights.md),
[`predict_density_ratio()`](https://max578.github.io/kernR/reference/predict_density_ratio.md)

## Examples

``` r
w <- runif(100, 0.5, 2)
effective_sample_size(w)
#> [1] 88.3534
```
