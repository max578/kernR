# Plot Weight Diagnostics

Plots the distribution of importance weights with effective sample size
annotation.

## Usage

``` r
plot_weights(weights, main = "Weight Distribution")
```

## Arguments

- weights:

  Numeric vector of importance weights.

- main:

  Title. Default is "Weight Distribution".

## Value

Invisibly returns `weights`.

## See also

Other density ratio and propensity:
[`assess_overlap()`](https://max578.github.io/kernR/reference/assess_overlap.md),
[`effective_sample_size()`](https://max578.github.io/kernR/reference/effective_sample_size.md),
[`estimate_density_ratio()`](https://max578.github.io/kernR/reference/estimate_density_ratio.md),
[`estimate_propensity()`](https://max578.github.io/kernR/reference/estimate_propensity.md),
[`fit_density_ratio()`](https://max578.github.io/kernR/reference/fit_density_ratio.md),
[`predict_density_ratio()`](https://max578.github.io/kernR/reference/predict_density_ratio.md)

## Examples

``` r
set.seed(1L)
weights <- rgamma(200L, shape = 2, rate = 2)
plot_weights(weights)

```
