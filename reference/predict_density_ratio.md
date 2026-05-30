# Predict from a Fitted Density-Ratio Model

Applies a `density_ratio_fit` object (from
[`fit_density_ratio()`](https://max578.github.io/kernR/reference/fit_density_ratio.md))
to new `(x, z)` rows. All four backends compute ratios in log-space
internally for numerical stability; `type` controls the returned
representation.

## Usage

``` r
predict_density_ratio(
  object,
  new_x,
  new_z,
  type = c("log_ratio", "weight", "ratio")
)
```

## Arguments

- object:

  A `density_ratio_fit`.

- new_x:

  Numeric vector or matrix. Treatment values to evaluate.

- new_z:

  Numeric matrix or data.frame. Confounders to evaluate.

- type:

  Character. Return type: `"log_ratio"` – natural-log density ratio
  (default; preferred for downstream calculation); `"ratio"` – raw
  density ratio (`exp(log_ratio)`); `"weight"` – IES-compatible
  normalised weights (positive, sum-to-`n_new`).

## Value

Numeric vector of length `nrow(new_x)`.

## See also

[`fit_density_ratio()`](https://max578.github.io/kernR/reference/fit_density_ratio.md).

## Examples

``` r
set.seed(1L)
n <- 200L
z <- matrix(rnorm(n * 2L), n, 2L)
x <- z[, 1L] + rnorm(n)
fit <- fit_density_ratio(x, z, method = "logistic", seed = 1L)
weights <- predict_density_ratio(fit, new_x = x, new_z = z,
                                 type = "weight")
summary(weights)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.8512  0.9721  1.0000  1.0000  1.0270  1.1256 
```
