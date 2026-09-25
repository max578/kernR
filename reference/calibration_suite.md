# Calibration suite for a continuous prediction

Quantifies whether predictions are right in *level*, not merely in
order. A model can rank cases correctly while being systematically wrong
about magnitude, and it is magnitude that a loss-optimal decision
prices.

## Usage

``` r
calibration_suite(x, ...)

# Default S3 method
calibration_suite(
  x,
  observed,
  n_boot = 1000L,
  tolerance = 0.1,
  conf_level = 0.95,
  seed = NULL,
  ...
)
```

## Arguments

- x:

  Numeric vector of predictions, or an object with a method.

- ...:

  Passed to methods.

- observed:

  Numeric vector of realised outcomes, the same length as `x`.

- n_boot:

  Bootstrap replicates for the percentile intervals.

- tolerance:

  The deviation from a slope of 1 the caller would want to detect.
  Drives the power-aware verdict; it does not affect the estimates.

- conf_level:

  Interval level for every reported statistic.

- seed:

  Optional integer seed for the bootstrap, for reproducibility.

## Value

An object of class `calibration_suite` (additionally `kernR_abstention`
when the verdict is `"undetermined"`), a list with `slope`, `intercept`,
`citl` and `ici` (each a named numeric of `estimate`, `lower`, `upper`),
`curve` (a data.frame for plotting), `detectable`, `verdict`,
`abstained`, and the inputs `n`, `n_boot`, `tolerance`, `conf_level`.

## Details

Four statistics, each with a bootstrap percentile interval:

- slope:

  Coefficient of `observed ~ predicted`. Ideal 1. Below 1 means
  predictions are too spread out – the classic over-fitting signature,
  and the quantity that collapses first when a model moves to a new
  site.

- intercept:

  Intercept of `observed ~ offset(predicted)`. Ideal 0.

- citl:

  Calibration-in-the-large, `mean(observed) - mean(predicted)`. Ideal 0.
  The pure level shift, free of any slope effect.

- ici:

  Integrated calibration index: the mean absolute distance between a
  smooth of the observed outcome on the prediction and the prediction
  itself. In the units of the outcome.

## The verdict is power-aware

`verdict` is `"miscalibrated"` when the slope interval excludes 1;
`"calibrated"` when it contains 1 *and* is narrow enough to have
resolved a deviation of `tolerance`; and `"undetermined"` otherwise –
the test could not have detected the miscalibration the caller cares
about, so it declines to certify calibration. An `"undetermined"` result
is stamped `kernR_abstention`, so a caller gating on calibration cannot
mistake an under-powered pass for a real one. `detectable` reports the
smallest deviation from a slope of 1 that this sample could resolve.

## References

Van Calster, B., McLernon, D. J., van Smeden, M., Wynants, L., &
Steyerberg, E. W. (2019). Calibration: the Achilles heel of predictive
analytics. *BMC Medicine*, 17, 230.

Austin, P. C., & Steyerberg, E. W. (2019). The Integrated Calibration
Index (ICI) and related metrics for quantifying the calibration of
logistic regression models. *Statistics in Medicine*, 38(21), 4051-4065.

## See also

[`coverage_test()`](https://max578.github.io/kernR/reference/coverage_test.md)
for the width of a predictive distribution.

Other goodness-of-fit tests:
[`concordance_test()`](https://max578.github.io/kernR/reference/concordance_test.md),
[`concordance_test_nystrom()`](https://max578.github.io/kernR/reference/concordance_test_nystrom.md),
[`coverage_test()`](https://max578.github.io/kernR/reference/coverage_test.md),
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
# A calibrated prediction is one the outcome varies around: generate the
# prediction first, then the outcome from it. Note that `truth + noise` is
# NOT calibrated as a prediction of `truth` -- it carries noise the outcome
# does not, so its slope attenuates below 1 by construction.
pred <- stats::rnorm(300L, mean = 5, sd = 2)
obs <- pred + stats::rnorm(300L, sd = 0.5)
calibration_suite(pred, obs, n_boot = 200L)
#> 
#>   Calibration suite (continuous outcome)
#> 
#>   n = 300 observations; 200 bootstrap replicates (200 usable)
#>   slope           1.006  [  0.978,   1.041]   (ideal 1)
#>   intercept      -0.005  [ -0.070,   0.048]   (ideal 0)
#>   CITL           -0.005  [ -0.070,   0.048]   (ideal 0)
#>   ICI             0.040  [  0.025,   0.094]   (ideal 0, outcome units)
#> 
#>   Detectable:   0.032  (smallest slope departure this n resolves)
#>   Tolerance:    0.100  (what the caller asked to detect)
#>   Verdict:      calibrated
#> 

# Over-spread predictions: slope below 1
calibration_suite(5 + 1.6 * (pred - 5), pred, n_boot = 200L)
#> 
#>   Calibration suite (continuous outcome)
#> 
#>   n = 300 observations; 200 bootstrap replicates (200 usable)
#>   slope           0.625  [  0.625,   0.625]   (ideal 1)
#>   intercept      -0.040  [ -0.160,   0.085]   (ideal 0)
#>   CITL           -0.040  [ -0.160,   0.085]   (ideal 0)
#>   ICI             0.912  [  0.825,   0.985]   (ideal 0, outcome units)
#> 
#>   Detectable:   0.000  (smallest slope departure this n resolves)
#>   Tolerance:    0.100  (what the caller asked to detect)
#>   Verdict:      miscalibrated
#> 

# A noisy handful cannot certify anything: verdict is "undetermined"
calibration_suite(pred[1:6], pred[1:6] + stats::rnorm(6L, sd = 3),
                  n_boot = 200L)
#> 
#>   Calibration suite (continuous outcome)
#> 
#>   n = 6 observations; 200 bootstrap replicates (200 usable)
#>   slope           0.525  [ -0.516,   2.703]   (ideal 1)
#>   intercept       2.147  [  0.499,   4.067]   (ideal 0)
#>   CITL            2.147  [  0.499,   4.067]   (ideal 0)
#>   ICI                NA  [     NA,      NA]   (ideal 0, outcome units)
#> 
#>   Detectable:   1.609  (smallest slope departure this n resolves)
#>   Tolerance:    0.100  (what the caller asked to detect)
#>   Verdict:      undetermined
#>                 the interval is wider than the tolerance, so this
#>                 sample cannot certify calibration -- not a pass
#> 
```
