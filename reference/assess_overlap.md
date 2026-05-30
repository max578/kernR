# Assess Propensity Score Overlap

Diagnoses overlap (positivity) between treated and control groups by
summarising the propensity score distributions.

## Usage

``` r
assess_overlap(propensity, treatment = NULL)
```

## Arguments

- propensity:

  A `propensity_fit` object or a numeric vector of scores.

- treatment:

  Binary treatment vector. Required if `propensity` is a numeric vector.

## Value

A list of class `"overlap_diagnostic"` with:

- treated:

  Summary statistics of propensity scores for treated.

- control:

  Summary statistics for controls.

- overlap_warning:

  Logical. TRUE if overlap is poor.

## Examples

``` r
set.seed(1L)
n <- 200L
treatment <- rbinom(n, 1L, 0.5)
scores <- plogis(rnorm(n) + 0.6 * treatment)
assess_overlap(scores, treatment)
#> Propensity Score Overlap Diagnostic
#> 
#> Treated:   min = 0.141, q25 = 0.484, median = 0.627, q75 = 0.780, max = 0.914 
#> Control:   min = 0.053, q25 = 0.368, median = 0.482, q75 = 0.651, max = 0.934 
#> Overlap:   87.7 %
```
