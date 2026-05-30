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

## Examples

``` r
w <- runif(100, 0.5, 2)
effective_sample_size(w)
#> [1] 88.3534
```
