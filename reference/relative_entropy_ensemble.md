# Relative Entropy Between Two Sample Ensembles (Gaussian approximation)

Convenience wrapper around
[`relative_entropy()`](https://max578.github.io/kernR/reference/relative_entropy.md)
for ensemble inputs (for example a PESTO posterior ensemble): each
ensemble is summarised by its sample mean and covariance and the
Gaussian relative entropy is returned. This is a moment-matched
approximation – exact when the ensembles are Gaussian and a second-order
surrogate otherwise.

## Usage

``` r
relative_entropy_ensemble(p_draws, q_draws)
```

## Arguments

- p_draws:

  Numeric matrix `n_p x k` (or vector for `k = 1`): draws from `p` (the
  smoother ensemble, for ACI).

- q_draws:

  Numeric matrix `n_q x k` (or vector): draws from `q` (the filter
  ensemble). Must have the same number of columns as `p_draws`.

## Value

A single non-negative numeric: the moment-matched relative entropy.

## See also

[`relative_entropy()`](https://max578.github.io/kernR/reference/relative_entropy.md)

## Examples

``` r
set.seed(1)
p <- matrix(rnorm(2000), ncol = 2)
q <- matrix(rnorm(2000), ncol = 2) + 1
relative_entropy_ensemble(p, q)
#> [1] 0.9801008
```
