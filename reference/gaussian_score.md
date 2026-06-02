# Score function for a multivariate normal target

Builds a score function – the gradient of the log density, \\\nabla_x
\log p(x) = -\Sigma^{-1}(x - \mu)\\ – for a multivariate normal target,
suitable for the `score` argument of
[`ksd_test()`](https://max578.github.io/kernR/reference/ksd_test.md).

## Usage

``` r
gaussian_score(mean = NULL, sigma = NULL)
```

## Arguments

- mean:

  Numeric vector of length `d`, the target mean. `NULL` (default) uses
  the zero vector of the sample dimension.

- sigma:

  Numeric `d x d` covariance matrix. `NULL` (default) uses the identity
  matrix of the sample dimension. Must be symmetric and invertible.

## Value

A function of one argument (an `n x d` numeric matrix) returning the
`n x d` matrix of scores. The mean and covariance are captured by the
closure.

## Details

The returned closure accepts the sample matrix and returns the score
evaluated row-wise. Leaving `mean` or `sigma` at `NULL` defaults them to
the zero vector and the identity matrix of the dimension seen at call
time, so `gaussian_score()` with no arguments is the standard-normal
score in any dimension.

## References

Liu, Q., Lee, J. D., & Jordan, M. I. (2016). A kernelized Stein
discrepancy for goodness-of-fit tests. *Proceedings of the 33rd
International Conference on Machine Learning*, PMLR 48, 276-284.

## See also

[`ksd_test()`](https://max578.github.io/kernR/reference/ksd_test.md)

Other goodness-of-fit tests:
[`concordance_test()`](https://max578.github.io/kernR/reference/concordance_test.md),
[`coverage_test()`](https://max578.github.io/kernR/reference/coverage_test.md),
[`ksd_test()`](https://max578.github.io/kernR/reference/ksd_test.md),
[`numeric_score()`](https://max578.github.io/kernR/reference/numeric_score.md)

## Author

Max Moldovan, <max.moldovan@adelaide.edu.au>

## Examples

``` r
set.seed(1)
x <- matrix(stats::rnorm(200L), ncol = 2L)

# Standard-normal target in two dimensions
s0 <- gaussian_score()
str(s0(x))
#>  num [1:100, 1:2] 0.626 -0.184 0.836 -1.595 -0.33 ...

# Correlated-normal target
sig <- matrix(c(1, 0.5, 0.5, 1), nrow = 2L)
s1 <- gaussian_score(mean = c(0, 0), sigma = sig)
```
