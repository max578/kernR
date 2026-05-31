# Print a Conditional Mean Embedding Fit

Prints a compact summary of a fitted conditional mean embedding: the
training-sample size, the input and output dimensions, the resolved
kernels, and the ridge regularisation parameter that was used.

## Usage

``` r
# S3 method for class 'cme_fit'
print(x, ...)
```

## Arguments

- x:

  A `cme_fit` object returned by
  [`fit_cme()`](https://max578.github.io/kernR/reference/fit_cme.md).

- ...:

  Unused; present for S3 generic compatibility.

## Value

The `cme_fit` object `x`, invisibly.

## See also

[`fit_cme()`](https://max578.github.io/kernR/reference/fit_cme.md)

## Examples

``` r
set.seed(1L)
x <- matrix(rnorm(60L), ncol = 2L)
y <- matrix(x[, 1L] + rnorm(30L, sd = 0.2), ncol = 1L)
print(fit_cme(x, y, lambda = 1e-2))
#> Conditional mean embedding (kernel ridge regression)
#>   Training points: 30
#>   Input dim:       2
#>   Output dim:      30
#>   Kernel (x):      rbf
#>   Kernel (y):      rbf
#>   Ridge lambda:     0.01
```
