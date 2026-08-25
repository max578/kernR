# Predict Conditional Mean Embedding Weights at New Points

Returns the row weights \\\alpha(x^\*) = k(x^\*, X\_{\mathrm{train}})
W\\ from a fitted
[`fit_cme()`](https://max578.github.io/kernR/reference/fit_cme.md)
object. Each row of the result is the weight vector for combining
training-`y` quantities (kernel values or values themselves) to produce
a CME prediction at `x_new`.

## Usage

``` r
# S3 method for class 'cme_fit'
predict(object, x_new, ...)
```

## Arguments

- object:

  A `cme_fit` object.

- x_new:

  Numeric matrix of new conditioning points.

- ...:

  Currently ignored.

## Value

An n_new x n_train matrix of embedding weights.

## Details

For a typical "predict Y at new X" workflow use
[`kernel_downscale()`](https://max578.github.io/kernR/reference/kernel_downscale.md),
which combines this with the training Y matrix to return predictions
directly.

## Examples

``` r
set.seed(1L)
x <- matrix(rnorm(60L), ncol = 2L)
y <- matrix(x[, 1L] + rnorm(30L, sd = 0.2), ncol = 1L)
fit <- fit_cme(x, y, lambda = 1e-2)
x_new <- matrix(rnorm(6L), ncol = 2L)
w <- predict(fit, x_new)
dim(w)  # 3 new points x 30 training points
#> [1]  3 30
```
