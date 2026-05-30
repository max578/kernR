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
