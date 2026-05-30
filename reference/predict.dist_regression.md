# Predict from a Fitted Distribution Regression Model

Predict from a Fitted Distribution Regression Model

## Usage

``` r
# S3 method for class 'dist_regression'
predict(object, newdata, ...)
```

## Arguments

- object:

  A `dist_regression` fit.

- newdata:

  A list of new bags (matrices with the same number of columns as the
  training bags). Bag sizes may differ from training.

- ...:

  Currently ignored.

## Value

Numeric vector (length `length(newdata)`) when the model was fitted with
a scalar `y`; numeric matrix (`length(newdata) x d_y`) for multivariate
`y`.
