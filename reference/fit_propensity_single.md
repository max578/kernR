# Fit Propensity Model (Single Fold)

Fit Propensity Model (Single Fold)

## Usage

``` r
fit_propensity_single(treatment, covariates, method, newdata = NULL)
```

## Arguments

- treatment:

  Binary vector.

- covariates:

  Numeric matrix.

- method:

  Classification method.

- newdata:

  Optional matrix for prediction. If NULL, predicts on training data.

## Value

Vector of predicted probabilities.
