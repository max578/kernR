# Cross-Validate Ridge Parameter for KRR

Simple LOO-CV for the ridge parameter in kernel ridge regression. Tests
a grid of lambda values and picks the one minimising LOO error.

## Usage

``` r
cv_ridge_lambda(Kx, Ky)
```

## Arguments

- Kx:

  n x n kernel matrix.

- Ky:

  n x n kernel matrix.

## Value

Optimal lambda (numeric scalar).
