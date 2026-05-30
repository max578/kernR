# Compute DR-DETT Test Statistic (augmented IPW, effect on the treated)

Compute DR-DETT Test Statistic (augmented IPW, effect on the treated)

## Usage

``` r
compute_dr_dett_stat(Ky, treatment, e_hat, C0)
```

## Arguments

- Ky:

  n x n outcome kernel matrix.

- treatment:

  Binary treatment vector.

- e_hat:

  Propensity scores.

- C0:

  n x n control-arm conditional-mean-embedding coefficient matrix (zero
  matrix gives the inverse-probability-weighted statistic).

## Value

Scalar test statistic.
