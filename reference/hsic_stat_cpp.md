# Compute the HSIC statistic (unweighted)

Biased HSIC estimator: (1/n^2) \* trace(K_x_c %\*% K_y_c) where K_c = H
K H with H = I - 1/n.

## Usage

``` r
hsic_stat_cpp(Kx, Ky)
```

## Arguments

- Kx:

  n x n kernel matrix for X.

- Ky:

  n x n kernel matrix for Y.

## Value

Scalar HSIC value.
