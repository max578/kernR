# Compute the weighted HSIC statistic (for bd-HSIC)

Weighted version \$\$\sum\_{i,j} w_i w_j (K_x^c)\_{ij}
(K_y^c)\_{ij}\$\$.

## Usage

``` r
weighted_hsic_stat_cpp(Kx, Ky, w)
```

## Arguments

- Kx:

  n x n kernel matrix for X.

- Ky:

  n x n kernel matrix for Y.

- w:

  Weight vector of length n.

## Value

Scalar weighted HSIC value.
