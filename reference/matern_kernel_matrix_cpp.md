# Matern kernel matrix

Supports nu = 0.5, 1.5, 2.5, and Inf (RBF).

## Usage

``` r
matern_kernel_matrix_cpp(x, y, bandwidth, nu)
```

## Arguments

- x:

  Numeric matrix (n x d).

- y:

  Numeric matrix (m x d).

- bandwidth:

  Positive scalar.

- nu:

  Smoothness parameter.

## Value

n x m kernel matrix.
