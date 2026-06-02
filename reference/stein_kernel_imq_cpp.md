# Inverse multi-quadric Stein kernel matrix

Builds the n x n Stein-kernel matrix u_p(x_i, x_j) for the IMQ base
kernel k(x, y) = (c^2 + \|\|x - y\|\|^2)^beta under the Langevin Stein
operator, given the score (gradient of the log target density) evaluated
at each sample point.

## Usage

``` r
stein_kernel_imq_cpp(X, S, beta, c2)
```

## Arguments

- X:

  Numeric matrix (n x d): the sample.

- S:

  Numeric matrix (n x d): the score evaluated row-wise at X.

- beta:

  Negative scalar exponent in (-1, 0).

- c2:

  Squared offset c^2 (positive scalar).

## Value

Symmetric n x n Stein-kernel matrix.
