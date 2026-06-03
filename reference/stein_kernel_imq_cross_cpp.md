# Inverse multi-quadric Stein cross-kernel matrix (n x m)

Rectangular counterpart of stein_kernel_imq_cpp: builds the n x m block
u_p(x_i, z_j) of the IMQ Stein kernel between the full sample X (with
score S) and a set of landmark points Z (with score Sm). Used to
assemble the Nystrom factorisation of the Stein-kernel matrix for
ksd_test_nystrom(). The per-pair formula is identical to the symmetric
builder; only the index ranges differ (rows over X, columns over Z).

## Usage

``` r
stein_kernel_imq_cross_cpp(X, S, Z, Sm, beta, c2)
```

## Arguments

- X:

  Numeric matrix (n x d): the full sample.

- S:

  Numeric matrix (n x d): the score evaluated row-wise at X.

- Z:

  Numeric matrix (m x d): the landmark points.

- Sm:

  Numeric matrix (m x d): the score evaluated row-wise at Z.

- beta:

  Negative scalar exponent in (-1, 0).

- c2:

  Squared offset c^2 (positive scalar).

## Value

n x m Stein cross-kernel matrix.
