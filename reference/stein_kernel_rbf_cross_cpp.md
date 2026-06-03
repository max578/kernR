# RBF Stein cross-kernel matrix (n x m)

Rectangular counterpart of stein_kernel_rbf_cpp: the n x m block
u_p(x_i, z_j) of the Gaussian Stein kernel between the full sample X
(with score S) and landmark points Z (with score Sm). Bandwidth
convention matches stein_kernel_rbf_cpp.

## Usage

``` r
stein_kernel_rbf_cross_cpp(X, S, Z, Sm, h2)
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

- h2:

  Squared bandwidth h^2 (positive scalar).

## Value

n x m Stein cross-kernel matrix.
