# RBF Stein kernel matrix

Builds the n x n Stein-kernel matrix u_p(x_i, x_j) for the Gaussian base
kernel k(x, y) = exp(-\|\|x - y\|\|^2 / (2 h^2)) under the Langevin
Stein operator, given the score evaluated at each sample point. The
bandwidth convention matches rbf_kernel_matrix_cpp.

## Usage

``` r
stein_kernel_rbf_cpp(X, S, h2)
```

## Arguments

- X:

  Numeric matrix (n x d): the sample.

- S:

  Numeric matrix (n x d): the score evaluated row-wise at X.

- h2:

  Squared bandwidth h^2 (positive scalar).

## Value

Symmetric n x n Stein-kernel matrix.
