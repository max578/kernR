# RBF (Gaussian) kernel matrix

K(x, y) = exp(-\|\|x - y\|\|^2 / (2 \* bandwidth^2))

## Usage

``` r
rbf_kernel_matrix_cpp(x, y, bandwidth)
```

## Arguments

- x:

  Numeric matrix (n x d).

- y:

  Numeric matrix (m x d).

- bandwidth:

  Positive scalar.

## Value

n x m kernel matrix.
