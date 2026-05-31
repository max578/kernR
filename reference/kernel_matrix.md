# Compute a Kernel Matrix

Computes the kernel (Gram) matrix between two sets of observations.

## Usage

``` r
kernel_matrix(x, y = NULL, kernel = kernel_spec())
```

## Arguments

- x:

  Numeric matrix (n x d) or vector.

- y:

  Numeric matrix (m x d) or vector. If `NULL` (default), computes the
  kernel matrix of `x` with itself.

- kernel:

  A `kernel_spec` object. Default is RBF with median heuristic
  bandwidth.

## Value

An n x m numeric matrix.

## See also

Other kernel primitives:
[`kernel_spec()`](https://max578.github.io/kernR/reference/kernel_spec.md),
[`select_bandwidth()`](https://max578.github.io/kernR/reference/select_bandwidth.md)

## Examples

``` r
x <- matrix(rnorm(100), 50, 2)
K <- kernel_matrix(x)
dim(K)  # 50 x 50
#> [1] 50 50
```
