# Create a Kernel Specification

Constructs a kernel specification object used throughout `kernR` for
computing kernel matrices. Supports RBF (Gaussian), Matern, linear, and
polynomial kernels.

## Usage

``` r
kernel_spec(
  type = c("rbf", "matern", "linear", "polynomial"),
  bandwidth = "median",
  nu = 2.5,
  degree = 2L,
  offset = 1
)
```

## Arguments

- type:

  Character. Kernel type: `"rbf"` (default), `"matern"`, `"linear"`, or
  `"polynomial"`.

- bandwidth:

  Numeric or `"median"`. Lengthscale parameter for RBF and Matern
  kernels. If `"median"` (default), the median heuristic is used to
  select bandwidth automatically from the data.

- nu:

  Numeric. Smoothness parameter for the Matern kernel. Common choices:
  0.5 (Laplace), 1.5, 2.5, Inf (RBF). Default is 2.5.

- degree:

  Integer. Degree for polynomial kernel. Default is 2.

- offset:

  Numeric. Offset for polynomial kernel. Default is 1.

## Value

An object of class `"kernel_spec"`.

## See also

Other kernel primitives:
[`kernel_matrix()`](https://max578.github.io/kernR/reference/kernel_matrix.md),
[`resolve_bandwidth()`](https://max578.github.io/kernR/reference/resolve_bandwidth.md),
[`select_bandwidth()`](https://max578.github.io/kernR/reference/select_bandwidth.md),
[`weighted_hsic_stat()`](https://max578.github.io/kernR/reference/weighted_hsic_stat.md)

## Examples

``` r
# Default RBF kernel with median heuristic bandwidth
k <- kernel_spec()

# RBF with fixed bandwidth
k <- kernel_spec("rbf", bandwidth = 1.0)

# Matern kernel
k <- kernel_spec("matern", nu = 1.5)

# Linear kernel (no bandwidth needed)
k <- kernel_spec("linear")
```
