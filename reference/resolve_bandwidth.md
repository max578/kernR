# Resolve Kernel Bandwidth

Resolves a `kernel_spec`'s bandwidth against data: when
`kernel$bandwidth` is `"median"` (the default for RBF and Matern
kernels) the median heuristic is computed from `x`; otherwise the fixed
bandwidth is returned unchanged. This is the helper
[`kernel_matrix()`](https://max578.github.io/kernR/reference/kernel_matrix.md)
uses internally, exposed so that callers building kernel matrices by
hand – for example a custom HSIC statistic via
[`weighted_hsic_stat()`](https://max578.github.io/kernR/reference/weighted_hsic_stat.md)
– resolve the bandwidth exactly as the package does.

## Usage

``` r
resolve_bandwidth(kernel, x)
```

## Arguments

- kernel:

  A `kernel_spec` object.

- x:

  Numeric matrix (`n` by `d`) of data used to compute the median
  heuristic when `kernel$bandwidth` is `"median"`.

## Value

A `kernel_spec` with a resolved numeric bandwidth.

## See also

[`kernel_spec()`](https://max578.github.io/kernR/reference/kernel_spec.md),
[`kernel_matrix()`](https://max578.github.io/kernR/reference/kernel_matrix.md),
[`weighted_hsic_stat()`](https://max578.github.io/kernR/reference/weighted_hsic_stat.md)

Other kernel primitives:
[`kernel_matrix()`](https://max578.github.io/kernR/reference/kernel_matrix.md),
[`kernel_spec()`](https://max578.github.io/kernR/reference/kernel_spec.md),
[`select_bandwidth()`](https://max578.github.io/kernR/reference/select_bandwidth.md),
[`weighted_hsic_stat()`](https://max578.github.io/kernR/reference/weighted_hsic_stat.md)

## Examples

``` r
x <- matrix(rnorm(40), ncol = 2)
resolve_bandwidth(kernel_spec(), x)
#> Kernel specification:
#>   Type: rbf 
#>   Bandwidth: 1.914354 
```
