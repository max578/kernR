# Select Kernel Bandwidth

Computes a bandwidth (lengthscale) for RBF or Matern kernels using the
specified method.

## Usage

``` r
select_bandwidth(x, method = "median")
```

## Arguments

- x:

  Numeric matrix or vector.

- method:

  Character. `"median"` (default), `"scott"`, or a positive number for
  fixed bandwidth.

## Value

A positive numeric scalar.

## Details

- `"median"`: The median heuristic sets bandwidth to the square root of
  the median of pairwise squared distances. Robust default for most
  kernel tests (Gretton et al., 2012).

- `"scott"`: Scott's rule: `n^(-1/(d+4)) * sd_pooled`. Good for density
  estimation but may undersmooth for testing.

## See also

Other kernel primitives:
[`kernel_matrix()`](https://max578.github.io/kernR/reference/kernel_matrix.md),
[`kernel_spec()`](https://max578.github.io/kernR/reference/kernel_spec.md),
[`resolve_bandwidth()`](https://max578.github.io/kernR/reference/resolve_bandwidth.md),
[`weighted_hsic_stat()`](https://max578.github.io/kernR/reference/weighted_hsic_stat.md)

## Examples

``` r
x <- matrix(rnorm(200), 100, 2)
select_bandwidth(x, "median")
#> [1] 1.749091
```
