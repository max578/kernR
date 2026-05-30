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

## Examples

``` r
x <- matrix(rnorm(200), 100, 2)
select_bandwidth(x, "median")
#> [1] 1.749091
```
