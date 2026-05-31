# Nystrom Low-Rank Kernel Factorisation

Builds a Nystrom approximation \\K \approx F F^\top\\ of an `n x n`
kernel matrix using `m << n` uniformly-sampled landmarks. The returned
factor `F` is an `n x m` matrix; downstream kernel computations (HSIC,
MMD, etc.) that can be expressed in `F` reduce from `O(n^2)` to `O(n m)`
cost.

## Usage

``` r
nystrom_factor(
  x,
  kernel = kernel_spec(),
  m = 100L,
  regularise = 1e-06,
  seed = NULL
)
```

## Arguments

- x:

  Numeric matrix `n x d` (or vector). Data points.

- kernel:

  A
  [`kernel_spec()`](https://max578.github.io/kernR/reference/kernel_spec.md).
  Default RBF with median heuristic.

- m:

  Integer. Number of landmark points. Capped at `nrow(x) - 1`. Default
  `100L`.

- regularise:

  Small positive numeric. Ridge added to `W` before Cholesky for
  numerical stability. Default `1e-6`.

- seed:

  Integer or `NULL`. Random seed for landmark sampling.

## Value

An object of class `"kernel_factor"` with components:

- F:

  The `n x m_eff` factor matrix.

- method:

  `"nystrom"`.

- m:

  Effective rank (`<= m`).

- kernel:

  Resolved kernel spec.

- n:

  Number of rows in the input.

## Details

The construction is:

1.  Sample `m` landmark indices uniformly without replacement.

2.  Compute the landmark Gram \\W = K\_{mm}\\ (`m x m`) and the
    cross-Gram \\C = K\_{nm}\\ (`n x m`).

3.  Stabilise: \\W\_\epsilon = W + \epsilon I\\.

4.  Cholesky factor \\W\_\epsilon = L L^\top\\.

5.  Return \\F = C L^{-\top}\\, so that \\F F^\top = C W\_\epsilon^{-1}
    C^\top \approx K\\.

Any
[`kernel_spec()`](https://max578.github.io/kernR/reference/kernel_spec.md)
is supported. Bandwidth selection (median heuristic) is performed on the
full dataset before landmark sampling.

## References

Williams, C. K. I., & Seeger, M. (2001). Using the Nystrom method to
speed up kernel machines. *NeurIPS*, 13.

## See also

[`rff_features()`](https://max578.github.io/kernR/reference/rff_features.md),
[`hsic_test_nystrom()`](https://max578.github.io/kernR/reference/hsic_test_nystrom.md)

Other low-rank acceleration:
[`hsic_test_nystrom()`](https://max578.github.io/kernR/reference/hsic_test_nystrom.md),
[`rff_features()`](https://max578.github.io/kernR/reference/rff_features.md)

## Examples

``` r
set.seed(1)
x <- matrix(stats::rnorm(2000L), ncol = 2L)
f <- nystrom_factor(x, m = 80L, seed = 1L)
dim(f$F)
#> [1] 1000   80
```
