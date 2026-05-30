# Random Fourier Features for the RBF Kernel

Builds an `n x D` feature map \\\Phi\\ such that \\\Phi \Phi^\top
\approx K\\ for the RBF kernel \\K(x, y) = \exp(-\\x-y\\^2 /
(2\sigma^2))\\, via Rahimi & Recht (2007): draw \\\omega_k \sim N(0,
\sigma^{-2} I_d)\\ and \\b_k \sim U\[0, 2\pi\]\\, then \\\phi_k(x) =
\sqrt{2/D} \cos(\omega_k^\top x + b_k)\\.

## Usage

``` r
rff_features(x, kernel = kernel_spec("rbf"), D = 200L, seed = NULL)
```

## Arguments

- x:

  Numeric matrix `n x d` (or vector).

- kernel:

  A
  [`kernel_spec()`](https://max578.github.io/kernR/reference/kernel_spec.md)
  with `type = "rbf"`. Bandwidth may be `"median"` (resolved against
  `x`) or a positive numeric.

- D:

  Integer. Number of random features. Larger D -\> better approximation
  but higher memory / compute. Default `200L`.

- seed:

  Integer or `NULL`. Random seed.

## Value

An object of class `"kernel_factor"` with components:

- F:

  The `n x D` random-feature matrix.

- method:

  `"rff"`.

- m:

  Equal to `D`.

- kernel:

  Resolved kernel spec.

- n:

  Number of rows in the input.

- omega, b:

  Random draws (kept for reproducible re-encoding).

## Details

Currently RBF-only; other shift-invariant kernels (Matern with specific
`nu`) require their own Fourier spectra and are not yet implemented.

## References

Rahimi, A., & Recht, B. (2007). Random features for large-scale kernel
machines. *NeurIPS*, 20.

## See also

[`nystrom_factor()`](https://max578.github.io/kernR/reference/nystrom_factor.md),
[`hsic_test_nystrom()`](https://max578.github.io/kernR/reference/hsic_test_nystrom.md)

## Examples

``` r
set.seed(1)
x <- matrix(stats::rnorm(2000L), ncol = 2L)
f <- rff_features(x, D = 150L, seed = 1L)
dim(f$F)
#> [1] 1000  150
```
