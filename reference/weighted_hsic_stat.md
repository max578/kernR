# Weighted HSIC Statistic

Computes the weighted Hilbert–Schmidt Independence Criterion (HSIC)
between two pre-computed kernel matrices, \$\$\frac{1}{W^2} \sum\_{i,j}
w_i w_j (K_x^c)\_{ij} (K_y^c)\_{ij},\$\$ where `K_x^c` and `K_y^c` are
weight-centred and \\W = \sum_i w_i\\. This is the exact statistic
[`bd_hsic_test()`](https://max578.github.io/kernR/reference/bd_hsic_test.md)
accumulates over its permutation null, exposed so that downstream
methods – for example theory-anchored causal inference – can build a
bd-HSIC-compatible statistic against a custom reference distribution
without re-implementing the engine.

## Usage

``` r
weighted_hsic_stat(Kx, Ky, w = rep(1, nrow(Kx)))
```

## Arguments

- Kx:

  Numeric `n` by `n` kernel matrix for the first variable.

- Ky:

  Numeric `n` by `n` kernel matrix for the second variable, with the
  same dimensions as `Kx`.

- w:

  Non-negative numeric weight vector of length `n`. Defaults to uniform
  weights, giving the unweighted HSIC.

## Value

A single numeric value: the weighted HSIC statistic.

## Details

With uniform weights this reduces to the (biased) HSIC estimator;
supplying density-ratio weights \\w(t, z)\\ yields the backdoor-adjusted
statistic used for causal testing.

## See also

[`bd_hsic_test()`](https://max578.github.io/kernR/reference/bd_hsic_test.md)
for the full causal test that uses this statistic;
[`resolve_bandwidth()`](https://max578.github.io/kernR/reference/resolve_bandwidth.md)
and
[`kernel_matrix()`](https://max578.github.io/kernR/reference/kernel_matrix.md)
for building the kernel inputs.

Other kernel primitives:
[`kernel_matrix()`](https://max578.github.io/kernR/reference/kernel_matrix.md),
[`kernel_spec()`](https://max578.github.io/kernR/reference/kernel_spec.md),
[`resolve_bandwidth()`](https://max578.github.io/kernR/reference/resolve_bandwidth.md),
[`select_bandwidth()`](https://max578.github.io/kernR/reference/select_bandwidth.md)

## Examples

``` r
set.seed(1)
x <- matrix(rnorm(40), ncol = 2)
y <- matrix(rnorm(40), ncol = 2)
Kx <- kernel_matrix(x, kernel = resolve_bandwidth(kernel_spec(), x))
Ky <- kernel_matrix(y, kernel = resolve_bandwidth(kernel_spec(), y))
weighted_hsic_stat(Kx, Ky)
#> [1] 0.007846936
```
