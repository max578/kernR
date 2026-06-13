# Scaling HSIC with Nystrom and Random Fourier Features

## The scaling wall

Exact kernel methods carry an `O(n^2)` storage cost (the Gram matrix)
and `O(n^2)` or `O(n^3)` compute per operation. For ag-systems ensembles
in the low thousands of members, that wall hits around `n = 5000`: a
5000 x 5000 double-precision Gram matrix is 200 MB, and permutation HSIC
iterates over it `B` times.

Two well-known low-rank approximations of a kernel matrix `K` (`n x n`)
drop both costs to `O(n m)` for some rank `m << n`:

- **Nystrom** (Williams & Seeger, 2001) – sample `m` landmark
  observations, compute their Gram `W` (`m x m`) and cross-Gram `C`
  (`n x m`), then `K \approx C W^{-1} C^\top = F F^\top` with
  `F = C \cdot \mathrm{chol}(W)^{-1}`. Works for any kernel.
- **Random Fourier Features** (Rahimi & Recht, 2007) – for
  shift-invariant kernels (currently RBF in kernR), draw
  `\omega_k \sim N(0, \sigma^{-2} I)` and `b_k \sim U[0, 2\pi]`; the
  feature map `\phi(x) = \sqrt{2/D} \cos(\omega^\top x + b)` gives
  `K \approx \Phi \Phi^\top` for `D` features. Data-independent (the
  same random projection applies to any point cloud with the same
  bandwidth).

[`hsic_test_nystrom()`](https://max578.github.io/kernR/reference/hsic_test_nystrom.md)
wires both into a drop-in accelerated HSIC test.

## A tiny correctness check

Both approximations recover the exact kernel as the rank approaches `n`
(Nystrom) or as `D` grows (RFF).

``` r

library(kernR)

n <- 80L
x <- matrix(stats::rnorm(n * 2L), n, 2L)
k <- kernel_spec("rbf", bandwidth = 1.0)
K_full <- kernel_matrix(x, kernel = k)

# Nystrom at near-full rank
ny <- nystrom_factor(x, kernel = k, m = n - 1L, seed = 1L)
K_ny <- tcrossprod(ny$F)

# RFF at a moderate D
rf <- rff_features(x, kernel = k, D = 1500L, seed = 1L)
K_rf <- tcrossprod(rf$F)

rel_err <- function(K_hat) {
  sqrt(sum((K_full - K_hat)^2)) / sqrt(sum(K_full^2))
}
data.frame(
  method = c("nystrom", "rff"),
  rank   = c(ncol(ny$F), ncol(rf$F)),
  rel_err = c(rel_err(K_ny), rel_err(K_rf))
)
#>    method rank      rel_err
#> 1 nystrom   79 1.849946e-07
#> 2     rff 1500 4.728235e-02
```

Both should produce small relative Frobenius error on this size.

## Scaling demonstration

``` r

benchmark_hsic <- function(n, m, B = 49L, seed = 1L) {
  set.seed(seed)
  xx <- stats::rnorm(n)
  y  <- xx + stats::rnorm(n, sd = 0.5)

  t_exact <- system.time(
    res_exact <- hsic_test(xx, y, n_permutations = B, seed = seed)
  )["elapsed"]
  t_nystrom <- system.time(
    res_nystrom <- hsic_test_nystrom(xx, y, m = m,
                                     n_permutations = B, seed = seed)
  )["elapsed"]
  t_rff <- system.time(
    res_rff <- hsic_test_nystrom(xx, y, method = "rff", m = m,
                                 n_permutations = B, seed = seed)
  )["elapsed"]
  data.frame(
    n         = n,
    m         = m,
    method    = c("hsic_test (exact)", "nystrom", "rff"),
    elapsed_s = round(c(t_exact, t_nystrom, t_rff), 3),
    p_value   = round(c(res_exact$p_value, res_nystrom$p_value,
                        res_rff$p_value), 3)
  )
}

benchmark_hsic(n = 500L,  m = 60L)
#>     n  m            method elapsed_s p_value
#> 1 500 60 hsic_test (exact)     0.096    0.02
#> 2 500 60           nystrom     0.020    0.02
#> 3 500 60               rff     0.017    0.02
```

For larger `n`, the gap widens:

``` r

benchmark_hsic(n = 1500L, m = 100L)
#>      n   m            method elapsed_s p_value
#> 1 1500 100 hsic_test (exact)     1.217    0.02
#> 2 1500 100           nystrom     0.133    0.02
#> 3 1500 100               rff     0.131    0.02
```

The verdict (reject vs accept) agrees across exact and approximate tests
at moderate `m`; the approximation is for speed, not for detecting a
different signal.

## When to use which

| Scenario | Pick |
|----|----|
| Small `n` (\< 500) | [`hsic_test()`](https://max578.github.io/kernR/reference/hsic_test.md) – exact, fast enough. |
| Medium `n` (500-5000), any kernel | `hsic_test_nystrom(method = "nystrom")`. Choose `m` in `[50, 200]`. |
| Large `n` (\> 5000), RBF kernel | `hsic_test_nystrom(method = "rff")` with `D` in `[200, 1000]`. Data-independent projection; cheapest at large n. |
| Repeated tests on different `y` with same `x` | Cache `nystrom_factor(x)` or `rff_features(x)` once; pass the factor in (future API). |

## Notes on practice

- **Bias vs unbiased estimator.** The exact
  [`hsic_test()`](https://max578.github.io/kernR/reference/hsic_test.md)
  uses the Gretton-2008 *unbiased* HSIC. The Nystrom/RFF version uses
  the *biased* estimator `(1/n^2) tr(H K_x H K_y)`, because that form
  factors cleanly through low-rank approximations. The bias is `O(1/n)`
  – negligible in the large-`n` regime where these approximations are
  useful.
- **Choosing `m` for Nystrom.** Rule of thumb: `m = ceiling(2 sqrt(n))`
  is a defensible starting point. Diagnose via the Frobenius
  reconstruction error on a small subsample if precision matters.
- **Choosing `D` for RFF.** RFF variance scales as `O(1/D)`; the bias is
  zero in expectation. Larger `D` is always better in accuracy but costs
  memory; `D = 200-500` is usually adequate for verdict-style tests.
- **Reproducibility.** Both factorisations are randomised; pass `seed`
  for deterministic output.

## References

- Williams, C. K. I., & Seeger, M. (2001). Using the Nystrom method to
  speed up kernel machines. *NeurIPS*, 13.
- Rahimi, A., & Recht, B. (2007). Random features for large-scale kernel
  machines. *NeurIPS*, 20.
- Gretton, A., Bousquet, O., Smola, A., & Schölkopf, B. (2005).
  Measuring statistical dependence with Hilbert-Schmidt norms. *ALT*,
  63-77.
