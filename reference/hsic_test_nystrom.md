# HSIC Independence Test via Low-Rank Factorisation

Accelerated HSIC test using either Nystrom (default) or random Fourier
features (RFF) factorisations of the input kernel matrices. For large
`n` this scales as `O(n m)` per permutation instead of `O(n^2)`, with
`m << n` controlling the speed / accuracy trade-off.

## Usage

``` r
hsic_test_nystrom(
  x,
  y,
  kernel_x = kernel_spec(),
  kernel_y = kernel_spec(),
  method = c("nystrom", "rff"),
  m = 100L,
  m_x = NULL,
  m_y = NULL,
  n_permutations = 500L,
  alpha = 0.05,
  seed = NULL,
  regularise = 1e-06
)
```

## Arguments

- x, y:

  Numeric matrices (or vectors). Same number of rows.

- kernel_x, kernel_y:

  [`kernel_spec()`](https://max578.github.io/kernR/reference/kernel_spec.md)s
  for the two factors.

- method:

  Character. `"nystrom"` (default) or `"rff"`.

- m:

  Integer. Rank used for the approximation (number of Nystrom landmarks
  or RFF features). Used for both factors unless `m_x` / `m_y` are
  supplied.

- m_x, m_y:

  Optional integers overriding `m` per factor.

- n_permutations:

  Integer. Default `500L`.

- alpha:

  Numeric in `(0, 1)`. Default `0.05`.

- seed:

  Integer or `NULL`.

- regularise:

  Ridge for Nystrom Cholesky (ignored under `method = "rff"`). Default
  `1e-6`.

## Value

An object of class `"kernel_test_result"` with the standard fields
(`statistic`, `p_value`, `method`, `n`, `n_permutations`,
`null_distribution`, `kernel_x`, `kernel_y`, `call`) plus:

- approximation:

  `"nystrom"` or `"rff"`.

- m_x, m_y:

  Effective ranks used.

## Details

The test uses the *biased* HSIC estimator \\(1/n^2) \mathrm{tr}(H K_x H
K_y)\\, which is the form that factorises cleanly through low-rank
approximations. The bias is `O(1/n)` and negligible in the large-`n`
regime where Nystrom / RFF are useful.

The permutation null is built by row-permuting the (centred) `y` factor;
per-permutation cost is `O(n m_x m_y)`.

## References

Williams, C. K. I., & Seeger, M. (2001). Using the Nystrom method to
speed up kernel machines. *NeurIPS*, 13.

Rahimi, A., & Recht, B. (2007). Random features for large-scale kernel
machines. *NeurIPS*, 20.

Gretton, A., Bousquet, O., Smola, A., & Scholkopf, B. (2005). Measuring
statistical dependence with Hilbert-Schmidt norms. *ALT*, 63-77.

## See also

[`hsic_test()`](https://max578.github.io/kernR/reference/hsic_test.md),
[`nystrom_factor()`](https://max578.github.io/kernR/reference/nystrom_factor.md),
[`rff_features()`](https://max578.github.io/kernR/reference/rff_features.md)

Other low-rank acceleration:
[`concordance_test_nystrom()`](https://max578.github.io/kernR/reference/concordance_test_nystrom.md),
[`ksd_test_nystrom()`](https://max578.github.io/kernR/reference/ksd_test_nystrom.md),
[`nystrom_factor()`](https://max578.github.io/kernR/reference/nystrom_factor.md),
[`rff_features()`](https://max578.github.io/kernR/reference/rff_features.md)

## Examples

``` r
# \donttest{
set.seed(1)
n <- 1500L
x <- stats::rnorm(n)
y <- x^2 + stats::rnorm(n, sd = 0.5)
fit <- hsic_test_nystrom(x, y, m = 80L,
                         n_permutations = 99L, seed = 1L)
fit
#> 
#>    HSIC (nystrom) Test
#> 
#> Statistic: 0.0211099 
#> P-value:   0.0100 
#> N:         1500 
#> Perms:     99 
#> Kernel X:  rbf (bw = 0.9793)
#> Kernel Y:  rbf (bw =  1.02)
#> 
# }
```
