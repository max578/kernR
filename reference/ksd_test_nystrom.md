# Accelerated Kernel Stein Discrepancy Goodness-of-Fit Test (Nystrom)

Low-rank counterpart to
[`ksd_test()`](https://max578.github.io/kernR/reference/ksd_test.md) for
large samples. The `n x n` Stein-kernel matrix is never materialised: it
is replaced by a Nystrom factor `F` (`n x m`, `m << n`) with \\F F^\top
\approx U\\, and both the KSD U-statistic and its wild-bootstrap null
are computed from `F` in `O(n m)` rather than `O(n^2)`. Because the
Langevin Stein kernel is itself positive semi-definite, \\F F^\top\\ is
a valid Stein kernel in its own right, so this is exactly the
[`ksd_test()`](https://max578.github.io/kernR/reference/ksd_test.md)
procedure applied to the rank-`m` kernel \\F F^\top\\: the
wild-bootstrap calibration of the degenerate U-statistic null is
preserved, and the approximation trades statistical power – not test
validity – for speed.

## Usage

``` r
ksd_test_nystrom(
  x,
  score = NULL,
  kernel = c("imq", "rbf"),
  beta = -0.5,
  bandwidth = "median",
  method = c("nystrom"),
  m = 100L,
  n_boot = 1000L,
  alpha = 0.05,
  seed = NULL,
  regularise = 1e-06
)
```

## Arguments

- x:

  Numeric vector, matrix, or data.frame. The `n x d` sample to test, one
  observation per row. At least five rows are required.

- score:

  Either `NULL` or a function. When `NULL` (default) the target is the
  standard multivariate normal, with score \\-x\\. When a function, it
  must accept the `n x d` sample matrix and return the `n x d` matrix of
  scores \\\nabla_x \log p(x)\\ evaluated row-wise; see
  [`gaussian_score()`](https://max578.github.io/kernR/reference/gaussian_score.md)
  for the multivariate-normal factory.

- kernel:

  Character. Base kernel for the Stein kernel: `"imq"` (inverse
  multi-quadric, default) or `"rbf"` (Gaussian).

- beta:

  Numeric in `(-1, 0)`. Exponent of the IMQ kernel; ignored when
  `kernel = "rbf"`. Default `-0.5`.

- bandwidth:

  Numeric or `"median"`. The IMQ offset `c` or the RBF bandwidth `h`.
  `"median"` (default) uses the median-heuristic bandwidth of the
  sample.

- method:

  Character. Factorisation method; currently `"nystrom"` only.

- m:

  Integer. Number of Nystrom landmarks (the approximation rank); capped
  at `n - 1`. Larger `m` improves power at higher cost. Default `100L`.

- n_boot:

  Integer. Number of wild-bootstrap replicates for the null. Default
  `1000`.

- alpha:

  Numeric in `(0, 1)`. Significance level for the verdict. Default
  `0.05`.

- seed:

  Integer or `NULL`. Random seed for reproducibility.

- regularise:

  Small positive numeric. Ridge added to the landmark Stein block before
  its Cholesky, for numerical stability. Default `1e-6`.

## Value

An object of class `c("ksd_test", "kernel_test_result")` carrying the
same fields as
[`ksd_test()`](https://max578.github.io/kernR/reference/ksd_test.md)
plus:

- approximation:

  `"nystrom"`.

- m:

  Effective rank used for the factorisation.

## Details

Currently Nystrom-only. Random-Fourier-feature factorisation of the
Stein kernel requires the analytic Fourier derivatives of the base
kernel and is deferred; the `method` argument is reserved for that
extension.

Use [`ksd_test()`](https://max578.github.io/kernR/reference/ksd_test.md)
for exact results at moderate `n`; reach for this when the `O(n^2)`
Stein matrix is the bottleneck. The verdict object and its fields match
[`ksd_test()`](https://max578.github.io/kernR/reference/ksd_test.md)
plus `approximation` and `m`.

## References

Liu, Q., Lee, J. D., & Jordan, M. I. (2016). A kernelized Stein
discrepancy for goodness-of-fit tests. *ICML*, PMLR 48, 276-284.

Chwialkowski, K., Strathmann, H., & Gretton, A. (2016). A kernel test of
goodness of fit. *ICML*, PMLR 48, 2606-2615.

Williams, C. K. I., & Seeger, M. (2001). Using the Nystrom method to
speed up kernel machines. *NeurIPS*, 13.

## See also

[`ksd_test()`](https://max578.github.io/kernR/reference/ksd_test.md),
[`nystrom_factor()`](https://max578.github.io/kernR/reference/nystrom_factor.md),
[`hsic_test_nystrom()`](https://max578.github.io/kernR/reference/hsic_test_nystrom.md),
[`concordance_test_nystrom()`](https://max578.github.io/kernR/reference/concordance_test_nystrom.md)

Other goodness-of-fit tests:
[`concordance_test()`](https://max578.github.io/kernR/reference/concordance_test.md),
[`concordance_test_nystrom()`](https://max578.github.io/kernR/reference/concordance_test_nystrom.md),
[`coverage_test()`](https://max578.github.io/kernR/reference/coverage_test.md),
[`gaussian_score()`](https://max578.github.io/kernR/reference/gaussian_score.md),
[`ksd_test()`](https://max578.github.io/kernR/reference/ksd_test.md),
[`numeric_score()`](https://max578.github.io/kernR/reference/numeric_score.md)

Other low-rank acceleration:
[`concordance_test_nystrom()`](https://max578.github.io/kernR/reference/concordance_test_nystrom.md),
[`hsic_test_nystrom()`](https://max578.github.io/kernR/reference/hsic_test_nystrom.md),
[`nystrom_factor()`](https://max578.github.io/kernR/reference/nystrom_factor.md),
[`rff_features()`](https://max578.github.io/kernR/reference/rff_features.md)

## Author

Max Moldovan, <max.moldovan@adelaide.edu.au>

## Examples

``` r
# \donttest{
set.seed(1)
x_ok  <- matrix(stats::rnorm(4000L), ncol = 2L)
fit_ok <- ksd_test_nystrom(x_ok, m = 80L, n_boot = 199L, seed = 1L)
fit_ok
#> 
#>    KSD (nystrom) Test
#> 
#> Statistic: 0.000282506 
#> P-value:   0.2050 
#> N:         2000 
#> Perms:     199 
#> Kernel X:  imq
#> 
#> Goodness-of-fit verdict
#>   Stein kernel: imq (beta = -0.5)
#>   Bandwidth:    1.726 (median heuristic)
#>   Bootstrap:    wild, B = 199
#>   Surprise:     2.286 bits
#>   Verdict:      consistent with target
#> 

x_bad <- x_ok + 1
ksd_test_nystrom(x_bad, m = 80L, n_boot = 199L, seed = 1L)
#> 
#>    KSD (nystrom) Test
#> 
#> Statistic: 0.82391 
#> P-value:   0.0050 
#> N:         2000 
#> Perms:     199 
#> Kernel X:  imq
#> 
#> Goodness-of-fit verdict
#>   Stein kernel: imq (beta = -0.5)
#>   Bandwidth:    1.726 (median heuristic)
#>   Bootstrap:    wild, B = 199
#>   Surprise:     7.644 bits
#>   Verdict:      REJECT (sample inconsistent with target)
#> 
# }
```
