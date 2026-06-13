# Kernel Stein Discrepancy Goodness-of-Fit Test

Tests whether a sample is consistent with a target distribution `p`,
where `p` is supplied through its score \\\nabla_x \log p(x)\\ rather
than a reference sample. The statistic is the (unbiased) kernel Stein
discrepancy (KSD); calibration uses a wild bootstrap of the degenerate
U-statistic null. Because only the score enters, the target may be known
up to an unknown normalising constant.

## Usage

``` r
ksd_test(
  x,
  score = NULL,
  kernel = c("imq", "rbf"),
  beta = -0.5,
  bandwidth = "median",
  n_boot = 1000L,
  alpha = 0.05,
  seed = NULL,
  n_exact_max = 5000L
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

- n_boot:

  Integer. Number of wild-bootstrap replicates for the null. Default
  `1000`.

- alpha:

  Numeric in `(0, 1)`. Significance level for the verdict. Default
  `0.05`.

- seed:

  Integer or `NULL`. Random seed for reproducibility.

- n_exact_max:

  Integer or `Inf`. Sample-size ceiling for the exact `O(n^2)` test.
  Above it, the call is delegated to
  [`ksd_test_nystrom()`](https://max578.github.io/kernR/reference/ksd_test_nystrom.md)
  (with a message; the verdict object records
  `approximation = "nystrom"`). `Inf` forces the exact test at any size.
  Default `5000L`.

## Value

An object of class `c("ksd_test", "kernel_test_result")` carrying the
standard `kernel_test_result` fields plus:

- statistic:

  The unbiased KSD U-statistic.

- p_value:

  Upper-tail wild-bootstrap p-value (with `+1` correction).

- stein_kernel:

  Base kernel used (`"imq"` or `"rbf"`).

- beta:

  IMQ exponent, or `NA` for the RBF kernel.

- bandwidth:

  Resolved IMQ offset or RBF bandwidth.

- surprise_bits:

  Shannon-information surprise `-log2(p_value)`.

- alpha, reject:

  Verdict level and `p_value <= alpha`.

## Details

KSD is the score-based, one-sample complement to
[`mmd_test()`](https://max578.github.io/kernR/reference/mmd_test.md):
where MMD compares two samples, KSD compares a sample against a
*density*. The calibration framing is direct – given posterior or
ensemble draws and the score of the distribution they claim to
represent, KSD asks whether the draws actually follow that distribution.
It is sensitive to mean, variance, and tail mis-specification.

The default base kernel is the inverse multi-quadric (IMQ), \\k(x, y) =
(c^2 + \lVert x - y \rVert^2)^\beta\\ with \\\beta \in (-1, 0)\\. Gorham
& Mackey (2017) show the IMQ Stein discrepancy detects non-convergence
in regimes where the Gaussian (RBF) Stein discrepancy is blind,
particularly as dimension grows; the RBF base kernel remains available
via `kernel = "rbf"`. The offset `c` (IMQ) and bandwidth `h` (RBF)
default to the median heuristic over the sample.

Reproducibility: the wild-bootstrap multipliers are drawn through R's
RNG, so a non-`NULL` `seed` makes the p-value reproducible under the
active RNG kind (the R default Mersenne-Twister unless changed by the
caller).

The exact test materialises the `n x n` Stein-kernel matrix, so memory
and compute scale as `O(n^2)`. To keep large samples tractable without a
silent loss of exactness, a sample with more than `n_exact_max` rows is
delegated to
[`ksd_test_nystrom()`](https://max578.github.io/kernR/reference/ksd_test_nystrom.md)
– a low-rank approximation that is *announced* by a message and
*recorded* in the returned object's `approximation` and `m` fields, so
the result stays reproducible from the object. Set `n_exact_max = Inf`
to force the exact `O(n^2)` test at any size, or call
[`ksd_test_nystrom()`](https://max578.github.io/kernR/reference/ksd_test_nystrom.md)
directly to control the approximation rank `m`.

## References

Liu, Q., Lee, J. D., & Jordan, M. I. (2016). A kernelized Stein
discrepancy for goodness-of-fit tests. *Proceedings of the 33rd
International Conference on Machine Learning*, PMLR 48, 276-284.

Chwialkowski, K., Strathmann, H., & Gretton, A. (2016). A kernel test of
goodness of fit. *Proceedings of the 33rd International Conference on
Machine Learning*, PMLR 48, 2606-2615.

Gorham, J., & Mackey, L. (2017). Measuring sample quality with kernels.
*Proceedings of the 34th International Conference on Machine Learning*,
PMLR 70, 1292-1301.

## See also

[`mmd_test()`](https://max578.github.io/kernR/reference/mmd_test.md),
[`mmd_ppc()`](https://max578.github.io/kernR/reference/mmd_ppc.md),
[`gaussian_score()`](https://max578.github.io/kernR/reference/gaussian_score.md)

Other goodness-of-fit tests:
[`concordance_test()`](https://max578.github.io/kernR/reference/concordance_test.md),
[`concordance_test_nystrom()`](https://max578.github.io/kernR/reference/concordance_test_nystrom.md),
[`coverage_test()`](https://max578.github.io/kernR/reference/coverage_test.md),
[`gaussian_score()`](https://max578.github.io/kernR/reference/gaussian_score.md),
[`joint_coverage_test()`](https://max578.github.io/kernR/reference/joint_coverage_test.md),
[`ksd_test_nystrom()`](https://max578.github.io/kernR/reference/ksd_test_nystrom.md),
[`numeric_score()`](https://max578.github.io/kernR/reference/numeric_score.md)

## Author

Max Moldovan, <max.moldovan@adelaide.edu.au>

## Examples

``` r
set.seed(1)

# Well-specified: standard-normal sample against standard-normal target
x_ok <- matrix(stats::rnorm(400L), ncol = 2L)
fit_ok <- ksd_test(x_ok, n_boot = 199L, seed = 1L)
fit_ok
#> 
#>    KSD Test
#> 
#> Statistic: -0.00225082 
#> P-value:   0.6300 
#> N:         200 
#> Perms:     199 
#> Kernel X:  imq
#> 
#> Goodness-of-fit verdict
#>   Stein kernel: imq (beta = -0.5)
#>   Bandwidth:    1.608 (median heuristic)
#>   Bootstrap:    wild, B = 199
#>   Surprise:     0.667 bits
#>   Verdict:      consistent with target
#> 

# Mis-specified: mean-shifted sample against the same target
x_bad <- x_ok + 1
fit_bad <- ksd_test(x_bad, n_boot = 199L, seed = 1L)
fit_bad
#> 
#>    KSD Test
#> 
#> Statistic: 0.946639 
#> P-value:   0.0050 
#> N:         200 
#> Perms:     199 
#> Kernel X:  imq
#> 
#> Goodness-of-fit verdict
#>   Stein kernel: imq (beta = -0.5)
#>   Bandwidth:    1.608 (median heuristic)
#>   Bootstrap:    wild, B = 199
#>   Surprise:     7.644 bits
#>   Verdict:      REJECT (sample inconsistent with target)
#> 

# Explicit non-standard target via the Gaussian score factory
sig <- matrix(c(1, 0.6, 0.6, 1), nrow = 2L)
x_cor <- x_ok %*% chol(sig)
ksd_test(x_cor, score = gaussian_score(sigma = sig),
         n_boot = 199L, seed = 1L)
#> 
#>    KSD Test
#> 
#> Statistic: -0.00438649 
#> P-value:   0.6050 
#> N:         200 
#> Perms:     199 
#> Kernel X:  imq
#> 
#> Goodness-of-fit verdict
#>   Stein kernel: imq (beta = -0.5)
#>   Bandwidth:    1.472 (median heuristic)
#>   Bootstrap:    wild, B = 199
#>   Surprise:     0.725 bits
#>   Verdict:      consistent with target
#> 
```
