# Accelerated Kernel k-sample Concordance Test (Nystrom / RFF)

Low-rank counterpart to
[`concordance_test()`](https://max578.github.io/kernR/reference/concordance_test.md)
for large ensembles. The pooled sample is factorised once – by the
Nystrom method (default) or random Fourier features – into an `n x m`
factor `F` with \\F F^\top \approx K\\; the summed pairwise unbiased
MMD-squared and its joint-permutation null are then computed from `F` in
`O(n m)` per permutation rather than `O(n^2)`, with `m << n` controlling
the speed / accuracy trade-off. The verdict object and its
interpretation are identical to
[`concordance_test()`](https://max578.github.io/kernR/reference/concordance_test.md);
only the cost scales differently.

## Usage

``` r
concordance_test_nystrom(
  x,
  kernel = kernel_spec(),
  method = c("nystrom", "rff"),
  m = 100L,
  n_permutations = 500L,
  alpha = 0.05,
  seed = NULL,
  regularise = 1e-06
)
```

## Arguments

- x:

  A list of two or more samples to compare. Each element is a numeric
  vector, matrix, or data.frame with `n_k` observations (rows) over a
  shared `d` columns. A named list labels the sources in the output; an
  unnamed list is labelled `Source 1`, `Source 2`, and so on. Each
  sample needs at least five rows.

- kernel:

  Kernel specification. Default is RBF with the median heuristic over
  the pooled sample.

- method:

  Character. `"nystrom"` (default) or `"rff"`. RFF requires an RBF
  `kernel`.

- m:

  Integer. Rank of the approximation: the number of Nystrom landmarks or
  RFF features. Larger `m` improves accuracy at higher cost. Default
  `100L`.

- n_permutations:

  Integer. Number of joint permutations for the null. Default `500`.

- alpha:

  Numeric in `(0, 1)`. Significance level for the verdict. Default
  `0.05`.

- seed:

  Integer or `NULL`. Random seed for reproducibility.

- regularise:

  Small positive numeric. Ridge added before the Nystrom Cholesky for
  numerical stability; ignored under `method = "rff"`. Default `1e-6`.

## Value

An object of class `c("concordance_test", "kernel_test_result")`
carrying the same fields as
[`concordance_test()`](https://max578.github.io/kernR/reference/concordance_test.md)
plus:

- approximation:

  `"nystrom"` or `"rff"`.

- m:

  Effective rank used for the factorisation.

## Details

The factorisation of the pooled sample preserves the per-source mean
embeddings (per-source column sums of `F`), so the pairwise discrepancy
matrix still localises which source departs. The joint-permutation null
is built by relabelling the rows of `F` – the low-rank analogue of
permuting the pooled-sample labels in
[`concordance_test()`](https://max578.github.io/kernR/reference/concordance_test.md).

Use
[`concordance_test()`](https://max578.github.io/kernR/reference/concordance_test.md)
for exact results at moderate `n`; reach for this function when the
pooled sample is large enough that the `O(n^2)` kernel matrix is the
bottleneck. RFF (`method = "rff"`) requires an RBF kernel; Nystrom
supports any
[`kernel_spec()`](https://max578.github.io/kernR/reference/kernel_spec.md).

## References

Gretton, A., Borgwardt, K. M., Rasch, M. J., Scholkopf, B., & Smola, A.
(2012). A kernel two-sample test. *Journal of Machine Learning
Research*, 13, 723-773.

Williams, C. K. I., & Seeger, M. (2001). Using the Nystrom method to
speed up kernel machines. *NeurIPS*, 13.

Rahimi, A., & Recht, B. (2007). Random features for large-scale kernel
machines. *NeurIPS*, 20.

## See also

[`concordance_test()`](https://max578.github.io/kernR/reference/concordance_test.md),
[`nystrom_factor()`](https://max578.github.io/kernR/reference/nystrom_factor.md),
[`rff_features()`](https://max578.github.io/kernR/reference/rff_features.md),
[`hsic_test_nystrom()`](https://max578.github.io/kernR/reference/hsic_test_nystrom.md)

Other goodness-of-fit tests:
[`concordance_test()`](https://max578.github.io/kernR/reference/concordance_test.md),
[`coverage_test()`](https://max578.github.io/kernR/reference/coverage_test.md),
[`gaussian_score()`](https://max578.github.io/kernR/reference/gaussian_score.md),
[`ksd_test()`](https://max578.github.io/kernR/reference/ksd_test.md),
[`ksd_test_nystrom()`](https://max578.github.io/kernR/reference/ksd_test_nystrom.md),
[`numeric_score()`](https://max578.github.io/kernR/reference/numeric_score.md)

Other low-rank acceleration:
[`hsic_test_nystrom()`](https://max578.github.io/kernR/reference/hsic_test_nystrom.md),
[`ksd_test_nystrom()`](https://max578.github.io/kernR/reference/ksd_test_nystrom.md),
[`nystrom_factor()`](https://max578.github.io/kernR/reference/nystrom_factor.md),
[`rff_features()`](https://max578.github.io/kernR/reference/rff_features.md)

## Author

Max Moldovan, <max.moldovan@adelaide.edu.au>

## Examples

``` r
# \donttest{
set.seed(1)
big <- list(
  engine_a = matrix(stats::rnorm(4000L), ncol = 2L),
  engine_b = matrix(stats::rnorm(4000L), ncol = 2L),
  engine_c = matrix(stats::rnorm(4000L), ncol = 2L) + 0.4
)
fit <- concordance_test_nystrom(big, m = 80L,
                                n_permutations = 199L, seed = 1L)
fit
#> 
#>    Concordance (nystrom) Test
#> 
#> Statistic: 0.0747907 
#> P-value:   0.0050 
#> N:         6000 
#> Perms:     199 
#> Kernel X:  rbf (bw = 1.683)
#> 
#> Concordance verdict
#>   Sources:    3 (engine_a, engine_b, engine_c)
#>   Verdict:    REJECT (sources are not mutually concordant)
#>   Pairwise MMD-squared:
#>          engine_a engine_b engine_c
#> engine_a    0     6.36e-05 0.0348  
#> engine_b 6.36e-05    0     0.0399  
#> engine_c 0.0348   0.0399      0    
#> 
fit$pairwise
#>              engine_a     engine_b   engine_c
#> engine_a 0.000000e+00 6.362778e-05 0.03479769
#> engine_b 6.362778e-05 0.000000e+00 0.03992936
#> engine_c 3.479769e-02 3.992936e-02 0.00000000
# }
```
