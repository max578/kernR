# Kernel k-sample Concordance Test

Tests whether two or more samples come from a common distribution, using
the summed pairwise Maximum Mean Discrepancy with a joint-permutation
null. The samples are typically posterior draws from different inference
engines, or scenario ensembles from different simulators; the test asks
whether they are mutually concordant. Unlike repeated two-sample
[`mmd_test()`](https://max578.github.io/kernR/reference/mmd_test.md)
calls, the null is a single shared relabeling of the pooled sample, so
the family-wise error is controlled and the overall verdict is one
calibrated p-value.

## Usage

``` r
concordance_test(
  x,
  kernel = kernel_spec(),
  n_permutations = 500L,
  alpha = 0.05,
  seed = NULL,
  n_exact_max = 5000L
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

- n_permutations:

  Integer. Number of joint permutations for the null. Default `500`.

- alpha:

  Numeric in `(0, 1)`. Significance level for the verdict. Default
  `0.05`.

- seed:

  Integer or `NULL`. Random seed for reproducibility.

- n_exact_max:

  Integer or `Inf`. Pooled-sample-size ceiling for the exact `O(n^2)`
  test. Above it, the call is delegated to
  [`concordance_test_nystrom()`](https://max578.github.io/kernR/reference/concordance_test_nystrom.md)
  (with a message; the verdict object records
  `approximation = "nystrom"`). `Inf` forces the exact test at any size.
  Default `5000L`.

## Value

An object of class `c("concordance_test", "kernel_test_result")`
carrying the standard `kernel_test_result` fields plus:

- statistic:

  Summed pairwise unbiased MMD-squared.

- p_value:

  Upper-tail joint-permutation p-value (with `+1` correction).

- n_groups:

  Number of sources compared.

- group_sizes:

  Named integer vector of per-source sample sizes.

- pairwise:

  Symmetric `K x K` matrix of pairwise unbiased MMD-squared,
  row/column-named by source.

- alpha, reject:

  Verdict level and `p_value <= alpha`.

## Details

The returned object carries the full pairwise MMD discrepancy matrix, so
a rejection can be read down to the offending pair: convergence across
sources is corroborating evidence, and divergence localises which source
departs and on which margin. This is the cross-engine concordance role –
a sample-based complement to the score-based
[`ksd_test()`](https://max578.github.io/kernR/reference/ksd_test.md).

The exact test materialises the pooled `n x n` kernel matrix (`O(n^2)`).
To keep large ensembles tractable without a silent loss of exactness, a
pooled sample with more than `n_exact_max` rows is delegated to
[`concordance_test_nystrom()`](https://max578.github.io/kernR/reference/concordance_test_nystrom.md)
– a low-rank approximation that is *announced* by a message and
*recorded* in the returned object's `approximation` and `m` fields. Set
`n_exact_max = Inf` to force the exact test, or call
[`concordance_test_nystrom()`](https://max578.github.io/kernR/reference/concordance_test_nystrom.md)
directly to control the approximation rank `m`.

## References

Gretton, A., Borgwardt, K. M., Rasch, M. J., Scholkopf, B., & Smola, A.
(2012). A kernel two-sample test. *Journal of Machine Learning
Research*, 13, 723-773.

## See also

[`mmd_test()`](https://max578.github.io/kernR/reference/mmd_test.md),
[`ksd_test()`](https://max578.github.io/kernR/reference/ksd_test.md),
[`mmd_ppc()`](https://max578.github.io/kernR/reference/mmd_ppc.md)

Other goodness-of-fit tests:
[`concordance_test_nystrom()`](https://max578.github.io/kernR/reference/concordance_test_nystrom.md),
[`coverage_test()`](https://max578.github.io/kernR/reference/coverage_test.md),
[`gaussian_score()`](https://max578.github.io/kernR/reference/gaussian_score.md),
[`joint_coverage_test()`](https://max578.github.io/kernR/reference/joint_coverage_test.md),
[`ksd_test()`](https://max578.github.io/kernR/reference/ksd_test.md),
[`ksd_test_nystrom()`](https://max578.github.io/kernR/reference/ksd_test_nystrom.md),
[`numeric_score()`](https://max578.github.io/kernR/reference/numeric_score.md)

## Author

Max Moldovan, <max.moldovan@adelaide.edu.au>

## Examples

``` r
set.seed(1)

# Three concordant sources (same distribution): not rejected
draws <- list(
  engine_a = matrix(stats::rnorm(400L), ncol = 2L),
  engine_b = matrix(stats::rnorm(400L), ncol = 2L),
  engine_c = matrix(stats::rnorm(400L), ncol = 2L)
)
fit_ok <- concordance_test(draws, n_permutations = 199L, seed = 1L)
fit_ok
#> 
#>    Concordance Test
#> 
#> Statistic: 0.000908946 
#> P-value:   0.3800 
#> N:         600 
#> Perms:     199 
#> Kernel X:  rbf (bw = 1.715)
#> 
#> Concordance verdict
#>   Sources:    3 (engine_a, engine_b, engine_c)
#>   Verdict:    concordant
#>   Pairwise MMD-squared:
#>          engine_a  engine_b engine_c 
#> engine_a    0      0.00341  -0.000463
#> engine_b 0.00341      0     -0.00204 
#> engine_c -0.000463 -0.00204    0     
#> 

# One source departs (mean-shifted): rejected, and the pairwise matrix
# localises engine_c
draws$engine_c <- draws$engine_c + 1
fit_bad <- concordance_test(draws, n_permutations = 199L, seed = 1L)
fit_bad$pairwise
#>             engine_a    engine_b  engine_c
#> engine_a 0.000000000 0.002898312 0.1882500
#> engine_b 0.002898312 0.000000000 0.2157871
#> engine_c 0.188249967 0.215787070 0.0000000
```
