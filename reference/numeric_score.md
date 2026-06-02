# Finite-difference score from a log-density

Builds a score function – the gradient of the log density – by central
finite differences, for use as the `score` argument of
[`ksd_test()`](https://max578.github.io/kernR/reference/ksd_test.md).
The target need only be expressible as a function returning a (possibly
unnormalised) log density; any additive normalising constant cancels in
the gradient, so the constant may be omitted.

## Usage

``` r
numeric_score(log_density, h = 1e-04)
```

## Arguments

- log_density:

  A function accepting an `n x d` numeric matrix and returning a numeric
  vector of length `n`: the log density (up to an additive constant)
  evaluated row-wise.

- h:

  Positive numeric. Finite-difference step. Default `1e-4`. Too large
  biases the gradient; too small amplifies floating-point noise.

## Value

A function of one argument (an `n x d` numeric matrix) returning the
`n x d` matrix of finite-difference scores, suitable as the `score`
argument of
[`ksd_test()`](https://max578.github.io/kernR/reference/ksd_test.md).

## Details

This makes
[`ksd_test()`](https://max578.github.io/kernR/reference/ksd_test.md)
usable against any target a user can write down as a log density, not
only targets with a hand-derived score. It is also the kernR-side
adapter for a log-posterior contract: wrapping an exported log-posterior
evaluator in `numeric_score()` yields the score
[`ksd_test()`](https://max578.github.io/kernR/reference/ksd_test.md)
needs to check whether posterior draws are calibrated against that
posterior, with no dependency added on the producer – the evaluator is
passed in as an ordinary function.

Central differences cost `2 * d` evaluations of `log_density` per call,
where `d` is the dimension. For a cheap closed-form target prefer
[`gaussian_score()`](https://max578.github.io/kernR/reference/gaussian_score.md)
or a hand-written score; `numeric_score()` is the general fallback.

## See also

[`ksd_test()`](https://max578.github.io/kernR/reference/ksd_test.md),
[`gaussian_score()`](https://max578.github.io/kernR/reference/gaussian_score.md)

Other goodness-of-fit tests:
[`concordance_test()`](https://max578.github.io/kernR/reference/concordance_test.md),
[`coverage_test()`](https://max578.github.io/kernR/reference/coverage_test.md),
[`gaussian_score()`](https://max578.github.io/kernR/reference/gaussian_score.md),
[`ksd_test()`](https://max578.github.io/kernR/reference/ksd_test.md)

## Author

Max Moldovan, <max.moldovan@adelaide.edu.au>

## Examples

``` r
set.seed(1)
x <- matrix(stats::rnorm(200L), ncol = 2L)

# Standard-normal log density (unnormalised): -||x||^2 / 2
ld <- function(z) -0.5 * rowSums(z^2)
s <- numeric_score(ld)

# Matches the closed-form standard-normal score -x to finite-difference order
max(abs(s(x) - (-x)))
#> [1] 5.821565e-12

# Use directly in a goodness-of-fit test
ksd_test(x, score = numeric_score(ld), n_boot = 199L, seed = 1L)
#> 
#>    KSD Test
#> 
#> Statistic: 7.48743e-05 
#> P-value:   0.4000 
#> N:         100 
#> Perms:     199 
#> Kernel X:  imq
#> 
#> Goodness-of-fit verdict
#>   Stein kernel: imq (beta = -0.5)
#>   Bandwidth:     1.57 (median heuristic)
#>   Bootstrap:    wild, B = 199
#>   Surprise:     1.322 bits
#>   Verdict:      consistent with target
#> 
```
