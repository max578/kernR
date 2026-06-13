# Tidy a kernR Test Result Into a One-Row-Per-Term Data Frame

Turns a kernR verdict object into a flat, `broom`-style summary
`data.frame`: one row per tested term, with stable, documented columns.
The method is registered against the
[`tidy()`](https://generics.r-lib.org/reference/tidy.html) generic
(re-exported by `broom`), so `broom::tidy(result)` and
`generics::tidy(result)` both dispatch here.

## Usage

``` r
# S3 method for class 'kernel_test_result'
tidy(x, ...)
```

## Arguments

- x:

  A `kernel_test_result` (the class returned by
  [`hsic_test()`](https://max578.github.io/kernR/reference/hsic_test.md),
  [`mmd_test()`](https://max578.github.io/kernR/reference/mmd_test.md),
  [`bd_hsic_test()`](https://max578.github.io/kernR/reference/bd_hsic_test.md),
  and friends) or a subclass such as `mmd_ppc`.

- ...:

  Currently unused; present for generic compatibility.

## Value

A one-row `data.frame` with columns:

- term:

  Character. The quantity tested – the test `method` (for example
  `"bd-HSIC"`, `"MMD"`).

- statistic:

  Numeric. The observed test statistic.

- p.value:

  Numeric. The permutation p-value. Note the `broom`-canonical dot,
  distinct from the result's native `p_value` field.

- n:

  Integer. The sample size the statistic was computed on.

- n_permutations:

  Integer. Number of permutations in the null.

- ess:

  Numeric. Effective sample size of the importance weights, or `NA` for
  unweighted tests.

Subclasses contribute extra columns where they carry extra fields: an
`mmd_ppc` result adds `surprise_bits` and `reject`.

## Details

This is the supported accessor for downstream code. Reaching into the
result list by field name is fragile: kernR's native field is `p_value`
(with an underscore), whereas the `broom` convention – followed here –
is `p.value` (with a dot). Use
[`tidy()`](https://generics.r-lib.org/reference/tidy.html) and read
`p.value` rather than guessing at the raw field name.

## See also

[`hsic_test()`](https://max578.github.io/kernR/reference/hsic_test.md),
[`mmd_test()`](https://max578.github.io/kernR/reference/mmd_test.md),
[`bd_hsic_test()`](https://max578.github.io/kernR/reference/bd_hsic_test.md),
[`mmd_ppc()`](https://max578.github.io/kernR/reference/mmd_ppc.md)

## Examples

``` r
set.seed(1)
x <- matrix(rnorm(200L), 100L, 2L)
y <- x[, 1L] + rnorm(100L)
res <- hsic_test(x, y, n_permutations = 199L, seed = 1L)
tidy(res)
#>   term  statistic p.value   n n_permutations ess
#> 1 HSIC 0.01021759   0.005 100            199  NA
```
