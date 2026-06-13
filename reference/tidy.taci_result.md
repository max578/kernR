# Tidy a TACI Mechanism-Consistency Result

Turns a
[`taci_test()`](https://max578.github.io/kernR/reference/taci_test.md)
verdict into a flat, `broom`-style one-row summary. Registered against
the [`tidy()`](https://generics.r-lib.org/reference/tidy.html) generic
(re-exported by `broom`), so `broom::tidy()` and
[`generics::tidy()`](https://generics.r-lib.org/reference/tidy.html)
both dispatch here.

## Usage

``` r
# S3 method for class 'taci_result'
tidy(x, ...)
```

## Arguments

- x:

  A `taci_result` object from
  [`taci_test()`](https://max578.github.io/kernR/reference/taci_test.md).

- ...:

  Currently unused; present for generic compatibility.

## Value

A one-row `data.frame` with columns:

- term:

  Character. Always `"taci"`.

- statistic:

  Numeric. The observed bd-HSIC statistic.

- p.value:

  Numeric. The H0 tail p-value (`p_h0`).

- decision:

  Character. The three-way decision enum (`"no_effect"`,
  `"mechanism_consistent_effect"`, `"mechanism_inconsistent_effect"`).

- grounding:

  Character. `"grounded"` when the mechanism's provenance was declared,
  else `"[unverified]"` (Independent Oracle Principle).

- h1_percentile:

  Numeric. Where the observed statistic sits in the model-implied H1
  band, in `[0, 1]`.

- borderline:

  Logical. Whether the consistency label is fragile.

- n:

  Integer. The sample size.

## Details

The `p.value` column carries the broom-canonical dot and reports the H0
tail p-value (`p_h0`). The `grounding` column threads through the
Independent Oracle Principle label, so a downstream summary cannot
silently drop the fact that a verdict built on an un-declared mechanism
is `"[unverified]"`.

## See also

[`taci_test()`](https://max578.github.io/kernR/reference/taci_test.md)

## Examples

``` r
set.seed(1)
n <- 60L
nrate <- runif(n, 0, 200)
yield <- 1.1 + 4.2 * (1 - exp(-0.018 * nrate)) + rnorm(n, 0, 0.25)
post <- cbind(ymax = rnorm(150L, 4.2, 0.30),
              rate = rnorm(150L, 0.018, 0.004),
              y0   = rnorm(150L, 1.1, 0.15))
mitscherlich <- function(theta, X, t) {
  theta[3L] + theta[1L] * (1 - exp(-theta[2L] * t))
}
res <- taci_test(post, mitscherlich, X = matrix(1, n, 1L),
                 treatment = nrate, outcome = yield,
                 n_perm = 80L, seed = 1L)
tidy(res)
#>   term  statistic    p.value                    decision    grounding
#> 1 taci 0.06587821 0.01234568 mechanism_consistent_effect [unverified]
#>   h1_percentile borderline  n
#> 1         0.425      FALSE 60
```
