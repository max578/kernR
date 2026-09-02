# Did a kernR result decline to stand behind its own verdict?

`TRUE` when a reliability gate failed (the ESS floor, the proxymix
density-ratio fit-quality gate, or the TACI posterior-adequacy guard)
and the result was stamped with the typed `kernR_abstention` marker.
This is kernR's local mirror of the federation's
`is_orchestra_decline()` predicate
(`ORCHESTRA_dev/integration/refusal_contract.R`), scoped to the member's
own result classes so calling code inside kernR does not need to depend
on the composition-layer contract file.

## Usage

``` r
is_kernR_abstention(x)
```

## Arguments

- x:

  Any object.

## Value

A length-1 logical.

## Examples

``` r
set.seed(42L)
n <- 60L
z <- matrix(rnorm(n * 2L), n, 2L)
x <- z[, 1L] + rnorm(n, sd = 0.5)
y <- 0.5 * x + rnorm(n, sd = 0.5)
res <- suppressWarnings(bd_hsic_test(x, y, z, min_ess_fraction = 0.999))
is_kernR_abstention(res)  # TRUE -- the ESS floor tripped
#> [1] TRUE
```
