# Verify an orchestra manifest's payload integrity

Recomputes the payload hash from the object's own load-bearing slots and
compares it to the stored `data_hash`. A mismatch means the payload was
modified after emission. Matches the reference contract's
`verify_manifest()`.

## Usage

``` r
verify_manifest(m)
```

## Arguments

- m:

  An `orchestra_manifest` object.

## Value

A list with logical `ok` and a human-readable `message`.

## See also

[`as_orchestra_manifest()`](https://max578.github.io/kernR/reference/as_orchestra_manifest.md)

Other orchestra manifest:
[`as_orchestra_manifest()`](https://max578.github.io/kernR/reference/as_orchestra_manifest.md),
[`orchestra_manifest()`](https://max578.github.io/kernR/reference/orchestra_manifest.md)

## Examples

``` r
mechanism <- function(theta, X, t) theta[3] + theta[1] * (1 - exp(-theta[2] * t))
set.seed(1L)
n <- 60L
nrate <- runif(n, 0, 200)
yield <- 1.1 + 4.2 * (1 - exp(-0.018 * nrate)) + rnorm(n, 0, 0.25)
post <- cbind(ymax = rnorm(100L, 4.2, 0.30),
             rate = rnorm(100L, 0.018, 0.004),
             y0   = rnorm(100L, 1.1, 0.15))
verdict <- taci_test(post, mechanism, X = matrix(1, n, 1),
                     treatment = nrate, outcome = yield,
                     n_perm = 99L, seed = 1L)
m <- as_orchestra_manifest(verdict)
verify_manifest(m)$ok  # TRUE -- payload untampered
#> Error: Can't find method for `verify_manifest(<orchestra_manifest>)`.
```
