# Emit a kernR verdict as an orchestra manifest

The federation's emit / adapt generic. Both methods set
`inferential_target = "treatment_effects"`: kernR's verdicts are
hypothesis-test or mechanism-consistency calls about a treatment effect,
never a point estimate or a prediction. `abstained` in the typed
`summary` is set from the result's own reliability flags –
`posterior_adequacy$ok` for a TACI verdict, `ess_warning` /
`density_ratio_warning` for a kernel test – so a downstream consumer
reading only the manifest sees the same reliability signal the raw
result already carried, without having to know which flag belongs to
which verb.

## Usage

``` r
as_orchestra_manifest(x, ...)

# S3 method for class 'taci_result'
as_orchestra_manifest(x, ..., run_id = NULL)

# S3 method for class 'kernel_test_result'
as_orchestra_manifest(x, ..., run_id = NULL)
```

## Arguments

- x:

  A `taci_result` or `kernel_test_result` object.

- ...:

  Method-specific arguments.

- run_id:

  Optional character run identifier; a content hash is derived when
  omitted, so two identical verdicts get the same identifier.

## Value

An `orchestra_manifest` S7 object.

## See also

[`verify_manifest()`](https://max578.github.io/kernR/reference/verify_manifest.md),
[`taci_test()`](https://max578.github.io/kernR/reference/taci_test.md)

Other orchestra manifest:
[`orchestra_manifest()`](https://max578.github.io/kernR/reference/orchestra_manifest.md),
[`verify_manifest()`](https://max578.github.io/kernR/reference/verify_manifest.md)

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
m@summary$headline
#> [1] "mechanism_consistent_effect [unverified]"
m@summary$abstained  # FALSE: the posterior-adequacy gate passed
#> [1] FALSE
```
