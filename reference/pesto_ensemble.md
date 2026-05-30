# PESTO Ensemble Manifest (Constructor)

Lightweight constructor for a posterior-predictive ensemble produced by
PESTO (or any compatible UQ engine). Bundles a posterior-predictive
sample matrix with optional held-out observations and free-form
metadata, providing a stable interface for
[`mmd_ppc()`](https://max578.github.io/kernR/reference/mmd_ppc.md) and
other downstream verdict-layer tools.

## Usage

``` r
pesto_ensemble(posterior, observed = NULL, metadata = list())
```

## Arguments

- posterior:

  Numeric matrix `M x d`: `M` posterior-predictive draws over `d` output
  dimensions (e.g. yield, biomass).

- observed:

  Optional numeric matrix `n_obs x d` of held-out observations. May be
  `NULL` when observations are supplied later to
  [`mmd_ppc()`](https://max578.github.io/kernR/reference/mmd_ppc.md).

- metadata:

  Optional named list of free-form metadata (run id, ensemble seed,
  holdout year, etc.).

## Value

An object of class `"pesto_ensemble"`.

## Details

This is the kernR-side schema for the cross-package contract; until
PESTO ships its native manifest emitter, callers can construct the
object directly from in-memory matrices.

## See also

[`mmd_ppc()`](https://max578.github.io/kernR/reference/mmd_ppc.md)

## Examples

``` r
set.seed(1)
post <- matrix(stats::rnorm(200L), ncol = 2L)
obs  <- matrix(stats::rnorm(20L),  ncol = 2L)
ens  <- pesto_ensemble(post, obs, metadata = list(holdout_year = 2018))
ens
#> 
#>   PESTO ensemble manifest
#> 
#> Posterior: 100 draws x 2 dims
#> Observed:  10 obs x 2 dims
#> Metadata:  holdout_year
#> 
```
