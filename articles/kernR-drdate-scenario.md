# DR-DATE for APSIM Scenario Counterfactuals

## The ag-scenario question

Given an APSIM model calibrated on observational data, we want to ask:
*does an intervention (e.g. stubble retention) shift the distribution of
simulated outputs?* The canonical setup is **one PESTO posterior,
forward-simulated under two scenarios** — baseline management and the
intervention. With identical parameter posteriors, the question reduces
to a balanced two-sample distributional comparison on outputs.

[`dr_date_scenario()`](https://max578.github.io/kernR/reference/dr_date_scenario.md)
runs the DR-DATE statistic of Fawkes-Hu-Evans- Sejdinovic (2024) over a
pair of PESTO 0.3.0 `pesto_ensemble_manifest` objects. With identical
parameter columns, the doubly-robust correction collapses to a balanced
MMD on the outputs; with different parameter posteriors, the DR
adjustment absorbs the covariate shift (provided positivity holds).

``` r

library(kernR)
library(PESTO)
```

## One calibration, two scenarios

We use a synthetic linear-Gaussian forward model in place of APSIM so
the vignette runs without an APSIM install. Replace `forward` with
`PESTO::apsim_callback(...)` to drive a real APSIM ensemble.

``` r

npar  <- 2L
nobs  <- 4L
nreal <- 60L
sigma <- 0.05

G <- matrix(rnorm(nobs * npar), nobs, npar)
y_obs <- as.numeric(G %*% c(1.0, -0.5)) + rnorm(nobs, sd = sigma)
names(y_obs) <- paste0("o", seq_len(nobs))

prior <- matrix(rnorm(nreal * npar), nreal, npar,
                dimnames = list(NULL, c("p1", "p2")))
forward <- function(theta) theta %*% t(G)
```

Calibrate **once** with PESTO’s in-process IES driver:

``` r

fit <- pesto_ies_callback(
  forward_model  = forward,
  prior_ensemble = prior,
  obs            = y_obs,
  obs_sd         = sigma,
  noptmax        = 3L,
  verbose        = FALSE
)
```

Now forward-simulate that posterior under two scenarios. The
intervention here is a deterministic shift on the outputs — in a real ag
application this would be e.g. a stubble-retention management rule that
changes APSIM’s predicted yield trajectory.

``` r

par_post   <- as.matrix(fit$par_ensemble[, c("p1", "p2"), with = FALSE])
out_base   <- par_post %*% t(G)
out_intv   <- out_base + 0.6                       # the intervention
colnames(out_base) <- colnames(out_intv) <- names(y_obs)
```

Wrap each scenario as a manifest. **The two manifests share the same
parameter posterior** — that’s the canonical structure.

``` r

real_names <- fit$par_ensemble$real_name

m_baseline <- pesto_ensemble_manifest(
  run_id          = "wagga_baseline_2026",
  params          = data.frame(real_name = real_names, par_post,
                               check.names = FALSE),
  outputs         = data.frame(real_name = real_names, out_base,
                               check.names = FALSE),
  weights         = setNames(rep(1 / sigma, nobs), names(y_obs)),
  obs_target      = setNames(y_obs, names(y_obs)),
  data_hash       = "sha256:vignette_baseline",
  pesto_version   = as.character(packageVersion("PESTO")),
  timestamp       = Sys.time(),
  method          = "ies_callback",
  noptmax         = 3L,
  lambda_schedule = 1
)

m_intervention <- pesto_ensemble_manifest(
  run_id          = "wagga_stubble_2026",
  params          = data.frame(real_name = real_names, par_post,
                               check.names = FALSE),     # SAME params
  outputs         = data.frame(real_name = real_names, out_intv,
                               check.names = FALSE),
  weights         = setNames(rep(1 / sigma, nobs), names(y_obs)),
  obs_target      = setNames(y_obs, names(y_obs)),
  data_hash       = "sha256:vignette_intervention",
  pesto_version   = as.character(packageVersion("PESTO")),
  timestamp       = Sys.time(),
  method          = "ies_callback",
  noptmax         = 3L,
  lambda_schedule = 1
)
m_baseline
#> <pesto_ensemble_manifest> schema 1.0.0
#>   run_id        : wagga_baseline_2026
#>   method        : ies_callback  (noptmax=3)
#>   ensemble      : 60 realisations x 2 parameters | 4 observations
#>   failure rate  : 0.00%
#>   pesto version : 0.6.0.9000  apsim: NA
#>   timestamp     : 2026-06-13T22:24:01+0000
#>   data hash     : sha256:vignette_baseline
```

In production you would build these manifests via
[`as_manifest()`](https://rdrr.io/pkg/PESTO/man/as_manifest.html) on two
separate IES runs *only when the scenarios genuinely produced different
calibration data*. For the “did the intervention shift the forward
outputs?” question, the same-posterior construction above is the right
one.

## Run the DR-DATE scenario test

``` r

res <- dr_date_scenario(
  baseline       = m_baseline,
  intervention   = m_intervention,
  n_permutations = 299L,
  seed           = 1L
)
print(res)
#> 
#>    DR-DATE (scenario) Test
#> 
#> Statistic: 0.788964 
#> P-value:   0.0033 
#> N:         120 
#> Perms:     299 
#> Kernel Y:  rbf (bw = 1.199)
#> ESS:       59.4 
#> 
#> Scenario contrast
#>   baseline      : wagga_baseline_2026 (n=60)
#>   intervention  : wagga_stubble_2026 (n=60)
#>   outputs tested: o1, o2, o3, o4
#>   PESTO versions: baseline=0.6.0.9000, intervention=0.6.0.9000
#>   fidelity      : baseline=single, intervention=single
#>   Verdict:        REJECT (distributions differ; intervention has effect)
```

With a clear distributional shift, the test should reject and return a
low `p_value`. The `Verdict` line gives the directly actionable read.

## Null case — same outputs

If the intervention has no effect (shift = 0), the two output
distributions are identical and the test should fail to reject:

``` r

out_null <- out_base                       # no intervention effect
m_null <- pesto_ensemble_manifest(
  run_id          = "wagga_baseline_replicate",
  params          = data.frame(real_name = real_names, par_post,
                               check.names = FALSE),
  outputs         = data.frame(real_name = real_names, out_null,
                               check.names = FALSE),
  weights         = setNames(rep(1 / sigma, nobs), names(y_obs)),
  obs_target      = setNames(y_obs, names(y_obs)),
  data_hash       = "sha256:vignette_null",
  pesto_version   = as.character(packageVersion("PESTO")),
  timestamp       = Sys.time(),
  method          = "ies_callback",
  noptmax         = 3L,
  lambda_schedule = 1
)
res_null <- dr_date_scenario(
  baseline       = m_baseline,
  intervention   = m_null,
  n_permutations = 299L,
  seed           = 2L
)
res_null$p_value
#> [1] 0.96
```

## What gets validated

[`dr_date_scenario()`](https://max578.github.io/kernR/reference/dr_date_scenario.md)
hard-stops if the two manifests don’t agree on:

- parameter column names (you can’t compare incompatible priors),
- observation column names (you can’t compare different outputs),
- PESTO major.minor version (forward / backward incompatibility safety).

This is by design — silent comparison of incompatible scenarios is a
worse failure mode than a noisy abort.

## Output sub-selection

For ensembles with many output columns, focus the test on the outputs of
scientific interest with `output =`:

``` r

dr_date_scenario(
  m_baseline, m_intervention,
  output         = c("o1", "o3"),     # subset
  n_permutations = 99L, seed = 3L
)$outputs_tested
#> [1] "o1" "o3"
```

## When parameter posteriors genuinely differ

If the two scenarios came from *different calibration data* (so the two
PESTO runs produced different posteriors), then `params` will differ
between the two manifests. DR-DATE will then run its full doubly-robust
correction, with one important caveat: it assumes **positivity** — every
parameter region with non-zero baseline density also has non-zero
intervention density (and vice versa). If the two posteriors are
completely separated in parameter space, the test will see “perfect
separation” in the propensity-model fit and become uninformative
(`p_value ~ 1`). For that regime, use
[`PESTO::pesto_ies_callback()`](https://rdrr.io/pkg/PESTO/man/pesto_ies_callback.html)
with a shared prior + overlapping calibration data, or restrict to
outputs whose causal pathway is independent of the separating
parameters.

## Where the cross-package plumbing lives

- Forward-model run:
  [`PESTO::pesto_ies_callback()`](https://rdrr.io/pkg/PESTO/man/pesto_ies_callback.html)
  (in-process R callback) or
  [`PESTO::apsim_callback()`](https://rdrr.io/pkg/PESTO/man/apsim_callback.html)
  (apsimx adapter).
- Run capture:
  [`PESTO::as_manifest()`](https://rdrr.io/pkg/PESTO/man/as_manifest.html)
  → `pesto_ensemble_manifest`.
- Persistence:
  [`PESTO::write_manifest()`](https://rdrr.io/pkg/PESTO/man/write_manifest.html)
  /
  [`PESTO::read_manifest()`](https://rdrr.io/pkg/PESTO/man/read_manifest.html)
  /
  [`PESTO::verify_manifest()`](https://rdrr.io/pkg/PESTO/man/verify_manifest.html).
- Distributional verdict: this function.
- (Future) proxymix density-ratio backend plug-in: tracked under §C1 of
  the roadmap; will become a fourth propensity-model option once it
  lands.

## Reproducibility

``` r

sessionInfo()
#> R version 4.6.0 (2026-04-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] PESTO_0.6.0.9000 kernR_0.8.0     
#> 
#> loaded via a namespace (and not attached):
#>  [1] vctrs_0.7.3        cli_3.6.6          knitr_1.51         rlang_1.2.0       
#>  [5] xfun_0.58          otel_0.2.0         generics_0.1.4     S7_0.2.2          
#>  [9] textshaping_1.0.5  jsonlite_2.0.0     data.table_1.18.4  glue_1.8.1        
#> [13] htmltools_0.5.9    ragg_1.5.2         sass_0.4.10        scales_1.4.0      
#> [17] rmarkdown_2.31     grid_4.6.0         evaluate_1.0.5     jquerylib_0.1.4   
#> [21] fastmap_1.2.0      yaml_2.3.12        lifecycle_1.0.5    compiler_4.6.0    
#> [25] RColorBrewer_1.1-3 fs_2.1.0           Rcpp_1.1.1-1.1     farver_2.1.2      
#> [29] systemfonts_1.3.2  digest_0.6.39      R6_2.6.1           bslib_0.11.0      
#> [33] gtable_0.3.6       tools_4.6.0        ggplot2_4.0.3      pkgdown_2.2.0     
#> [37] cachem_1.1.0       desc_1.4.3
```
