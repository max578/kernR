# Density-Ratio Backends and the proxymix Binding

## The cross-package contract

[`kernR::estimate_density_ratio()`](https://max578.github.io/kernR/reference/estimate_density_ratio.md)
is the entry point for every kernR estimator that needs to reweight
observational samples to an interventional distribution — most notably
[`bd_hsic_test()`](https://max578.github.io/kernR/reference/bd_hsic_test.md)
for the backdoor-HSIC causal test. The function offers four
density-ratio backends behind a single signature:

| Backend | Family | When to use |
|----|----|----|
| `logistic` | Noise-contrastive classifier | Default; robust on smooth, unimodal densities |
| `ranger` | Random-forest classifier | Flexible non-linear NCE; needs `ranger` |
| `xgboost` | Gradient-boosted classifier | Strong on tabular interactions; needs `xgboost` |
| `proxymix` | Gaussian-mixture density ratio | Multimodal/skewed densities; parametric alternative |

The `proxymix` backend is the cross-package wedge between kernR (the
distributional verdict layer of the UQ ag stack) and proxymix (the
Gaussian-mixture proxy / KL-density-ratio bridge, Hoek & Elliott 2024).
It fits one GMM to the joint sample cloud `(x, z)`, a second GMM to the
product-of-marginals cloud `(x_perm, z)`, then evaluates the analytic
ratio of the two mixture densities at each observation. No classifier
calibration step — the ratio is closed-form in the fitted parameters.

## Four backends, one problem

A toy confounded design: `z` is a 2-D Gaussian confounder; `x` is a
linear-Gaussian function of `z`; `y` carries a real causal effect from
`x` plus a confounded path through `z`.

``` r

suppressPackageStartupMessages(library(kernR))
set.seed(2026L)
n <- 200L
z <- matrix(rnorm(n * 2L), n, 2L)
x <- z[, 1L] + rnorm(n, sd = 0.5)
y <- 0.7 * x + z[, 2L] + rnorm(n, sd = 0.4)
```

Each backend produces a `density_ratio_fit` object exposing ESS and a
weight vector. We tabulate them side-by-side.

``` r

dr_logistic <- estimate_density_ratio(x, z, method = "logistic", seed = 1L)
dr_ranger <- if (requireNamespace("ranger", quietly = TRUE)) {
  estimate_density_ratio(x, z, method = "ranger", seed = 1L)
} else NULL
dr_xgb <- if (requireNamespace("xgboost", quietly = TRUE)) {
  estimate_density_ratio(x, z, method = "xgboost", seed = 1L)
} else NULL
```

``` r

dr_proxymix <- estimate_density_ratio(
  x, z,
  method = "proxymix",
  proxymix_components = 2L,
  seed = 1L
)
```

| backend  |   ess | min_weight | max_weight |
|:---------|------:|-----------:|-----------:|
| logistic | 199.4 |     0.8860 |       1.20 |
| ranger   | 156.0 |     0.2380 |       3.21 |
| xgboost  |  81.3 |     0.0918 |       7.08 |
| proxymix |  14.3 |     0.0171 |      47.00 |

Density-ratio diagnostics across backends. {.table}

## bd-HSIC under each backend

``` r

res_logistic <- bd_hsic_test(
  x, y, z, density_ratio = "logistic",
  n_permutations = 199L, seed = 1L
)
```

``` r

res_proxymix <- bd_hsic_test(
  x, y, z, density_ratio = "proxymix",
  n_permutations = 199L, seed = 1L
)
#> Warning: bd_hsic_test(): ESS (3.4) is below 10% of n_test (100). The weighted
#> test statistic is dominated by a small number of high-weight observations; the
#> resulting p-value is not a reliable verdict. Increase n, switch density_ratio
#> backend, or tighten the design.
```

| backend  | statistic | p_value |  ess |
|:---------|----------:|--------:|-----:|
| logistic |    0.0211 |   0.155 | 99.9 |
| proxymix |    0.0062 |   0.320 |  3.4 |

bd-HSIC test under each density-ratio backend. {.table}

## When to reach for proxymix

The classifier-based backends (`logistic`, `ranger`, `xgboost`) are the
default for a reason: they tolerate misspecified densities, scale to
high-dimensional `z`, and have well-understood calibration. Reach for
`method = "proxymix"` when:

- the joint density is plausibly **multimodal** (multi-regime climate,
  paddock × variety designs with distinct production zones, animal
  cohorts with separable subpopulations) — a 2-component GMM represents
  this cleanly where a classifier would smear across the modes;
- you need a **parametric** density-ratio whose components you can
  inspect, hand off to a downstream Bayesian step (via
  [`proxymix::gmm_target_from_posterior()`](https://rdrr.io/pkg/proxymix/man/gmm_target_from_posterior.html)),
  or use as the seed of a KLD-EM refinement on a target you can evaluate
  but not sample from;
- classifier calibration is unreliable on the cohort at hand (small `n`,
  sharp class imbalance after the joint-vs-marginal split, pathological
  feature scaling).

The `proxymix` package is GRDC-firewalled (MIT, no GRDC IP flows in) and
ships its full Gaussian-mixture proxy API independently of kernR. kernR
consumes it as a soft dependency via
[`requireNamespace()`](https://rdrr.io/r/base/ns-load.html) — the
binding is one-way and rebuildable from the local `proxymix_*.tar.gz`
source.

## References

- Hoek, J. van der & Elliott, R. J. (2024). *Mixtures of multivariate
  Gaussians.* Stochastic Analysis and Applications. DOI:
  10.1080/07362994.2024.2372605.
- Hu, R., Sejdinovic, D. & Evans, R. J. (2024). A kernel test for causal
  association via noise contrastive backdoor adjustment. *Journal of
  Machine Learning Research*, 25(160), 1–56.

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
#> [1] kernR_0.3.1
#> 
#> loaded via a namespace (and not attached):
#>  [1] vctrs_0.7.3        cli_3.6.6          knitr_1.51         rlang_1.2.0       
#>  [5] xfun_0.57          textshaping_1.0.5  S7_0.2.2           jsonlite_2.0.0    
#>  [9] data.table_1.18.4  glue_1.8.1         ranger_0.18.0      proxymix_0.3.0    
#> [13] htmltools_0.5.9    PESTO_0.4.1        ragg_1.5.2         sass_0.4.10       
#> [17] scales_1.4.0       rmarkdown_2.31     grid_4.6.0         evaluate_1.0.5    
#> [21] jquerylib_0.1.4    fastmap_1.2.0      yaml_2.3.12        lifecycle_1.0.5   
#> [25] compiler_4.6.0     mvnfast_0.2.8      RColorBrewer_1.1-3 fs_2.1.0          
#> [29] Rcpp_1.1.1-1.1     lattice_0.22-9     farver_2.1.2       systemfonts_1.3.2 
#> [33] digest_0.6.39      xgboost_3.2.1.1    R6_2.6.1           Matrix_1.7-5      
#> [37] bslib_0.11.0       withr_3.0.2        gtable_0.3.6       tools_4.6.0       
#> [41] ggplot2_4.0.3      pkgdown_2.2.0      cachem_1.1.0       desc_1.4.3
```
