# Hierarchical and Nested Data

## Why Hierarchy Matters

Many real-world datasets have nested structure:

- Patients within hospitals
- Students within schools
- Plots within farms
- Repeated measures within subjects

Standard kernel tests assume independent observations. When observations
are clustered, **within-cluster correlation inflates type I error**.
[`hierarchical_test()`](https://max578.github.io/kernR/reference/hierarchical_test.md)
accounts for this by decomposing the test statistic and permuting within
clusters.

## Example: Agriculture Trial

Imagine a randomised fertiliser trial across 20 farms, each with 30
plots.

``` r

library(kernR)
set.seed(42)

n_farms <- 20
n_plots <- 30
n <- n_farms * n_plots
farm_id <- rep(1:n_farms, each = n_plots)

# Farm-level random effects
farm_effect <- rnorm(n_farms, sd = 2)[farm_id]
soil <- matrix(rnorm(n * 2), n, 2)

# Treatment assignment (partially confounded by soil)
treatment <- rbinom(n, 1, plogis(0.3 * soil[, 1]))

# Yield: treatment has a real effect + farm random effect
yield <- 0.8 * treatment + farm_effect + 0.5 * soil[, 1] + rnorm(n)

result <- hierarchical_test(
  y = yield,
  treatment = treatment,
  covariates = soil,
  cluster_id = farm_id,
  method = "dr-date",
  n_permutations = 100,
  weight_method = "icc",
  seed = 1
)
result
#> 
#>    Hierarchical-DRDATE Test
#> 
#> Statistic: 0.00874629 
#> P-value:   0.0099 
#> N:         600 
#> Perms:     100 
#> Kernel Y:  rbf (bw = 2.642)
```

The test detects the treatment effect while correctly accounting for the
farm-level clustering.

## Decomposition: Within vs Between

The test provides both components:

``` r

cat("Within-cluster average statistic:",
  mean(result$hierarchical$within_stats, na.rm = TRUE), "\n")
#> Within-cluster average statistic: 0.09414643
cat("Between-cluster statistic:",
  result$hierarchical$between_stat, "\n")
#> Between-cluster statistic: -0.004092593
cat("Combined statistic:",
  result$statistic, "\n")
#> Combined statistic: 0.008746285
cat("Weight method:", result$hierarchical$weight_method, "\n")
#> Weight method: icc
```

## Weight Methods

| Method          | Behaviour                                                |
|-----------------|----------------------------------------------------------|
| `"equal"`       | Equal weight to within and between components            |
| `"icc"`         | Weight by ICC (more between-weight when clusters differ) |
| `"within_only"` | Ignore between-cluster variation                         |

## Example: No Treatment Effect

``` r

set.seed(42)
yield_null <- farm_effect + 0.5 * soil[, 1] + rnorm(n) # No treatment effect

result_null <- hierarchical_test(
  y = yield_null,
  treatment = treatment,
  covariates = soil,
  cluster_id = farm_id,
  method = "dr-date",
  n_permutations = 100,
  seed = 1
)
cat("P-value under null:", result_null$p_value, "\n")
#> P-value under null: 0.3663366
```

## When to Use Hierarchical Tests

| Scenario | Recommendation |
|----|----|
| Independent observations | Use standard [`dr_date_test()`](https://max578.github.io/kernR/reference/dr_date_test.md) / [`bd_hsic_test()`](https://max578.github.io/kernR/reference/bd_hsic_test.md) |
| Clustered data (known groups) | Use [`hierarchical_test()`](https://max578.github.io/kernR/reference/hierarchical_test.md) with `cluster_id` |
| Few large clusters | `weight_method = "icc"` |
| Many small clusters | `weight_method = "equal"` |
| Unsure about between-cluster effects | `weight_method = "within_only"` (conservative) |

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
#> [1] kernR_0.3.0
#> 
#> loaded via a namespace (and not attached):
#>  [1] vctrs_0.7.3        cli_3.6.6          knitr_1.51         rlang_1.2.0       
#>  [5] xfun_0.57          S7_0.2.2           textshaping_1.0.5  jsonlite_2.0.0    
#>  [9] data.table_1.18.4  glue_1.8.1         htmltools_0.5.9    PESTO_0.4.1       
#> [13] ragg_1.5.2         sass_0.4.10        scales_1.4.0       rmarkdown_2.31    
#> [17] grid_4.6.0         evaluate_1.0.5     jquerylib_0.1.4    fastmap_1.2.0     
#> [21] yaml_2.3.12        lifecycle_1.0.5    compiler_4.6.0     RColorBrewer_1.1-3
#> [25] fs_2.1.0           Rcpp_1.1.1-1.1     farver_2.1.2       systemfonts_1.3.2 
#> [29] digest_0.6.39      R6_2.6.1           bslib_0.11.0       gtable_0.3.6      
#> [33] tools_4.6.0        pkgdown_2.2.0      ggplot2_4.0.3      cachem_1.1.0      
#> [37] desc_1.4.3
```
