# Causal Association Testing with bd-HSIC

## The Problem

Standard independence tests cannot distinguish *causal* association from
*confounded* association. If treatment X and outcome Y share a common
cause Z, they will appear dependent even if X has no causal effect on Y.

The **bd-HSIC** test (Hu, Sejdinovic & Evans, 2024) solves this by
testing the *do-null* hypothesis:

> H_0: p(y \| do(x)) = p\*(y) for all x

This uses Pearl’s do-operator: after intervening on X, is Y still
associated with X?

## How It Works

1.  **Density ratio estimation**: Estimate w(x, z) = p\*(x) / p(x\|z) to
    reweight observational samples to the interventional distribution.
2.  **Weighted HSIC**: Compute HSIC between X and Y under the reweighted
    (interventional) distribution.
3.  **Cluster-based permutation**: Obtain p-values by permuting Y within
    clusters of similar conditional densities p(x\|z).

## Example: Linear Causal Effect

``` r

library(kernR)
set.seed(42)

n <- 300
z <- matrix(rnorm(n * 2), n, 2)
x <- 0.5 * z[, 1] + rnorm(n)           # X depends on confounder Z
y <- 0.8 * x + 0.5 * z[, 2] + rnorm(n) # Y depends causally on X and on Z

result <- bd_hsic_test(x, y, z,
  n_permutations = 200,
  seed = 1
)
result
#> 
#>    bd-HSIC Test
#> 
#> Statistic: 0.0182464 
#> P-value:   0.0050 
#> N:         150 
#> Perms:     200 
#> Kernel X:  rbf (bw = 1.062)
#> Kernel Y:  rbf (bw = 1.425)
#> ESS:       149.7
```

The test detects the causal association between X and Y.

## Example: No Causal Effect (Confounding Only)

``` r

set.seed(42)
n <- 300
z <- matrix(rnorm(n * 2), n, 2)
x <- 0.5 * z[, 1] + rnorm(n)
y <- 0.5 * z[, 1] + z[, 2] + rnorm(n) # Y depends on Z, not on X

result_null <- bd_hsic_test(x, y, z,
  n_permutations = 200,
  seed = 1
)
result_null
#> 
#>    bd-HSIC Test
#> 
#> Statistic: 0.0022688 
#> P-value:   0.4826 
#> N:         150 
#> Perms:     200 
#> Kernel X:  rbf (bw = 1.062)
#> Kernel Y:  rbf (bw = 1.553)
#> ESS:       149.7
```

The large p-value correctly indicates no causal effect.

## Example: Non-Linear Causal Effect

A key advantage of bd-HSIC is detecting non-linear effects that linear
methods (PDS, Double ML) completely miss:

``` r

set.seed(42)
n <- 400
z <- matrix(rnorm(n * 2), n, 2)
x <- z[, 1] + rnorm(n)
y <- x^2 + z[, 2] + rnorm(n, sd = 0.5) # Quadratic causal effect

result_nl <- bd_hsic_test(x, y, z,
  n_permutations = 200,
  seed = 1
)
result_nl
#> 
#>    bd-HSIC Test
#> 
#> Statistic: 0.0206771 
#> P-value:   0.0050 
#> N:         200 
#> Perms:     200 
#> Kernel X:  rbf (bw = 1.299)
#> Kernel Y:  rbf (bw = 2.029)
#> ESS:       199.6
```

## Diagnostic: Null Distribution

``` r

plot(result)
```

![bd-HSIC permutation null
distribution.](kernR-bdhsic_files/figure-html/plot-null-1.png)

bd-HSIC permutation null distribution.

## Using the Formula Interface

``` r

dat <- data.frame(y = y, x = x, z1 = z[, 1], z2 = z[, 2])
result_f <- kernel_causal_test(
  y ~ x | z1 + z2,
  data = dat,
  method = "bd-hsic",
  n_permutations = 100,
  seed = 1
)
result_f
#> 
#>    bd-HSIC Test
#> 
#> Statistic: 0.0206771 
#> P-value:   0.0099 
#> N:         200 
#> Perms:     100 
#> Kernel X:  rbf (bw = 1.299)
#> Kernel Y:  rbf (bw = 2.029)
#> ESS:       199.6
```

## When to Use bd-HSIC

| Scenario | Use bd-HSIC? |
|----|----|
| Testing if X causally affects Y (adjusting for Z) | Yes |
| Non-linear or non-monotone causal effects | Yes – key advantage |
| Continuous, binary, or mixed treatments | Yes |
| Very high-dimensional confounders | Consider using `density_ratio = "ranger"` |
| Extremely strong confounding | Caution: density ratio estimation may fail |

## References

- Hu, R., Sejdinovic, D., & Evans, R. J. (2024). A kernel test for
  causal association via noise contrastive backdoor adjustment. *JMLR*,
  25(160), 1-56.

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
#> [1] kernR_0.7.0.9000
#> 
#> loaded via a namespace (and not attached):
#>  [1] vctrs_0.7.3        cli_3.6.6          knitr_1.51         rlang_1.2.0       
#>  [5] xfun_0.58          S7_0.2.2           textshaping_1.0.5  jsonlite_2.0.0    
#>  [9] data.table_1.18.4  glue_1.8.1         htmltools_0.5.9    PESTO_0.6.0.9000  
#> [13] ragg_1.5.2         sass_0.4.10        scales_1.4.0       rmarkdown_2.31    
#> [17] grid_4.6.0         evaluate_1.0.5     jquerylib_0.1.4    fastmap_1.2.0     
#> [21] yaml_2.3.12        lifecycle_1.0.5    compiler_4.6.0     RColorBrewer_1.1-3
#> [25] fs_2.1.0           Rcpp_1.1.1-1.1     farver_2.1.2       systemfonts_1.3.2 
#> [29] digest_0.6.39      R6_2.6.1           bslib_0.11.0       gtable_0.3.6      
#> [33] tools_4.6.0        pkgdown_2.2.0      ggplot2_4.0.3      cachem_1.1.0      
#> [37] desc_1.4.3
```
