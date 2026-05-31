# Distributional Treatment Effect Tests (DR-DATE / DR-DETT)

## Beyond Mean Effects

Standard causal inference methods (Double ML, TMLE) test whether
treatment shifts the *mean* outcome. But many real treatment effects are
**distributional** – they change variance, shape, or modality without
necessarily changing the mean.

**DR-DATE** and **DR-DETT** (Fawkes, Hu, Evans & Sejdinovic, 2024) test
for *any* distributional difference between Y(1) and Y(0), using doubly
robust kernel embeddings.

## Key Concepts

- **DR-DATE** (Distributional Average Treatment Effect): Tests whether
  P(Y(1)) = P(Y(0)) over the entire population.
- **DR-DETT** (Distributional Effect on the Treated): Tests whether
  P(Y(1)\|T=1) = P(Y(0)\|T=1), focusing on the treated subgroup.
  Requires only one-sided overlap.
- **Double robustness**: Consistent if either the propensity model or
  the outcome model is correctly specified.

## Example: Mean Shift (Detectable by All Methods)

``` r

library(kernR)
set.seed(42)

n <- 300
x <- matrix(rnorm(n * 2), n, 2)
logit_p <- 0.3 * x[, 1] - 0.2 * x[, 2]
t <- rbinom(n, 1, plogis(logit_p))
y <- t * 1.0 + 0.5 * x[, 1] + rnorm(n, sd = 0.5) # Mean shift of 1.0

result <- dr_date_test(y, t, x,
  n_permutations = 200,
  seed = 1
)
result
#> 
#>    DR-DATE Test
#> 
#> Statistic: 0.37054 
#> P-value:   0.0050 
#> N:         300 
#> Perms:     200 
#> Kernel Y:  rbf (bw = 0.8972)
#> ESS:       139.7
```

## Example: Variance Effect Only (Invisible to Mean-Based Tests)

This is where DR-DATE shines. The treatment changes the *variance* of
the outcome but not the mean – DML and TMLE would have zero power here.

``` r

set.seed(42)
n <- 400
x <- matrix(rnorm(n * 2), n, 2)
t <- rbinom(n, 1, plogis(0.3 * x[, 1]))

# Treatment doubles the variance but does NOT shift the mean
y <- (1 - t) * rnorm(n, sd = 1) + t * rnorm(n, sd = 2.5) + 0.5 * x[, 1]

cat("Mean difference:", mean(y[t == 1]) - mean(y[t == 0]), "\n")
#> Mean difference: 0.5228389
cat("SD treated:", sd(y[t == 1]), "  SD control:", sd(y[t == 0]), "\n")
#> SD treated: 2.658814   SD control: 1.076454

result_var <- dr_date_test(y, t, x,
  n_permutations = 200,
  outcome_model = "zero",
  seed = 1
)
result_var
#> 
#>    DR-DATE Test
#> 
#> Statistic: 0.160419 
#> P-value:   0.0050 
#> N:         400 
#> Perms:     200 
#> Kernel Y:  rbf (bw = 1.684)
#> ESS:       176.9
```

## DR-DETT: Effect on the Treated

When overlap is imperfect (some covariate regions have nearly all
treated or all control units), DR-DETT is more robust because it
requires only one-sided overlap.

``` r

set.seed(42)
n <- 300
x <- matrix(rnorm(n * 2), n, 2)
t <- rbinom(n, 1, plogis(0.5 * x[, 1]))
y <- t * rnorm(n, mean = 0.5, sd = 1.5) + (1 - t) * rnorm(n) + x[, 1]

result_dett <- dr_dett_test(y, t, x,
  n_permutations = 200,
  seed = 1
)
result_dett
#> 
#>    DR-DETT Test
#> 
#> Statistic: 0.0194285 
#> P-value:   0.0597 
#> N:         300 
#> Perms:     200 
#> Kernel Y:  rbf (bw = 1.557)
#> ESS:       102.2
```

## Comparing the Tests

``` r

cat("DR-DATE p-value:", result_var$p_value, "\n")
#> DR-DATE p-value: 0.004975124
cat("DR-DETT p-value:", result_dett$p_value, "\n")
#> DR-DETT p-value: 0.05970149
```

## Using the Formula Interface

``` r

dat <- data.frame(y = y, treatment = t, x1 = x[, 1], x2 = x[, 2])
result_f <- kernel_causal_test(
  y ~ treatment | x1 + x2,
  data = dat,
  method = "dr-date",
  n_permutations = 100,
  seed = 1
)
result_f
#> 
#>    DR-DATE Test
#> 
#> Statistic: 0.0271132 
#> P-value:   0.0198 
#> N:         300 
#> Perms:     100 
#> Kernel Y:  rbf (bw = 1.557)
#> ESS:       133.7
```

## When to Use Which Test

| Test | Detects | Overlap Requirement | Best For |
|----|----|----|----|
| **DR-DATE** | Any distributional difference | Both sides | Population-level effects |
| **DR-DETT** | Distributional effect on treated | One-sided only | Imperfect overlap; policy questions about treated |
| **DML/TMLE** | Mean shifts only | Both sides | When only mean effects matter |

## References

- Fawkes, J., Hu, R., Evans, R. J., & Sejdinovic, D. (2024). Doubly
  robust kernel statistics for testing distributional treatment effects.
  *Transactions on Machine Learning Research*.

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
