# HSIC-Based Identifiability Diagnostic

Pre-PESTO (or pre-IES) screening: for each parameter `theta[, j]` and
each output `y[, k]`, computes an HSIC permutation test of independence
and flags parameters with no detectable association to any output as
unidentifiable. Useful for trimming the parameter space before
ensemble-smoother calibration of mechanistic ag-system models such as
APSIM.

## Usage

``` r
hsic_identifiability(
  theta,
  y,
  alpha = 0.05,
  p_adjust = c("BH", "holm", "hochberg", "bonferroni", "BY", "fdr", "none"),
  n_permutations = 500L,
  kernel_theta = kernel_spec(),
  kernel_y = kernel_spec(),
  seed = NULL
)
```

## Arguments

- theta:

  Numeric matrix `n x p` of parameter design points (one row per
  simulator run, one column per parameter). Vectors are coerced via
  [`as.matrix()`](https://rdrr.io/r/base/matrix.html).

- y:

  Numeric matrix `n x q` of simulator outputs. Vectors are coerced via
  [`as.matrix()`](https://rdrr.io/r/base/matrix.html).

- alpha:

  Numeric in `(0, 1)`. Identifiability threshold on the adjusted minimum
  p-value. Default `0.05`.

- p_adjust:

  Character. Across-grid p-value adjustment method passed to
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html). Default
  `"BH"`. Use `"none"` to disable.

- n_permutations:

  Integer. Permutations per HSIC test. Default 500.

- kernel_theta:

  A
  [`kernel_spec()`](https://max578.github.io/kernR/reference/kernel_spec.md)
  for parameter columns. Default RBF with per-column median heuristic.

- kernel_y:

  A
  [`kernel_spec()`](https://max578.github.io/kernR/reference/kernel_spec.md)
  for output columns. Default RBF with per-column median heuristic.

- seed:

  Integer or `NULL`. Random seed for reproducibility.

## Value

An object of class `"hsic_identifiability"` with components:

- statistic:

  `p x q` matrix of HSIC statistics.

- p_value:

  `p x q` matrix of raw permutation p-values.

- p_value_adjusted:

  `p x q` matrix of adjusted p-values (same as `p_value` when
  `p_adjust = "none"`).

- max_statistic:

  Length-`p` vector: per-parameter maximum HSIC across outputs.

- min_p_value:

  Length-`p` vector: per-parameter minimum adjusted p-value across
  outputs.

- identifiable:

  Length-`p` logical: `min_p_value <= alpha`.

- rank:

  Parameter indices ordered by descending `max_statistic`.

- alpha, p_adjust, n, n_permutations:

  Inputs / metadata.

- param_names, output_names:

  Character vectors.

- call:

  The matched call.

## Details

Kernel matrices are computed once per parameter and once per output, so
the total cost is `O((p + q) n^2)` (kernel construction) plus
`O(p q n_permutations n^2)` (permutation null), where `p` is the number
of parameters, `q` the number of outputs and `n` the design size.

A parameter is **identifiable** at level `alpha` when its smallest
(optionally adjusted) p-value across outputs satisfies `min_p <= alpha`.
Across-grid p-value adjustment defaults to Benjamini-Hochberg, which is
the natural FDR control for screening applications.

## References

Gretton, A., Fukumizu, K., Teo, C. H., Song, L., Scholkopf, B., & Smola,
A. J. (2008). A kernel statistical test of independence. *NeurIPS*, 20.

Da Veiga, S. (2015). Global sensitivity analysis with dependence
measures. *Journal of Statistical Computation and Simulation*, 85(7),
1283-1305.

## See also

[`lhs_design()`](https://max578.github.io/kernR/reference/lhs_design.md),
[`hsic_test()`](https://max578.github.io/kernR/reference/hsic_test.md)

Other sensitivity and identifiability:
[`hsic_sensitivity()`](https://max578.github.io/kernR/reference/hsic_sensitivity.md),
[`lhs_design()`](https://max578.github.io/kernR/reference/lhs_design.md)

## Examples

``` r
# \donttest{
set.seed(1)
n <- 60
# 3 active parameters + 1 inert
theta <- matrix(stats::runif(n * 4), nrow = n,
                dimnames = list(NULL, paste0("p", 1:4)))
y1 <- theta[, 1] + 0.5 * theta[, 2]^2 + stats::rnorm(n, sd = 0.1)
y2 <- sin(2 * pi * theta[, 3]) + stats::rnorm(n, sd = 0.1)
y <- cbind(yield = y1, biomass = y2)
fit <- hsic_identifiability(theta, y, n_permutations = 199, seed = 1)
print(fit)
#> 
#>   HSIC Identifiability Scan
#> 
#> Parameters:   4 
#> Outputs:      2 
#> N:            60 
#> Permutations: 199 
#> Alpha:        0.05 
#> P-adjust:     BH 
#> 
#> Identifiable (3): p1, p2, p3
#> Not identifiable (1): p4
#> 
#> Per-parameter ranking (descending max HSIC):
#>  parameter max_HSIC  min_p identifiable
#>         p3  0.07222 0.0200            *
#>         p1  0.05759 0.0200            *
#>         p2 0.009363 0.0400            *
#>         p4 0.001365 0.9543             
#> 
#>   (* = identifiable at alpha = 0.05 )
#> 
# }
```
