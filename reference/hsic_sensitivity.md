# HSIC-Based Distributional Sensitivity Index

Computes HSIC-based sensitivity indices (Da Veiga, 2015) for each input
parameter against each output. By default returns first-order indices
only; `total_order = TRUE` adds the complementary total-order
decomposition.

## Usage

``` r
hsic_sensitivity(
  theta,
  y,
  total_order = FALSE,
  p_value = TRUE,
  n_permutations = 500L,
  total_order_ci = FALSE,
  n_bootstrap = 200L,
  ci_level = 0.95,
  total_order_test = c("none", "cond_perm"),
  n_clusters_cp = "auto",
  p_adjust = c("BH", "holm", "hochberg", "bonferroni", "BY", "fdr", "none"),
  kernel_theta = kernel_spec(),
  kernel_y = kernel_spec(),
  seed = NULL,
  total_order_p_value = NULL
)
```

## Arguments

- theta:

  Numeric matrix `n x p` of input draws (one row per simulator run, one
  column per parameter). Vectors are coerced via
  [`as.matrix()`](https://rdrr.io/r/base/matrix.html).

- y:

  Numeric matrix `n x q` of simulator outputs. Vectors are coerced via
  [`as.matrix()`](https://rdrr.io/r/base/matrix.html).

- total_order:

  Logical. If `TRUE`, additionally compute the total-order
  HSIC-Sensitivity Index via Da Veiga's complement formulation. Requires
  `p >= 2`. Default `FALSE` (backwards- compatible).

- p_value:

  Logical. If `TRUE` (default), compute permutation p-values and
  BH-adjusted p-values for *first-order* indices. Set `FALSE` to skip
  the permutation null. Has no effect on total-order indices (see
  Details).

- n_permutations:

  Integer. Permutations per HSIC test. Default `500`.

- total_order_ci:

  Logical. If `TRUE`, compute pair-bootstrap percentile CIs for the
  total-order indices. Requires `total_order = TRUE`. Default `FALSE`
  (backwards-compatible). Replaces the misleading 0.0.0.9012
  `total_order_p_value` argument (removed in 0.0.0.9013; see Details).

- n_bootstrap:

  Integer. Pair-bootstrap resamples used for the total-order CI. Default
  `200`. Only consulted when `total_order_ci = TRUE`.

- ci_level:

  Numeric in `(0, 1)`. Two-sided percentile CI level for the total-order
  bootstrap. Default `0.95`. Only consulted when
  `total_order_ci = TRUE`.

- total_order_test:

  Character. `"none"` (default; backwards- compatible) or `"cond_perm"`
  (since 0.0.0.9014). The latter activates a conditional-permutation
  test for `H_0: X_j _||_ Y | X_{~j}` and populates
  `p_value_total_order` and `p_value_total_order_adjusted`. Requires
  `total_order = TRUE`. Different from (and replaces) the retracted
  0.0.0.9012 pair-bootstrap method.

- n_clusters_cp:

  Integer or `"auto"`. Number of conditioning bins for the
  conditional-permutation test (k-means on `X_{~j}`). `"auto"` chooses
  `min(floor(n / 5), 20)`. Only consulted when
  `total_order_test = "cond_perm"`.

- p_adjust:

  Character. Across-grid p-value adjustment method passed to
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html). Default
  `"BH"`.

- kernel_theta:

  A
  [`kernel_spec()`](https://max578.github.io/kernR/reference/kernel_spec.md)
  for parameter columns. Default RBF with median heuristic. Used for
  both first-order (per-column) and total-order (per `X_{~j}` subset).

- kernel_y:

  A
  [`kernel_spec()`](https://max578.github.io/kernR/reference/kernel_spec.md)
  for output columns. Default RBF with per-column median heuristic.

- seed:

  Integer or `NULL`. Random seed for reproducibility.

- total_order_p_value:

  Defunct as of 0.0.0.9013. Passing any non-`NULL` value errors with a
  pointer to `total_order_test` and `total_order_ci`.

## Value

An object of class `"hsic_sensitivity"` with components:

- index:

  `p x q` matrix of first-order HSIC-Sensitivity Indices in `[0, 1]`.

- statistic:

  `p x q` matrix of raw first-order HSIC statistics.

- p_value:

  `p x q` matrix of raw permutation p-values (first-order), or `NULL` if
  `p_value = FALSE`.

- p_value_adjusted:

  `p x q` matrix of adjusted p-values (first-order), or `NULL`.

- index_total_order:

  `p x q` matrix of total-order indices, or `NULL` when
  `total_order = FALSE`.

- statistic_total_order:

  `p x q` matrix of HSIC(`X_{~j}`, `Y_k`) statistics, or `NULL`.

- total_order:

  Logical flag: was total-order requested?

- total_index:

  Length-`p` vector: per-parameter maximum *first-order* index across
  outputs (used for ranking; kept under this historical name for
  backwards compatibility).

- rank:

  Parameter indices ordered by descending `total_index`.

- n, n_permutations, param_names, output_names, p_adjust:

  Metadata.

- call:

  The matched call.

## Details

**First-order index** (always computed): \$\$S^{HSIC}\_j =
\frac{HSIC(X_j, Y)}{\sqrt{HSIC(X_j, X_j) \cdot HSIC(Y, Y)}}\$\$ is the
direct contribution of `X_j` to `Y` – analogous to (but not equal to)
the Sobol first-order index.

**Total-order index** (when `total_order = TRUE`): \$\$T^{HSIC}\_j = 1 -
\frac{HSIC(X\_{\sim j}, Y)}{\sqrt{HSIC(X\_{\sim j}, X\_{\sim j}) \cdot
HSIC(Y, Y)}}\$\$ where `X_{~j}` is the parameter design with column `j`
removed. By construction the difference `T_j - S_j` captures the
contribution of `X_j` *through interactions* with other parameters. For
purely additive models `T_j = S_j`; in the presence of interaction
`T_j > S_j`.

Unlike variance-based Sobol indices, both versions of the
HSIC-Sensitivity Index capture non-linear and distributional effects: a
parameter that affects the variance, skewness, or tail of `Y` without
shifting its mean is invisible to Sobol but visible to HSIC. The
normalisation bounds each index to `[0, 1]`.

Optional permutation p-values are computed for first-order indices
(Benjamini-Hochberg-adjusted across the grid by default).

**Total-order uncertainty quantification (`total_order_ci = TRUE`, since
0.0.0.9013).** A pair-bootstrap percentile CI on the index `T_j` itself:
resample `(theta, y)` pairs with replacement `n_bootstrap` times,
recompute `T_j` on each resample, report a `ci_level` (default 95%)
two-sided percentile interval. This is uncertainty quantification, not a
hypothesis test.

**Total-order conditional-permutation significance test
(`total_order_test = "cond_perm"`, since 0.0.0.9014).** Tests
`H_0: X_j _||_ Y | X_{~j}` – under the null, `X_j` adds nothing beyond
what is already in `X_{~j}` and `T_j` is concentrated at zero. The null
is generated by **conditional permutation**: k-means-cluster the design
points by `X_{~j}` similarity, then within each cluster permute `Y`;
recompute `T_j` on each permuted design. p-value
`= (1 + #{T_perm >= T_obs}) / (1 + B)`. The `p_adjust` method applies
grid-wide.

The 0.0.0.9012 `total_order_p_value` pair-bootstrap implementation was
**not** this – it sampled the empirical joint, not a null-of-no-effect,
and was retracted in 0.0.0.9013. The 0.0.0.9014 conditional-permutation
test repopulates `p_value_total_order` with a properly-calibrated value;
the `total_order_test` flag on the result distinguishes the new mode
from the retracted one.

**Power caveat (empirical calibration 2026-05-16).** Repeated-seed
calibration (`orchestra_calibration_20260516.R`, B = 100, N = 80)
confirms the cond_perm test is null-calibrated (per-parameter type-I
rates 0.01-0.04 at nominal alpha = 0.05) but **conservative on additive
designs** (rejection rate ~3% on a strong-additive design where the
first-order test rejects above 90%). The reason is structural: on
additive `Y = sum_j f_j(X_j)`, the conditional permutation of `Y` within
`X_{~j}` bins preserves all of `X_{~j}`'s contribution and randomises
only the `X_j` component (which is independent of `X_{~j}`), so the
permuted `HSIC(X_{~j}, Y_perm)` is statistically close to the observed
`HSIC(X_{~j}, Y)`. Use cond_perm for **interaction detection**; use the
first-order permutation `p_value` (always computed when
`p_value = TRUE`) for **additive contributions**. The two are
complementary, not interchangeable.

Total-order CIs cost `B * p * q` kernel re-evaluations on resampled
designs; set `n_bootstrap` accordingly. Use the CI for honest
uncertainty bars on the index magnitude; do not interpret it as a
significance verdict.

## References

Da Veiga, S. (2015). Global sensitivity analysis with dependence
measures. *Journal of Statistical Computation and Simulation*, 85(7),
1283-1305.

Gretton, A., Bousquet, O., Smola, A., & Scholkopf, B. (2005). Measuring
statistical dependence with Hilbert-Schmidt norms. *ALT*, 63-77.

## See also

[`hsic_identifiability()`](https://max578.github.io/kernR/reference/hsic_identifiability.md),
[`hsic_test()`](https://max578.github.io/kernR/reference/hsic_test.md)

Other sensitivity and identifiability:
[`hsic_identifiability()`](https://max578.github.io/kernR/reference/hsic_identifiability.md),
[`lhs_design()`](https://max578.github.io/kernR/reference/lhs_design.md)

## Examples

``` r
# \donttest{
set.seed(1)
n <- 80L
theta <- matrix(stats::runif(n * 3L), nrow = n,
                dimnames = list(NULL, c("active", "weak", "inert")))
y <- 2 * theta[, "active"] +
     0.5 * theta[, "weak"] +
     stats::rnorm(n, sd = 0.1)
fit <- hsic_sensitivity(theta, y, total_order = TRUE,
                        n_permutations = 199L, seed = 1L)
fit
#> 
#>   HSIC-Sensitivity Indices
#> 
#> Parameters:   3 
#> Outputs:      1 
#> N:            80 
#> Permutations: 199 
#> P-adjust:     BH 
#> Total-order: yes
#> 
#> Per-parameter ranking (descending S, max across outputs):
#>   S = first-order index   T = total-order index   interaction = T - S
#> 
#>  parameter S_first_max T_total_max interaction min_p_first
#>     active       0.879       0.956       0.077      0.0150
#>       weak       0.055       0.419       0.364      0.1200
#>      inert       0.014       0.328       0.314      0.7100
#> 
# }
```
