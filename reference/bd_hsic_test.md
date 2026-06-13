# Backdoor-HSIC Test for Causal Association

Tests the do-null hypothesis H_0: p(y \| do(x)) = p\*(y) using a
kernel-based test with backdoor adjustment via density ratio estimation.
Detects causal associations including non-linear effects that standard
linear methods miss.

## Usage

``` r
bd_hsic_test(
  x,
  y,
  z,
  kernel_x = kernel_spec(),
  kernel_y = kernel_spec(),
  density_ratio = c("logistic", "ranger", "xgboost", "proxymix", "rulsif"),
  n_permutations = 500L,
  n_clusters = "auto",
  split_ratio = 0.5,
  alpha = 0.05,
  seed = NULL,
  verbose = FALSE,
  cluster_id = NULL,
  permutation = c("auto", "within_cluster", "naive"),
  min_ess_fraction = 0.1
)
```

## Arguments

- x:

  Numeric vector or matrix. Treatment variable.

- y:

  Numeric vector or matrix. Outcome variable.

- z:

  Numeric matrix, data.frame, or data.table. Confounders.

- kernel_x:

  Kernel specification for treatment space. Default is RBF with median
  heuristic.

- kernel_y:

  Kernel specification for outcome space. Default is RBF with median
  heuristic.

- density_ratio:

  Character. Method for density ratio estimation: `"logistic"`
  (default), `"ranger"`, `"xgboost"`, `"proxymix"`, or `"rulsif"`. The
  `"proxymix"` backend fits Gaussian-mixture proxies to the joint and
  product-of-marginals sample clouds via classical EM (Hoek & Elliott,
  2024), giving a parametric alternative to NCE-based classifiers;
  useful for multimodal densities or when classifier calibration is
  unreliable. Requires the `proxymix` package (`>= 0.3.0`).

- n_permutations:

  Integer. Number of permutations for the null distribution. Default is
  500.

- n_clusters:

  Integer or `"auto"`. Number of *propensity* clusters for valid
  permutation when `cluster_id = NULL`. Default is `"auto"`.

- split_ratio:

  Numeric in (0, 1). Proportion of data for training the density ratio
  estimator. Default is 0.5.

- alpha:

  Numeric. Significance level. Default is 0.05.

- seed:

  Integer or `NULL`. Random seed for reproducibility.

- verbose:

  Logical. Print progress. Default is `FALSE`.

- cluster_id:

  Optional vector of length `nrow(x)` identifying external clusters
  (e.g. site, season, paddock, farm). When supplied, the permutation
  null is built by within-cluster reshuffling of `y`, which preserves
  cluster-level effects in the null. Coerced to factor; the test split
  inherits the cluster assignment. The result then carries a per-cluster
  stratified bd-HSIC alongside the pooled statistic.

- permutation:

  Character. Permutation scheme:

  - `"auto"` (default) – when `cluster_id` is supplied, equivalent to
    `"within_cluster"`; otherwise falls back to k-means clustering on
    propensity weights (the original Hu/Sejdinovic/Evans scheme).

  - `"within_cluster"` – requires `cluster_id`. Permutes `y` indices
    only within clusters; preserves cluster-level effects.

  - `"naive"` – unrestricted permutation across all observations. Use
    only when independence within clusters is plausible (rarely true in
    ag-systems data).

- min_ess_fraction:

  Numeric in `(0, 1)` or `0` / non-finite to disable. ESS-floor
  reliability gate (added 0.0.0.9013): if the weighted-HSIC effective
  sample size is below `min_ess_fraction * n_test`, a warning is emitted
  and `result$ess_warning` is `TRUE`. The default `0.1` (10%) is a
  conservative floor; tighten it for studies with strict reliability
  requirements.

## Value

An object of class `"kernel_test_result"`. When `cluster_id` is
supplied, the result additionally carries:

- permutation_scheme:

  Character: which scheme was used.

- cluster_id:

  Integer cluster assignment on the test split.

- cluster_levels:

  Character cluster labels.

- per_cluster_statistic:

  Per-cluster weighted HSIC (stratified contributions); `NA` for
  clusters with `< 2` test observations.

## Details

The bd-HSIC test (Hu, Sejdinovic & Evans, 2024) tests whether treatment
X has a causal effect on outcome Y after adjusting for confounders Z via
the backdoor criterion.

The test works by:

1.  Estimating density ratios w(x, z) = p\*(x) / p(x\|z) to reweight
    observational samples to the interventional distribution.

2.  Computing a weighted HSIC statistic between X and Y.

3.  Obtaining p-values via permutation of Y within exchangeability
    clusters – propensity-similarity clusters by default, or external
    design clusters (site / season / paddock) when `cluster_id` is
    supplied.

**Hierarchical extension.** When the design is naturally clustered
(multi-site agricultural trials, paddock x season factorial designs,
patient x hospital data), supplying `cluster_id` activates
within-cluster permutation: indices of `y` are reshuffled only within
each cluster, preserving cluster-level effects in the null. This is the
safer default for clustered data; naive permutation across clusters can
inflate Type I error when cluster effects exist.

Unlike PDS or Double ML, bd-HSIC can detect non-linear causal effects
(e.g., U-shaped relationships) where the treatment affects higher
moments of the outcome but not necessarily the mean.

**Small samples.** The minimum sample size is `6` (the train/test split
must leave at least two test observations for the weighted HSIC and two
propensity clusters). Small-N field trials are supported, but
reliability is not guaranteed by size alone: the ESS-floor gate
(`min_ess_fraction`) will warn when a small or poorly-overlapping sample
yields a weighted statistic dominated by a handful of high-weight
points.

## Train/test split (0.0.0.9014)

The density-ratio estimator is now **fit on the train split and
predicted on the held-out test split** via
[`fit_density_ratio()`](https://max578.github.io/kernR/reference/fit_density_ratio.md) +
[`predict_density_ratio()`](https://max578.github.io/kernR/reference/predict_density_ratio.md).
The documented `split_ratio` is honoured end-to-end for all four
classifier / proxymix backends. The 0.0.0.9013 sample-split leak warning
is therefore retired. RuLSIF, the kernel-based closed-form backend,
still uses
[`estimate_rulsif()`](https://max578.github.io/kernR/reference/estimate_rulsif.md)
on the train/test split natively.

The fitted density-ratio model is preserved on
`result$density_ratio_fit` for callers that want backend diagnostics
(see
[`?fit_density_ratio`](https://max578.github.io/kernR/reference/fit_density_ratio.md)
Value; proxymix exposes BIC, AIC, log-likelihood, convergence per GMM).

## proxymix fit-quality gate

For `density_ratio = "proxymix"`, the backend surfaces a single
`fit_quality` verdict from its per-GMM convergence diagnostics. If a
mixture proxy fails to converge, the density-ratio weights are
unreliable; the test then emits a warning and sets
`result$density_ratio_warning = TRUE` (it is `FALSE` for a clean fit and
for all other backends). This mirrors the ESS-floor gate: an
untrustworthy verdict is flagged, never reported silently.

## References

Hu, R., Sejdinovic, D., & Evans, R. J. (2024). A kernel test for causal
association via noise contrastive backdoor adjustment. *JMLR*, 25(160),
1-56.

## See also

Other causal association tests:
[`hierarchical_test()`](https://max578.github.io/kernR/reference/hierarchical_test.md),
[`kernel_causal_test()`](https://max578.github.io/kernR/reference/kernel_causal_test.md),
[`taci_test()`](https://max578.github.io/kernR/reference/taci_test.md)

## Examples

``` r
set.seed(42)
n <- 300
z <- matrix(rnorm(n * 2), n, 2)
x <- z[, 1] + rnorm(n)
y <- 0.5 * x + z[, 2] + rnorm(n, sd = 0.5)

result <- bd_hsic_test(x, y, z, n_permutations = 200, seed = 1)
print(result)
#> 
#>    bd-HSIC Test
#> 
#> Statistic: 0.0157602 
#> P-value:   0.0896 
#> N:         150 
#> Perms:     200 
#> Kernel X:  rbf (bw =  1.32)
#> Kernel Y:  rbf (bw = 1.338)
#> ESS:       149.1 
#> 
```
