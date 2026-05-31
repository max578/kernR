# Doubly Robust Distributional Average Treatment Effect Test (DR-DATE)

Tests whether the distributions of potential outcomes Y(1) and Y(0)
differ using a doubly robust kernel MMD statistic. Detects
distributional effects (variance, shape) that mean-based tests miss.

## Usage

``` r
dr_date_test(
  y,
  treatment,
  covariates,
  kernel_y = kernel_spec(),
  propensity_model = c("logistic", "ranger", "xgboost"),
  outcome_model = c("krr", "zero"),
  cross_fit = TRUE,
  n_folds = 5L,
  n_permutations = 500L,
  n_bins = 10L,
  regularisation = "cv",
  min_ess_fraction = 0.1,
  alpha = 0.05,
  seed = NULL,
  verbose = FALSE
)
```

## Arguments

- y:

  Numeric vector or matrix. Outcome variable.

- treatment:

  Binary vector (0/1). Treatment indicator.

- covariates:

  Numeric matrix, data.frame, or data.table. Confounders.

- kernel_y:

  Kernel specification for outcome space. Default is RBF.

- propensity_model:

  Character. Method for propensity estimation: `"logistic"` (default),
  `"ranger"`, or `"xgboost"`.

- outcome_model:

  Character. `"krr"` (kernel ridge regression, default) fits a
  conditional mean embedding for each arm and forms the doubly robust
  statistic; `"zero"` drops the outcome model and returns the
  inverse-probability-weighted (singly robust) statistic.

- cross_fit:

  Logical. If `TRUE` (default), both nuisances – the propensity score
  and the conditional mean embedding – are estimated by `n_folds`-fold
  cross-fitting and evaluated out-of-fold, as the doubly robust theory
  requires under flexible nuisance estimators. If `FALSE`, both are fit
  in-sample (faster, but the test can be anti-conservative).

- n_folds:

  Integer. Number of cross-fitting folds. Default is 5.

- n_permutations:

  Integer. Number of permutations. Default is 500.

- n_bins:

  Integer. Propensity score bins for permutation. Default is 10.

- regularisation:

  Numeric or `"cv"`. Ridge parameter for the CME. Default is `"cv"`.

- min_ess_fraction:

  Numeric in (0, 1). If the effective sample size of either arm's
  inverse-probability weights falls below this fraction of `n`, a
  reliability [`warning()`](https://rdrr.io/r/base/warning.html) is
  emitted. Default 0.1.

- alpha:

  Numeric. Significance level. Default is 0.05.

- seed:

  Integer or `NULL`. Random seed. Permutations are drawn through R's
  RNG, so a fixed `seed` makes the test fully reproducible.

- verbose:

  Logical. Print progress. Default is `FALSE`.

## Value

An object of class `"kernel_test_result"`. The `ess` element holds the
smaller of the two per-arm effective sample sizes and `ess_warning`
records whether the reliability floor was hit.

## Details

The DR-DATE test (Fawkes, Hu, Evans & Sejdinovic, 2024) constructs
doubly robust (augmented inverse-probability-weighted) estimators for
the counterfactual mean embeddings of Y(1) and Y(0) in a reproducing
kernel Hilbert space. For arm \\a\\ the augmented embedding is
\$\$\hat\mu_a = \frac{1}{n}\sum_i \tilde w\_{a,i} \bigl(k(y_i,\cdot) -
\hat m_a(x_i)\bigr) + \hat m_a(x_i),\$\$ where \\\hat m_a\\ is the
conditional mean embedding fitted on arm \\a\\ and \\\tilde w\_{a,i}\\
are stabilised inverse-probability weights. The statistic is
\\\\\hat\mu_1 - \hat\mu_0\\^2\\ in the RKHS. Setting
`outcome_model = "zero"` sets \\\hat m_a \equiv 0\\ and recovers the
inverse-probability-weighted statistic.

**Double robustness**: the test is consistent if *either* the propensity
model *or* the outcome (CME) model is correctly specified. Cross-fitting
(`cross_fit = TRUE`) makes this hold under flexible machine-learning
nuisances by removing own-observation overfitting bias (Chernozhukov et
al., 2018).

**Permutation null**: the reference distribution permutes treatment
labels *within* propensity-score bins, holding the fitted nuisances
fixed. This is valid under within-bin exchangeability of treatment given
the (binned) propensity score; calibration degrades as the bins coarsen
relative to the propensity variation inside them.

**Key advantage**: unlike DML or TMLE which test only for mean shifts,
DR-DATE detects *any* distributional difference including changes in
variance, skewness, or shape.

## References

Fawkes, J., Hu, R., Evans, R. J., & Sejdinovic, D. (2024). Doubly robust
kernel statistics for testing distributional treatment effects.
*Transactions on Machine Learning Research*.

Chernozhukov, V., Chetverikov, D., Demirer, M., Duflo, E., Hansen, C.,
Newey, W., & Robins, J. (2018). Double/debiased machine learning for
treatment and structural parameters. *The Econometrics Journal*, 21(1),
C1-C68.

## See also

Other distributional treatment effects:
[`dr_date_scenario()`](https://max578.github.io/kernR/reference/dr_date_scenario.md),
[`dr_dett_test()`](https://max578.github.io/kernR/reference/dr_dett_test.md)

## Examples

``` r
set.seed(42)
n <- 300
x <- matrix(rnorm(n * 2), n, 2)
logit_p <- 0.5 * x[, 1]
t <- rbinom(n, 1, plogis(logit_p))
y <- t * 1.0 + x[, 1] + rnorm(n, sd = 0.5)

result <- dr_date_test(y, t, x, n_permutations = 200, seed = 1)
print(result)
#> 
#>    DR-DATE Test
#> 
#> Statistic: 0.176911 
#> P-value:   0.0050 
#> N:         300 
#> Perms:     200 
#> Kernel Y:  rbf (bw = 1.249)
#> ESS:       133.7 
#> 
```
