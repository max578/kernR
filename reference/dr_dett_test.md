# Doubly Robust Distributional Effect on the Treated Test (DR-DETT)

Tests whether the distribution of the treated potential outcome Y(1)
differs from the control potential outcome Y(0) among the treated
subpopulation. Requires only one-sided overlap.

## Usage

``` r
dr_dett_test(
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
effective sample size of the control reconstruction weights and
`ess_warning` records whether the reliability floor was hit.

## Details

DR-DETT is analogous to DR-DATE but focuses on the **effect on the
treated** (ETT) rather than the average treatment effect. The treated
counterfactual Y(1) \| T = 1 is observed directly, so only the *control*
arm needs an outcome model. The control counterfactual Y(0) \| T = 1 is
reconstructed by augmented inverse-probability weighting, reweighting
controls by the treatment odds \\e(x) / (1 - e(x))\\ so that their
covariate distribution matches the treated. With `outcome_model = "krr"`
the control conditional mean embedding \\\hat m_0\\ supplies the doubly
robust augmentation; `outcome_model = "zero"` returns the
inverse-probability-weighted statistic. It requires only one-sided
overlap: P(T = 1 \| X) bounded away from 0 (not necessarily from 1), so
it applies where positivity fails for the control group.

## References

Fawkes, J., Hu, R., Evans, R. J., & Sejdinovic, D. (2024). Doubly robust
kernel statistics for testing distributional treatment effects.
*Transactions on Machine Learning Research*.

## See also

Other distributional treatment effects:
[`dr_date_scenario()`](https://max578.github.io/kernR/reference/dr_date_scenario.md),
[`dr_date_test()`](https://max578.github.io/kernR/reference/dr_date_test.md)

## Examples

``` r
set.seed(42)
n <- 300
x <- matrix(rnorm(n * 2), n, 2)
t <- rbinom(n, 1, plogis(0.5 * x[, 1]))
y <- t * rnorm(n, sd = 2) + (1 - t) * rnorm(n, sd = 1) + x[, 1]

result <- dr_dett_test(y, t, x, n_permutations = 200, seed = 1)
print(result)
#> 
#>    DR-DETT Test
#> 
#> Statistic: 0.0342116 
#> P-value:   0.0149 
#> N:         300 
#> Perms:     200 
#> Kernel Y:  rbf (bw = 1.716)
#> ESS:       102.2 
#> 
```
