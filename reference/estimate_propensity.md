# Estimate Propensity Scores

Estimates P(T = 1 \| X) using the specified model, with built-in
cross-fitting support.

## Usage

``` r
estimate_propensity(
  treatment,
  covariates,
  method = c("logistic", "ranger", "xgboost"),
  cross_fit = TRUE,
  n_folds = 5L,
  trim = 0.01,
  seed = NULL
)
```

## Arguments

- treatment:

  Binary vector (0/1). Treatment indicator.

- covariates:

  Numeric matrix or data.frame. Confounders.

- method:

  Character. `"logistic"` (default), `"ranger"`, or `"xgboost"`.

- cross_fit:

  Logical. If `TRUE`, uses 5-fold cross-fitting to produce out-of-sample
  propensity estimates. Default is `TRUE`.

- n_folds:

  Integer. Number of cross-fitting folds. Default is 5.

- trim:

  Numeric. Trim extreme propensity scores to `[trim, 1-trim]`. Default
  is 0.01.

- seed:

  Integer or `NULL`. Random seed for the cross-fitting fold assignment,
  so a fixed `seed` makes cross-fitted scores reproducible.

## Value

A list of class `"propensity_fit"` with components:

- scores:

  Estimated propensity scores P(T=1\|X).

- method:

  Method used.

- trim:

  Trimming threshold applied.

- n_trimmed:

  Number of scores that were trimmed.

## Examples

``` r
set.seed(42)
n <- 300
x <- matrix(rnorm(n * 3), n, 3)
logit_p <- 0.5 * x[, 1] - 0.3 * x[, 2]
t <- rbinom(n, 1, plogis(logit_p))
ps <- estimate_propensity(t, x)
summary(ps$scores)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.2537  0.4428  0.5121  0.5171  0.5989  0.7847 
```
