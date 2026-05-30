# Out-of-Fold Conditional-Mean-Embedding Coefficient Matrix

Builds the n x n matrix `C` whose row `i` holds the coefficients of the
fitted conditional mean embedding \\\hat m(x_i)\\ over the n outcome
embeddings \\k(y_l, \cdot)\\. The CME smoother weights are a function of
the conditioning (covariate) kernel only, so they are valid coordinates
over the global outcome basis. Columns are non-zero only for the arm's
training rows and, under cross-fitting, only for rows outside `i`'s
fold. An arm (or a fold's training arm) with fewer than 10 units leaves
a zero block (IPW-only for those rows); when `warn = TRUE` this is
surfaced via a single
[`warning()`](https://rdrr.io/r/base/warning.html).

## Usage

``` r
.cme_coef_matrix(
  y,
  covariates,
  arm_idx,
  n,
  lambda,
  cross_fit,
  n_folds,
  fold_id,
  kernel_y,
  warn = TRUE
)
```
