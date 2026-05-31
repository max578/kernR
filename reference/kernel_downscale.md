# Kernel-Based Statistical Downscaling

Predicts fine-resolution outputs at new coarse-resolution inputs via
conditional mean embedding (CME) regression in an RKHS. Given paired
coarse-fine training data `(coarse, fine)`, fits the operator \\E\[Y
\mid X = x\]\\ in closed form via kernel ridge regression and returns
predictions at `new_coarse`. This is the
Park-Muandet-Fukumizu-Sejdinovic conditional-mean-embedding scheme,
specialised to the regression form needed for spatial / temporal
downscaling.

## Usage

``` r
kernel_downscale(
  coarse,
  fine,
  new_coarse,
  kernel_coarse = kernel_spec(),
  kernel_fine = kernel_spec(),
  lambda = "cv",
  return_weights = FALSE
)
```

## Arguments

- coarse:

  Numeric matrix `n x d_coarse` of training coarse-resolution inputs.
  Vectors are coerced via
  [`as.matrix()`](https://rdrr.io/r/base/matrix.html).

- fine:

  Numeric matrix `n x d_fine` of training fine-resolution outputs.
  Multivariate outputs (`d_fine > 1`) are supported and predicted
  jointly.

- new_coarse:

  Numeric matrix `n_new x d_coarse` of coarse inputs at which to predict
  the fine outputs.

- kernel_coarse, kernel_fine:

  [`kernel_spec()`](https://max578.github.io/kernR/reference/kernel_spec.md)
  for the coarse and fine spaces. Defaults to RBF with median heuristic.
  `kernel_fine` is used only for the bandwidth-CV step; predictions are
  returned in the original `fine` units.

- lambda:

  Ridge regularisation parameter for the CME ridge regression. If `"cv"`
  (default), selected by leave-one-out cross-validation over
  `10^seq(-6, 1, length.out = 15)`.

- return_weights:

  Logical. If `TRUE`, the result carries the `n_new x n_train` weight
  matrix used to combine training fine values. Default `FALSE` (saves
  memory for large designs).

## Value

An object of class `"kernel_downscale"` with components:

- prediction:

  `n_new x d_fine` matrix of predicted fine outputs at `new_coarse`.

- n_train:

  Number of training pairs used.

- n_new:

  Number of prediction points.

- lambda:

  Regularisation used (CV-selected when `lambda = "cv"`).

- kernel_coarse, kernel_fine:

  Resolved kernel specs.

- weights:

  Optional `n_new x n_train` weight matrix.

- call:

  The matched call.

## Details

Typical ag-systems use: coarse climate-grid inputs (e.g. monthly
temperature, rainfall on a 25 km grid) -\> fine-resolution outputs
(paddock yield, biomass) at the same time index. Train on years where
both coarse and fine are observed; predict fine outputs at new coarse
inputs.

Compared to a linear regression baseline, the kernel approach captures
non-linear coarse-fine relationships without specifying the functional
form. Compared to deep-learning downscalers, it has a closed-form
solution, uses orders-of-magnitude less data, and carries an
interpretable kernel-bandwidth degrees-of-freedom knob.

## References

Park, J., Muandet, K., Fukumizu, K., & Sejdinovic, D. (2013). *Kernel
embeddings of conditional distributions.* IEEE Signal Processing
Magazine.

Muandet, K., Fukumizu, K., Sriperumbudur, B., & Scholkopf, B. (2017).
*Kernel mean embedding of distributions: A review and beyond.*
Foundations and Trends in Machine Learning, 10(1-2).

## See also

[`fit_cme()`](https://max578.github.io/kernR/reference/fit_cme.md),
[`dist_regression()`](https://max578.github.io/kernR/reference/dist_regression.md)
(for downscaling when each coarse "input" is a bag of points rather than
a single vector).

Other downscaling and embeddings:
[`aggregate_downscale()`](https://max578.github.io/kernR/reference/aggregate_downscale.md),
[`dist_regression()`](https://max578.github.io/kernR/reference/dist_regression.md),
[`fit_cme()`](https://max578.github.io/kernR/reference/fit_cme.md),
[`posterior_sample_aggregate()`](https://max578.github.io/kernR/reference/posterior_sample_aggregate.md)

## Examples

``` r
# \donttest{
set.seed(1)
n <- 80L
coarse <- matrix(stats::rnorm(n * 2L), n, 2L,
                 dimnames = list(NULL, c("temp", "rainfall")))
fine <- cbind(
  yield   = 2 * coarse[, "rainfall"] -
            0.5 * coarse[, "temp"]^2 +
            stats::rnorm(n, sd = 0.2),
  biomass = coarse[, "temp"] + coarse[, "rainfall"]^2 +
            stats::rnorm(n, sd = 0.2)
)
# Predict at a held-out grid
new_coarse <- matrix(stats::rnorm(20L * 2L), 20L, 2L,
                     dimnames = list(NULL, c("temp", "rainfall")))
fit <- kernel_downscale(coarse, fine, new_coarse)
print(fit)
#> 
#>   Kernel Downscaling (CME)
#> 
#> Training pairs:   80 
#> Prediction points: 20 
#> Output dims:      2 
#> Kernel (coarse):  rbf (bw =  1.49)
#> Lambda (ridge):   0.0001 
#> 
head(fit$prediction)
#>            yield    biomass
#> [1,]  0.65693954  2.0952777
#> [2,] -3.00126072  2.6748749
#> [3,]  0.08983194  0.9410209
#> [4,] -2.53991905 -1.3344578
#> [5,] -4.53189548  4.4139403
#> [6,] -2.81292492  2.4399332
# }
```
