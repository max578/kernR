# Estimate Conditional Mean Embedding via Kernel Ridge Regression

Estimates the conditional mean embedding \\\mu\_{Y\|X=x}\\ in the RKHS
using kernel ridge regression (Park, Muandet, Fukumizu & Sejdinovic,
2013; Muandet et al., 2017).

## Usage

``` r
fit_cme(
  x,
  y,
  kernel_x = kernel_spec(),
  kernel_y = kernel_spec(),
  lambda = "cv"
)
```

## Arguments

- x:

  Numeric matrix of conditioning variables (n x d_x).

- y:

  Numeric matrix of target variables (n x d_y).

- kernel_x:

  Kernel specification for `x`.

- kernel_y:

  Kernel specification for `y`.

- lambda:

  Ridge regularisation parameter. If `"cv"`, selected by leave-one-out
  cross-validation. Default is `"cv"`.

## Value

A list of class `"cme_fit"` with components:

- W:

  Operator matrix `(K_x + n lambda I)^{-1}` (n x n).

- Ky:

  Kernel matrix of `y`.

- x_train:

  Training `x` data.

- kernel_x, kernel_y:

  Resolved kernel specifications.

- lambda:

  Regularisation parameter used.

## Details

This is the lower-level building block used by
[`kernel_downscale()`](https://max578.github.io/kernR/reference/kernel_downscale.md).
Most users should call
[`kernel_downscale()`](https://max578.github.io/kernR/reference/kernel_downscale.md);
use `fit_cme()` directly when you need access to the trained operator
(the weight matrix `W`) for custom downstream computations – e.g.
constructing a kernel Bayes' rule update, plugging into a manuscript
figure pipeline, or composing with other RKHS operators.

## References

Park, J., Muandet, K., Fukumizu, K., & Sejdinovic, D. (2013). *Kernel
embeddings of conditional distributions: A unified kernel framework for
nonparametric inference in graphical models.* IEEE Signal Processing
Magazine.

Muandet, K., Fukumizu, K., Sriperumbudur, B., & Scholkopf, B. (2017).
*Kernel mean embedding of distributions: A review and beyond.*
Foundations and Trends in Machine Learning, 10(1-2).

## See also

[`predict.cme_fit()`](https://max578.github.io/kernR/reference/predict.cme_fit.md),
[`kernel_downscale()`](https://max578.github.io/kernR/reference/kernel_downscale.md)

## Examples

``` r
set.seed(1L)
x <- matrix(rnorm(60L), ncol = 2L)
y <- matrix(x[, 1L] + rnorm(30L, sd = 0.2), ncol = 1L)
fit <- fit_cme(x, y, lambda = 1e-2)
dim(fit$W)
#> [1] 30 30
```
