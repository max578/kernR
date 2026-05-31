# Kernel Distribution Regression

Regression where each input is a *bag* of points (a sample from a
distribution) rather than a single feature vector. Each bag is
implicitly mapped to its empirical mean embedding in the RKHS of an
inner kernel; the outer (between-bag) kernel acts on those embeddings,
and kernel ridge regression predicts a scalar or multivariate output.
This is the Szabó-Sriperumbudur-Póczos-Gretton (2016) "learning theory
for distribution regression" setup; kernel-mean-embedding background is
from Muandet et al. (2017).

## Usage

``` r
dist_regression(
  bags,
  y,
  inner_kernel = kernel_spec(),
  outer = c("linear", "rbf"),
  outer_bandwidth = "median",
  lambda = "cv"
)
```

## Arguments

- bags:

  A list of length `M` of numeric matrices. Each element is a bag
  (`n_i x d` points). All bags must have the same number of columns; bag
  sizes `n_i` may differ.

- y:

  Numeric vector of length `M`, or numeric matrix `M x d_y`, of training
  targets.

- inner_kernel:

  A
  [`kernel_spec()`](https://max578.github.io/kernR/reference/kernel_spec.md)
  applied between points within and across bags. Median bandwidth
  heuristic resolved on the pooled training points. Default RBF with
  median heuristic.

- outer:

  Character. Outer-kernel form: `"linear"` (default) or `"rbf"`.

- outer_bandwidth:

  Outer-kernel bandwidth (RBF only). `"median"` (default) resolves to
  the median embedding-space pairwise distance on the training Gram; a
  positive numeric overrides.

- lambda:

  Ridge regularisation for kernel ridge regression. If `"cv"` (default),
  selected by leave-one-out CV over `10^seq(-6, 1, length.out = 15)`.

## Value

An object of class `"dist_regression"` with components:

- alpha:

  Ridge weights (length `M` or `M x d_y`).

- G_train:

  Training bag-inner-mean Gram `(M x M)`.

- K_outer_train:

  Outer kernel matrix `(M x M)` actually used.

- bags_train:

  Bags used (kept for prediction).

- y_train:

  Targets used.

- inner_kernel:

  Resolved inner kernel.

- outer, outer_bandwidth:

  Outer kernel choice and resolved bandwidth (`NA` for linear outer).

- lambda:

  Ridge parameter used.

- M, d_y, call:

  Metadata.

## Details

**Outer kernels supported.**

- `"linear"`: \\K(P_i, P_j) = \langle \hat{\mu}\_{P_i}, \hat{\mu}\_{P_j}
  \rangle = \frac{1}{n_i n_j} \sum\_{k,l} k(x_i^{(k)}, x_j^{(l)})\\.

- `"rbf"`: \\K(P_i, P_j) = \exp(-\\ \hat{\mu}\_{P_i} - \hat{\mu}\_{P_j}
  \\^2 / (2\sigma^2))\\, where the embedding-space distance is recovered
  from the inner Gram via \\\\\hat\mu_i - \hat\mu_j\\^2 = G\_{ii} - 2
  G\_{ij} + G\_{jj}\\.

Typical ag-systems use: each paddock contributes a *bag* of soil-core
measurements (variable depth, multiple cores per paddock); the
regression predicts paddock-level yield from the distributional shape of
the soil profile. Distinct from
[`kernel_downscale()`](https://max578.github.io/kernR/reference/kernel_downscale.md)
in that the input is itself a distribution, not a fixed-length vector.

## References

Szabó, Z., Sriperumbudur, B. K., Póczos, B., & Gretton, A. (2016).
Learning theory for distribution regression. *Journal of Machine
Learning Research*, 17(152), 1-40.

Muandet, K., Fukumizu, K., Sriperumbudur, B., & Scholkopf, B. (2017).
*Kernel mean embedding of distributions: A review and beyond.*
Foundations and Trends in Machine Learning, 10(1-2).

## See also

[`predict.dist_regression()`](https://max578.github.io/kernR/reference/predict.dist_regression.md),
[`kernel_downscale()`](https://max578.github.io/kernR/reference/kernel_downscale.md)

Other downscaling and embeddings:
[`aggregate_downscale()`](https://max578.github.io/kernR/reference/aggregate_downscale.md),
[`fit_cme()`](https://max578.github.io/kernR/reference/fit_cme.md),
[`kernel_downscale()`](https://max578.github.io/kernR/reference/kernel_downscale.md),
[`posterior_sample_aggregate()`](https://max578.github.io/kernR/reference/posterior_sample_aggregate.md)

## Examples

``` r
# \donttest{
set.seed(1)
# 30 bags, each a sample from N(mu_i, 1); predict mu_i
M <- 30L
mu <- stats::rnorm(M)
bags <- lapply(mu, function(m)
  matrix(stats::rnorm(40, mean = m), ncol = 1L))
fit <- dist_regression(bags, y = mu, outer = "linear")
fit
#> 
#>   Kernel Distribution Regression
#> 
#> Training bags:     30 
#> Bag sizes:        40-40 (median 40)
#> Point dim:         1 
#> Output dim:        1 
#> Inner kernel:      rbf (bw = 1.325)
#> Outer kernel:      linear
#> Ridge lambda:      3.162e-05 
#> 

# Predict at new bags
new_mu <- stats::rnorm(5L)
new_bags <- lapply(new_mu, function(m)
  matrix(stats::rnorm(40, mean = m), ncol = 1L))
predict(fit, new_bags)
#> [1] -0.1377160  1.0718395  0.1172893 -0.5587325 -0.6755410
# }
```
