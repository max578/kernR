# Aggregate-Likelihood Downscaling

Inverts a known aggregation operator `Y = T(X) + eps` to recover the
fine-scale latent `X` from coarse / aggregate observations `Y`, using a
Gaussian-mixture prior on the latent. Implements the
aggregate-likelihood / kernel-downsizing framework (Sejdinovic et al.)
as a kernR-side method consuming an optional
[`proxymix::fit_proxymix()`](https://rdrr.io/pkg/proxymix/man/fit_proxymix.html)
latent prior.

## Usage

``` r
aggregate_downscale(
  y,
  aggregator,
  latent_prior,
  sigma_y = 0.1,
  n_samples_per_component = 200L,
  min_ess_fraction = 0.1,
  seed = NULL
)
```

## Arguments

- y:

  Numeric vector of length `dim_y`, or `1 x dim_y` matrix. The observed
  aggregate.

- aggregator:

  Either a numeric `dim_y x dim_x` matrix (treated as a linear
  aggregator and dispatched to the closed-form path) or a function
  `function(x) -> y_matrix` that maps an `n x dim_x` matrix of latent
  samples to an `n x dim_y` matrix of aggregates (non-linear IS path).

- latent_prior:

  Either (a) a list with elements `means` (`N x dim_x` matrix or list of
  `dim_x`-vectors), `covariances` (list of `N` `dim_x x dim_x`
  matrices), `weights` (length-`N` numeric, summing to 1) – or (b) a
  [`proxymix::gmm_fit`](https://rdrr.io/pkg/proxymix/man/gmm_fit.html)
  (any object exposing `@means`, `@covariances`, `@weights` slots).

- sigma_y:

  Numeric. Observation noise standard deviation (scalar;
  `eps ~ N(0, sigma_y^2 I)`). Default `0.1`.

- n_samples_per_component:

  Integer. Importance-sampling sample count per prior component
  (non-linear path only). Default `200L`.

- min_ess_fraction:

  Numeric in `(0, 1]`. ESS-floor reliability gate for the IS path: when
  per-component ESS drops below
  `min_ess_fraction * n_samples_per_component`, a warning is emitted.
  Default `0.1`. Set `0` to disable.

- seed:

  Integer or `NULL`. Random seed (non-linear path).

## Value

An object of class `"aggregate_downscale"` with components:

- `posterior_mean` – length-`dim_x` numeric, `E[X | y]`.

- `posterior_cov` – `dim_x x dim_x` matrix, `Cov[X | y]`
  (law-of-total-covariance over mixture components).

- `posterior_weights` – length-`N` numeric, posterior mixture weights
  (sum to 1).

- `posterior_components_means` – list of `N` length-`dim_x` posterior
  component means.

- `posterior_components_covariances` – list of `N` posterior component
  covariances.

- `aggregator_type` – `"linear"` or `"nonlinear"`.

- `method` – `"linear_closed_form"` or `"nonlinear_is"`.

- `ess_per_component` – length-`N` per-component ESS (IS path only; `NA`
  for closed form).

- `ess_warning` – `TRUE` if any per-component ESS fell below the floor
  (IS path only).

- `n_components`, `sigma_y`, `n_samples_per_component`, `call`.

## Details

This is the **third** downscaling method in kernR, structurally
different from the existing pair:

- [`kernel_downscale()`](https://max578.github.io/kernR/reference/kernel_downscale.md)
  (CME, Park-Muandet-Fukumizu-Sejdinovic 2013) – paired (coarse, fine)
  training data; supervised regression.

- [`dist_regression()`](https://max578.github.io/kernR/reference/dist_regression.md)
  (Szabo-Sriperumbudur-Poczos-Gretton 2016) –
  distribution-to-distribution mapping from a bag-of-points design.

- `aggregate_downscale()` (this function) – only aggregate
  observations + known aggregator + parametric latent prior. Used when
  no paired training data exists (only coarse-grid observations) but the
  aggregation operator is known (spatial averaging, temporal averaging,
  linear or non-linear projection).

Two computational paths, selected on `aggregator`'s class:

- **Linear-Gaussian closed form** (when `aggregator` is a matrix `A`):
  each prior component's posterior is a Kalman update;
  `K_k = Sigma_k A^T (A Sigma_k A^T + sigma_y^2 I)^{-1}`;
  `mu_k|y = mu_k + K_k (y - A mu_k)`, `Sigma_k|y = (I - K_k A) Sigma_k`;
  posterior mixture weights reweight by per-component evidence
  `N(y | A mu_k, A Sigma_k A^T + sigma_y^2 I)`.

- **Non-linear importance sampling** (when `aggregator` is a function):
  draw `n_samples_per_component` samples from each prior component,
  evaluate `T(.)`, weight by Gaussian likelihood
  `N(y | T(x), sigma_y^2 I)`, recover posterior moments + reweighted
  mixture weights from the importance-weighted samples. Reports
  per-component effective sample size; warns when below a stated floor.

The aggregate-likelihood / GMM-proxy direction is shared with the
companion proxymix Tier-2 stub
[`proxymix::from_aggregate_likelihood()`](https://rdrr.io/pkg/proxymix/man/from_aggregate_likelihood.html),
which targets the same problem from the prior-fitting side; this
function targets it from the consumption side (inversion given a fitted
prior).

## See also

[`kernel_downscale()`](https://max578.github.io/kernR/reference/kernel_downscale.md),
[`dist_regression()`](https://max578.github.io/kernR/reference/dist_regression.md).

## Examples

``` r
set.seed(1L)
# Linear-Gaussian: spatial averaging of two adjacent cells.
A <- matrix(c(0.5, 0.5), nrow = 1L)
prior <- list(
  means = list(c(0, 0), c(2, 2)),
  covariances = list(diag(2L), diag(2L)),
  weights = c(0.5, 0.5)
)
fit <- aggregate_downscale(y = 1.0, aggregator = A,
                           latent_prior = prior, sigma_y = 0.2)
fit$posterior_mean
#> [1] 1 1

# Non-linear: aggregator is sin of the sum.
agg_fn <- function(x) matrix(sin(rowSums(x)), ncol = 1L)
fit2 <- aggregate_downscale(y = 0.5, aggregator = agg_fn,
                            latent_prior = prior, sigma_y = 0.1,
                            n_samples_per_component = 300L, seed = 1L)
#> Warning: aggregate_downscale(): per-component IS effective sample size fell below 10% of n_samples_per_component for component(s) 2 (ESS = 27.9). Increase `n_samples_per_component`, sharpen the prior, or consider a coarser aggregator if posterior is locally narrow.
fit2$posterior_mean
#> [1] 0.8567869 0.8729094
```
