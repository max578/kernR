# Kernel-Based Downscaling: Vector and Distribution Inputs

## Three flavours of kernel downscaling

kernR ships three downscaling primitives, covering distinct
data-availability regimes:

1.  **[`kernel_downscale()`](https://max578.github.io/kernR/reference/kernel_downscale.md)**
    – vector-in / vector-out. Conditional mean embedding (CME)
    regression: each training observation is a single coarse-resolution
    vector; predict fine-resolution outputs at new coarse inputs.
    Standard kernel ridge regression in conditional-distribution form.
    Used when “coarse” and “fine” live on the same kind of object
    (e.g. one row per year).

2.  **[`dist_regression()`](https://max578.github.io/kernR/reference/dist_regression.md)**
    – bag-in / vector-out. Distribution regression: each training input
    is a *bag of points* (a sample from a distribution); predict a
    scalar or vector at the bag level. Used when the coarse covariate is
    itself a distribution rather than a single vector – e.g. a paddock’s
    soil-core measurements, an ensemble of climate-model realisations.

3.  **[`aggregate_downscale()`](https://max578.github.io/kernR/reference/aggregate_downscale.md)**
    – aggregate-likelihood inversion of a *known* aggregator. No paired
    (coarse, fine) training data — only the aggregate observation `y`
    and the operator `T` such that `Y = T(X) + eps`. A Gaussian-mixture
    prior on the latent `X` is inverted in closed form (linear-Gaussian
    `T`) or by importance sampling within each prior component
    (non-linear `T`). The Sejdinovic kernel-downsizing /
    aggregate-likelihood direction; companion to the proxymix
    `from_aggregate_likelihood()` Tier-2 stub.

The first two have closed-form solutions; the third has a closed form
for linear aggregators and a Monte-Carlo path for non-linear. All three
use orders of magnitude less data than deep-learning downscalers and
carry interpretable bandwidth / regularisation / component-count knobs.

## (1) `kernel_downscale()`: coarse climate → paddock yield

``` r

library(kernR)

# Training: 8 years of coarse climate -> paddock-level yield/biomass
n <- 80L
coarse <- matrix(
  c(stats::rnorm(n, mean = 18, sd = 2),    # mean monthly temperature
    stats::rnorm(n, mean = 450, sd = 80)), # cumulative rainfall (mm)
  ncol = 2L,
  dimnames = list(NULL, c("temp", "rainfall"))
)
truth <- function(z) {
  cbind(
    yield   = 0.02 * z[, "rainfall"] - 0.1 * (z[, "temp"] - 18)^2 +
              stats::rnorm(nrow(z), sd = 0.5),
    biomass = 0.03 * z[, "rainfall"] + 0.05 * z[, "temp"] +
              stats::rnorm(nrow(z), sd = 0.7)
  )
}
fine <- truth(coarse)
```

Predict at a held-out coarse grid:

``` r

new_coarse <- matrix(
  c(stats::rnorm(20L, 18, 2),
    stats::rnorm(20L, 450, 80)),
  ncol = 2L,
  dimnames = list(NULL, c("temp", "rainfall"))
)
fit <- kernel_downscale(coarse, fine, new_coarse)
fit
#> 
#>   Kernel Downscaling (CME)
#> 
#> Training pairs:   80 
#> Prediction points: 20 
#> Output dims:      2 
#> Kernel (coarse):  rbf (bw = 72.69)
#> Lambda (ridge):   0.0003162
```

``` r

head(fit$prediction)
#>           yield  biomass
#> [1,] 10.1869683 16.11029
#> [2,]  7.5849960 11.98351
#> [3,]  9.0823386 15.11778
#> [4,]  7.9127329 13.75706
#> [5,]  0.8226318  9.93738
#> [6,]  6.3918668 10.91197
```

The `lambda` was CV-selected via leave-one-out; the bandwidths follow
the median heuristic per kernel. Both can be overridden.

### When to reach for `kernel_downscale()`

- Both training and prediction have the same coarse-resolution *shape*:
  each row is one coarse vector.
- You want point predictions at the fine resolution, not a full
  conditional distribution.
- Sample size in the tens-to-thousands range; closed-form ridge scales
  as `O(n^3)` for the train-time solve but predicts in
  `O(n_new * n_train)`.

For very large training sets, plug
[`nystrom_factor()`](https://max578.github.io/kernR/reference/nystrom_factor.md)
outputs into a custom CME by writing `K_x = F F^T` and re-solving the
ridge in factor space – see
\[[`fit_cme()`](https://max578.github.io/kernR/reference/fit_cme.md)\]
for the operator-level access.

## (2) `dist_regression()`: distribution → scalar

The simplest illustrative case: predict the *mean* of a distribution
from the distribution itself, given a finite sample.

``` r

M <- 60L
mu_train <- stats::runif(M, -3, 3)
bags_train <- lapply(mu_train, function(m) {
  matrix(stats::rnorm(40L, mean = m, sd = 1), ncol = 1L)
})

fit_dr <- dist_regression(bags_train, y = mu_train, outer = "linear")
fit_dr
#> 
#>   Kernel Distribution Regression
#> 
#> Training bags:     60 
#> Bag sizes:        40-40 (median 40)
#> Point dim:         1 
#> Output dim:        1 
#> Inner kernel:      rbf (bw =  2.01)
#> Outer kernel:      linear
#> Ridge lambda:      0.0001
```

``` r

# Predict on fresh bags
M_new <- 10L
mu_new <- stats::runif(M_new, -3, 3)
bags_new <- lapply(mu_new, function(m) {
  matrix(stats::rnorm(40L, mean = m, sd = 1), ncol = 1L)
})
pred <- predict(fit_dr, bags_new)
data.frame(
  truth      = round(mu_new, 3L),
  prediction = round(pred,    3L),
  abs_error  = round(abs(pred - mu_new), 3L)
)
#>     truth prediction abs_error
#> 1  -2.258     -2.039     0.219
#> 2   0.256      0.323     0.067
#> 3  -0.717     -0.658     0.059
#> 4  -0.159     -0.066     0.094
#> 5   2.827      2.855     0.028
#> 6  -2.241     -2.323     0.082
#> 7   0.338      0.684     0.346
#> 8  -1.333     -1.112     0.221
#> 9   0.221      0.325     0.104
#> 10  1.103      1.140     0.036
```

### Bag-level vector targets

`y` can also be a matrix when each bag carries multiple scalar targets:

``` r

y_mat <- cbind(mean = mu_train, mean_sq = mu_train^2)
fit_mv <- dist_regression(bags_train, y = y_mat, outer = "linear")
predict(fit_mv, bags_new[1:3])
#>            mean   mean_sq
#> [1,] -2.0388518 4.4290140
#> [2,]  0.3229703 0.2034361
#> [3,] -0.6579623 0.6030948
```

### Outer kernel choice

- **`outer = "linear"`** – inner-product of empirical embeddings. Most
  appropriate when the target is well captured by the first moment of
  the bag (e.g. predicting bag mean).
- **`outer = "rbf"`** – Gaussian on embedding-space distance. Picks up
  higher-order distributional features. Use when the target depends on
  more than the bag mean (variance, skew, tails). Bandwidth follows the
  median heuristic on the embedding-space pairwise distances.

``` r

# Variance-targeting example: zero-mean bags with varying SDs
s_train <- stats::runif(M, 0.5, 2.5)
bags_var <- lapply(s_train, function(s) {
  matrix(stats::rnorm(60L, mean = 0, sd = s), ncol = 1L)
})

# Linear outer: insensitive to bag SD (means are all near zero)
fit_lin <- dist_regression(bags_var, y = s_train, outer = "linear")
fit_lin
#> 
#>   Kernel Distribution Regression
#> 
#> Training bags:     60 
#> Bag sizes:        60-60 (median 60)
#> Point dim:         1 
#> Output dim:        1 
#> Inner kernel:      rbf (bw = 1.518)
#> Outer kernel:      linear
#> Ridge lambda:      3.162e-05

# RBF outer: embedding-space distance differs even with equal means
fit_rbf <- dist_regression(bags_var, y = s_train, outer = "rbf")
fit_rbf
#> 
#>   Kernel Distribution Regression
#> 
#> Training bags:     60 
#> Bag sizes:        60-60 (median 60)
#> Point dim:         1 
#> Output dim:        1 
#> Inner kernel:      rbf (bw = 1.582)
#> Outer kernel:      rbf (bw = 0.1897)
#> Ridge lambda:      0.001
```

### Variable bag sizes

Bag sizes may differ across the design (e.g. paddocks with different
numbers of soil cores). The double-sum embedding handles this naturally
via [`mean()`](https://rdrr.io/r/base/mean.html) over the inner Gram:

``` r

mu_v <- stats::runif(20L, -2, 2)
bags_v <- lapply(mu_v, function(m) {
  n_i <- sample(20:80, 1L)
  matrix(stats::rnorm(n_i, mean = m), ncol = 1L)
})
fit_v <- dist_regression(bags_v, y = mu_v, outer = "linear")
range(vapply(fit_v$bags_train, nrow, integer(1L)))
#> [1] 23 80
```

### When to reach for `dist_regression()`

- Each training row is itself a *collection of measurements* rather than
  a single vector (paddock soil cores, ensemble realisations, spectral
  samples, repeated experimental measurements).
- The target is a property of the underlying distribution, not of a
  single point.
- Sample size in bags: 10s-100s of points per bag; number of bags
  10s-low 1000s. Cost is `O(M^2 \cdot \bar{m}^2)` for the inner Gram
  where `M` is bag count, `\bar{m}` mean bag size.

## (3) `aggregate_downscale()`: aggregate-likelihood inversion

Often the downscaling problem looks **nothing** like supervised
regression: no paired (coarse, fine) data is available. Instead, the
observation is `y = T(x) + eps` where `T` is a *known* operator (spatial
average, temporal average, satellite footprint convolution, a non-linear
sensor model), `eps` is observation noise, and a parametric prior on the
fine-scale latent `x` is available — for example, fitted via
[`proxymix::fit_proxymix()`](https://rdrr.io/pkg/proxymix/man/fit_proxymix.html)
on historical fine-scale data, or specified directly. The job is to
recover the posterior `p(x | y)`.

[`aggregate_downscale()`](https://max578.github.io/kernR/reference/aggregate_downscale.md)
dispatches on the aggregator’s class:

- **Linear matrix `A`** (`y = A x + eps`): closed-form per-component
  Kalman update plus reweighted mixture weights. Cost
  `O(N (dim_x^3 + dim_y^3))`.
- **Function `T(x)`** (non-linear): importance sampling within each
  prior component (`n_samples_per_component` draws, default 200).
  Per-component effective sample size is reported and an ESS-floor
  warning fires when IS collapses. Cost `O(N M (dim_x + dim_y))` per
  call.

Linear example — spatial averaging of two adjacent paddocks against a
two-cluster prior:

``` r

A <- matrix(c(0.5, 0.5), nrow = 1L)   # spatial average
prior <- list(
  means       = list(c(0, 0), c(2, 2)),
  covariances = list(diag(2L), diag(2L)),
  weights     = c(0.5, 0.5)
)
# Observe an aggregate that points to the second cluster
fit_lin <- aggregate_downscale(y = 1.8, aggregator = A,
                               latent_prior = prior, sigma_y = 0.15)
fit_lin
#> Aggregate-likelihood downscaling
#>   method:        linear_closed_form (linear)
#>   components:    2
#>   sigma_y:       0.15
#>   posterior mean:  1.8,  1.8 
#>   posterior wts: 0.045, 0.955
```

The posterior mixture weights shift toward the cluster consistent with
the observation (here, the `(2, 2)` cluster). The per-component
posterior is exact (Kalman update) and the cross-component law-of-
total-covariance is honest.

Non-linear example — a sinusoidal aggregator (sensor non-linearity):

``` r

agg_fn <- function(x) matrix(sin(rowSums(x)), ncol = 1L)
fit_nl <- aggregate_downscale(y = 0.5, aggregator = agg_fn,
                              latent_prior = prior, sigma_y = 0.1,
                              n_samples_per_component = 400L,
                              seed = 1L)
fit_nl$ess_per_component
#> [1] 53.49922 42.53140
fit_nl$posterior_mean
#> [1] 0.9601427 0.8983002
```

The per-component ESS is the reliability gate — if it collapses, the
function warns and the posterior should be treated as exploratory. Draw
posterior samples for downstream uncertainty propagation:

``` r

samples <- posterior_sample_aggregate(fit_nl, n = 500L, seed = 2L)
head(samples)
#>            [,1]       [,2]
#> [1,]  0.6648714  1.3170075
#> [2,] -0.2855020  0.9779909
#> [3,]  1.4959683  2.6373047
#> [4,] -0.1101655  0.1078616
#> [5,]  2.2575797  1.6941538
#> [6,]  0.3705832 -1.0080668
```

### When to reach for `aggregate_downscale()`

- Only **aggregate** observations are available (no paired fine-scale
  ground truth).
- The aggregation operator `T` is known (spatial averaging, temporal
  averaging, sensor convolution, mass-balance constraint).
- A parametric latent prior is available — e.g. from
  [`proxymix::fit_proxymix()`](https://rdrr.io/pkg/proxymix/man/fit_proxymix.html)
  on historical fine-scale data, or user-supplied from a domain model.
- You want a closed-form posterior when `T` is linear, or a clean
  importance-sampling pass when it is not — without committing to
  2000. 

Pair with
[`proxymix::fit_proxymix()`](https://rdrr.io/pkg/proxymix/man/fit_proxymix.html)
(the canonical prior fitter) and PESTO’s manifest contract (the
cross-package handoff) for the full APSIM → posterior-prior →
downscaled-fine-scale chain.

## Comparison

| Aspect | [`kernel_downscale()`](https://max578.github.io/kernR/reference/kernel_downscale.md) | [`dist_regression()`](https://max578.github.io/kernR/reference/dist_regression.md) | [`aggregate_downscale()`](https://max578.github.io/kernR/reference/aggregate_downscale.md) |
|----|----|----|----|
| Training data | Paired (coarse, fine) | Bags + scalar/vector targets | None paired; aggregate `y` only |
| Coarse-fine map | Learned (regression) | Learned (regression) | **Known** (`y = T(x) + eps`) |
| Latent prior | Implicit (empirical) | Implicit | Explicit (parametric GMM) |
| Closed form | Yes | Yes | Yes (linear `T`) / IS (non-linear) |
| Multi-output | Yes | Yes | Yes |
| Posterior uncertainty | Point + ridge CI | Point + ridge CI | Full mixture posterior |

## References

- Park, J., Muandet, K., Fukumizu, K., & Sejdinovic, D. (2013). *Kernel
  embeddings of conditional distributions.* IEEE Signal Processing
  Magazine, 30(4), 98-111.
- Szabó, Z., Sriperumbudur, B. K., Póczos, B., & Gretton, A. (2016).
  Learning theory for distribution regression. *Journal of Machine
  Learning Research*, 17(152), 1-40.
- Muandet, K., Fukumizu, K., Sriperumbudur, B., & Schölkopf, B. (2017).
  *Kernel mean embedding of distributions: A review and beyond.*
  Foundations and Trends in Machine Learning, 10(1-2), 1-141.
