# Sample from the posterior of an aggregate-downscale fit

Draws `n` samples from the per-component posterior mixture. Each draw
picks a component by `posterior_weights`, then samples from
`N(posterior_components_means[[k]], posterior_components_covariances[[k]])`.

## Usage

``` r
posterior_sample_aggregate(object, n = 1000L, seed = NULL)
```

## Arguments

- object:

  An `aggregate_downscale` fit.

- n:

  Integer. Number of posterior samples. Default `1000L`.

- seed:

  Integer or `NULL`. Random seed.

## Value

An `n x dim_x` numeric matrix.

## Examples

``` r
set.seed(1L)
A <- matrix(c(0.5, 0.5), nrow = 1L)
prior <- list(
  means = list(c(0, 0), c(2, 2)),
  covariances = list(diag(2L), diag(2L)),
  weights = c(0.5, 0.5)
)
fit <- aggregate_downscale(y = 1.0, aggregator = A,
                           latent_prior = prior, sigma_y = 0.2)
draws <- posterior_sample_aggregate(fit, n = 500L, seed = 1L)
colMeans(draws)
#> [1] 1.0163065 0.9701081
```
