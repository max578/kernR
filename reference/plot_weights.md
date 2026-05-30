# Plot Weight Diagnostics

Plots the distribution of importance weights with effective sample size
annotation.

## Usage

``` r
plot_weights(weights, main = "Weight Distribution")
```

## Arguments

- weights:

  Numeric vector of importance weights.

- main:

  Title. Default is "Weight Distribution".

## Value

Invisibly returns `weights`.

## Examples

``` r
set.seed(1L)
weights <- rgamma(200L, shape = 2, rate = 2)
plot_weights(weights)

```
