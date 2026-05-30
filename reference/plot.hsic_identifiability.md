# Plot an HSIC Identifiability Scan

Bar plot of per-parameter maximum HSIC across outputs, ordered by
magnitude. Bars for identifiable parameters are coloured; non-
identifiable parameters are shown in grey.

## Usage

``` r
# S3 method for class 'hsic_identifiability'
plot(x, col_yes = "#0072B2", col_no = "grey70", ...)
```

## Arguments

- x:

  An `hsic_identifiability` object.

- col_yes, col_no:

  Bar colours for identifiable / non-identifiable parameters.

- ...:

  Additional arguments passed to
  [`graphics::barplot()`](https://rdrr.io/r/graphics/barplot.html).

## Value

Invisibly returns `x`. Side effect: produces a base R plot.

## Examples

``` r
# \donttest{
set.seed(1)
n <- 50
theta <- matrix(stats::runif(n * 3), nrow = n,
                dimnames = list(NULL, c("active", "active2", "inert")))
y <- theta[, 1] + theta[, 2]^2 + stats::rnorm(n, sd = 0.1)
fit <- hsic_identifiability(theta, y, n_permutations = 199, seed = 1)
plot(fit)

# }
```
