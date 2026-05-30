# Plot a Kernel Test Result

Plots the permutation null distribution with the observed statistic.

## Usage

``` r
# S3 method for class 'kernel_test_result'
plot(x, ...)
```

## Arguments

- x:

  A `kernel_test_result` object.

- ...:

  Additional arguments (currently ignored).

## Value

Invisibly returns `x`. Side effect: produces a base R plot.

## Examples

``` r
set.seed(42)
x_data <- rnorm(100)
y_data <- x_data + rnorm(100, sd = 0.5)
res <- hsic_test(x_data, y_data)
plot(res)

```
