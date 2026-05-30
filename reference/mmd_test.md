# MMD Two-Sample Test

Tests whether two samples come from the same distribution using the
Maximum Mean Discrepancy (MMD). Uses a permutation test for inference.

## Usage

``` r
mmd_test(
  x,
  y,
  kernel = kernel_spec(),
  n_permutations = 500L,
  alpha = 0.05,
  seed = NULL
)
```

## Arguments

- x:

  Numeric vector, matrix, or data.frame. First sample.

- y:

  Numeric vector, matrix, or data.frame. Second sample.

- kernel:

  Kernel specification. Default is RBF with median heuristic.

- n_permutations:

  Integer. Number of permutations. Default is 500.

- alpha:

  Numeric. Significance level. Default is 0.05.

- seed:

  Integer or `NULL`. Random seed for reproducibility.

## Value

An object of class `"kernel_test_result"` with components:

- statistic:

  The observed MMD^2 statistic (unbiased).

- p_value:

  Permutation p-value.

- method:

  `"MMD"`.

- n:

  Total sample size (n_x + n_y).

- n_permutations:

  Number of permutations used.

- null_distribution:

  Vector of permuted MMD^2 values.

- kernel_x:

  Kernel specification used.

- call:

  The matched call.

## References

Gretton, A., Borgwardt, K. M., Rasch, M. J., Scholkopf, B., & Smola, A.
(2012). A kernel two-sample test. *JMLR*, 13, 723-773.

## Examples

``` r
set.seed(42)

# Same distribution
x <- matrix(rnorm(200), 100, 2)
y <- matrix(rnorm(200), 100, 2)
result <- mmd_test(x, y)
print(result)
#> 
#>    MMD Test
#> 
#> Statistic: -0.00292249 
#> P-value:   0.6327 
#> N:         200 
#> Perms:     500 
#> Kernel X:  rbf (bw = 1.601)
#> 

# Different distributions
y_shifted <- matrix(rnorm(200, mean = 1), 100, 2)
result <- mmd_test(x, y_shifted)
print(result)
#> 
#>    MMD Test
#> 
#> Statistic: 0.269115 
#> P-value:   0.0020 
#> N:         200 
#> Perms:     500 
#> Kernel X:  rbf (bw = 1.838)
#> 
```
