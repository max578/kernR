# HSIC Independence Test

Tests whether two variables are independent using the Hilbert-Schmidt
Independence Criterion (HSIC). Uses a permutation test for inference.

## Usage

``` r
hsic_test(
  x,
  y,
  kernel_x = kernel_spec(),
  kernel_y = kernel_spec(),
  n_permutations = 500L,
  alpha = 0.05,
  seed = NULL
)
```

## Arguments

- x:

  Numeric vector, matrix, or data.frame. First variable.

- y:

  Numeric vector, matrix, or data.frame. Second variable.

- kernel_x:

  Kernel specification for `x`. Default is RBF with median heuristic.

- kernel_y:

  Kernel specification for `y`. Default is RBF with median heuristic.

- n_permutations:

  Integer. Number of permutations for the null distribution. Default is
  500.

- alpha:

  Numeric. Significance level. Default is 0.05.

- seed:

  Integer or `NULL`. Random seed for reproducibility.

## Value

An object of class `"kernel_test_result"` with components:

- statistic:

  The observed HSIC test statistic.

- p_value:

  Permutation p-value.

- method:

  `"HSIC"`.

- n:

  Sample size.

- n_permutations:

  Number of permutations used.

- null_distribution:

  Vector of permuted HSIC values.

- kernel_x, kernel_y:

  Kernel specifications used (with resolved bandwidths).

- call:

  The matched call.

## References

Gretton, A., Fukumizu, K., Teo, C. H., Song, L., Scholkopf, B., & Smola,
A. J. (2008). A kernel statistical test of independence. *NeurIPS*, 20.

## See also

Other independence and two-sample tests:
[`mmd_test()`](https://max578.github.io/kernR/reference/mmd_test.md)

## Examples

``` r
set.seed(42)
n <- 200
x <- rnorm(n)

# Dependent case
y_dep <- x^2 + rnorm(n, sd = 0.5)
result <- hsic_test(x, y_dep)
print(result)
#> 
#>    HSIC Test
#> 
#> Statistic: 0.0220599 
#> P-value:   0.0020 
#> N:         200 
#> Perms:     500 
#> Kernel X:  rbf (bw = 0.9284)
#> Kernel Y:  rbf (bw = 0.9302)
#> 

# Independent case
y_ind <- rnorm(n)
result <- hsic_test(x, y_ind)
print(result)
#> 
#>    HSIC Test
#> 
#> Statistic: 0.000759494 
#> P-value:   0.5908 
#> N:         200 
#> Perms:     500 
#> Kernel X:  rbf (bw = 0.9284)
#> Kernel Y:  rbf (bw =  0.91)
#> 
```
