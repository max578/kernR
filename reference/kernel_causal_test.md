# Unified Kernel Causal Test

A convenience wrapper that dispatches to the appropriate test function
based on the `method` argument.

## Usage

``` r
kernel_causal_test(
  formula,
  data,
  method = c("dr-date", "dr-dett", "bd-hsic"),
  ...
)
```

## Arguments

- formula:

  A formula of the form `y ~ treatment | confounders`.

- data:

  A data.frame or data.table.

- method:

  Character. Test method: `"dr-date"` (default), `"dr-dett"`, or
  `"bd-hsic"`.

- ...:

  Additional arguments passed to the specific test function.

## Value

An object of class `"kernel_test_result"`.

## See also

Other causal association tests:
[`bd_hsic_test()`](https://max578.github.io/kernR/reference/bd_hsic_test.md),
[`hierarchical_test()`](https://max578.github.io/kernR/reference/hierarchical_test.md),
[`taci_test()`](https://max578.github.io/kernR/reference/taci_test.md)

## Examples

``` r
set.seed(42)
n <- 200
dat <- data.frame(
  y = rnorm(n),
  treatment = rbinom(n, 1, 0.5),
  x1 = rnorm(n),
  x2 = rnorm(n)
)
dat$y <- dat$y + 0.5 * dat$treatment + 0.3 * dat$x1

result <- kernel_causal_test(y ~ treatment | x1 + x2,
  data = dat, method = "dr-date",
  n_permutations = 100, seed = 1
)
print(result)
#> 
#>    DR-DATE Test
#> 
#> Statistic: 0.0801612 
#> P-value:   0.0099 
#> N:         200 
#> Perms:     100 
#> Kernel Y:  rbf (bw = 1.007)
#> ESS:       82.2 
#> 
```
