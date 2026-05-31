# Latin-Hypercube Design Over Bounded Parameters

Generates a Latin-hypercube sample of size `n` over the parameter bounds
in `bounds`. Each column is a random permutation of stratified uniform
draws (one per equal-width bin in `(0, 1]`), then scaled to the supplied
parameter range. The result is reproducible when `seed` is supplied.

## Usage

``` r
lhs_design(n, bounds, seed = NULL)
```

## Arguments

- n:

  Integer. Number of design points (rows). Must be `>= 2`.

- bounds:

  Two-column numeric matrix or data.frame of `[lower, upper]` bounds,
  with one row per parameter. If named, row names propagate to the
  column names of the returned design.

- seed:

  Integer or `NULL`. Random seed for reproducibility.

## Value

A numeric matrix of dimension `n x nrow(bounds)`. Column names are
inherited from `rownames(bounds)` when available, otherwise `theta1`,
`theta2`, ....

## Details

This is a lightweight helper aimed at pre-PESTO screening: produce a
design matrix to feed an APSIM (or any) simulator, then pass the
resulting input/output pairs to
[`hsic_identifiability()`](https://max578.github.io/kernR/reference/hsic_identifiability.md)
to flag unidentifiable parameters before ensemble-smoother calibration.

## References

McKay, M. D., Beckman, R. J., & Conover, W. J. (1979). A comparison of
three methods for selecting values of input variables in the analysis of
output from a computer code. *Technometrics*, 21(2), 239-245.

## See also

[`hsic_identifiability()`](https://max578.github.io/kernR/reference/hsic_identifiability.md)

Other sensitivity and identifiability:
[`hsic_identifiability()`](https://max578.github.io/kernR/reference/hsic_identifiability.md),
[`hsic_sensitivity()`](https://max578.github.io/kernR/reference/hsic_sensitivity.md)

## Examples

``` r
bounds <- rbind(
  slope     = c(0.1, 2.0),
  intercept = c(-1, 1),
  noise_sd  = c(0.05, 0.5)
)
design <- lhs_design(50, bounds, seed = 1)
head(design)
#>          slope  intercept   noise_sd
#> [1,] 0.2388260  0.5559651 0.18341284
#> [2,] 1.5693166  0.1229079 0.06067908
#> [3,] 0.1198987  0.3486707 0.05838398
#> [4,] 1.3580965 -0.6236229 0.07339595
#> [5,] 0.9411551 -0.9644144 0.12972970
#> [6,] 1.7191804  0.6863797 0.39026021
```
