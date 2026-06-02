# Wild bootstrap of the kernel Stein discrepancy null

Given the Stein-kernel matrix H, draws n_boot wild-bootstrap replicates
of the degenerate U-statistic null via independent Rademacher
multipliers \\W_i \in \\-1, +1\\\\: each replicate is \\(1 / (n (n -
1))) \sum\_{i \ne j} W_i W_j H\_{ij}\\. The diagonal is excluded to
match the unbiased U-statistic. Multipliers are drawn through R's RNG
(R::unif_rand within Rcpp's RNGScope), so callers honour set.seed().

## Usage

``` r
ksd_wild_bootstrap_cpp(H, n_boot)
```

## Arguments

- H:

  Symmetric n x n Stein-kernel matrix.

- n_boot:

  Number of bootstrap replicates.

## Value

Vector of n_boot bootstrap KSD statistics.
