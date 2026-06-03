# Factor-space wild bootstrap of the kernel Stein discrepancy null

Low-rank counterpart of ksd_wild_bootstrap_cpp. Given a Nystrom factor F
(n x m) with \\F F^\top \approx U\\ approximating the Stein-kernel
matrix, and the factor trace \\\mathrm{tr}(F F^\top)\\, draws n_boot
wild- bootstrap replicates of the degenerate U-statistic null via
Rademacher multipliers: each replicate is \\(\lVert F^\top w \rVert^2 -
\mathrm{tr}) / (n (n - 1))\\, computed in O(n m) rather than O(n^2).
Multipliers are drawn through R's RNG, so callers honour set.seed().

## Usage

``` r
ksd_wild_bootstrap_factor_cpp(F, tr, n_boot)
```

## Arguments

- F:

  Numeric n x m Nystrom factor of the Stein-kernel matrix.

- tr:

  Trace of F F^t (sum of squared factor entries).

- n_boot:

  Number of bootstrap replicates.

## Value

Vector of n_boot bootstrap KSD statistics.
