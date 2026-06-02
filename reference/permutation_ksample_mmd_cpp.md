# Permutation k-sample MMD: summed pairwise unbiased MMD^2 under joint relabel

Given the pooled (N x N) kernel matrix of K stacked groups and their
sizes, draws n_perm joint relabelings of the pooled sample into the
original group sizes and returns, for each, the summed pairwise unbiased
MMD^2 statistic \\\sum\_{a \< b} \mathrm{MMD}^2_u(a, b)\\. The single
shared relabeling per replicate (not independent per-pair permutation)
is what makes this a valid k-sample null. Relabeling uses r_randperm, so
callers honour set.seed().

## Usage

``` r
permutation_ksample_mmd_cpp(K_pool, sizes, n_perm)
```

## Arguments

- K_pool:

  (N x N) kernel matrix of the row-stacked groups.

- sizes:

  Integer vector of the K group sizes (summing to N).

- n_perm:

  Number of permutations.

## Value

Vector of n_perm summed-pairwise MMD^2 values under the null.
