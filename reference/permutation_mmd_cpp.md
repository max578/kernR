# Permutation MMD: compute MMD^2 for many permutations of pooled data

Given the full (n+m) x (n+m) kernel matrix of the pooled sample, permute
the assignment into two groups and compute MMD^2.

## Usage

``` r
permutation_mmd_cpp(K_pool, n, m, n_perm)
```

## Arguments

- K_pool:

  (n+m) x (n+m) kernel matrix.

- n:

  Size of first sample.

- m:

  Size of second sample.

- n_perm:

  Number of permutations.

## Value

Vector of n_perm MMD^2 values under permutation.
