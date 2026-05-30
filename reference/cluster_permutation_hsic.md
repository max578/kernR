# Cluster-Based Permutation for bd-HSIC

Permutes Y indices within clusters of similar conditional densities
p(x\|z), ensuring valid exchangeability under the null.

## Usage

``` r
cluster_permutation_hsic(Kx, Ky, weights, clusters, n_permutations)
```

## Arguments

- Kx:

  n x n kernel matrix for X.

- Ky:

  n x n kernel matrix for Y.

- weights:

  Density ratio weights.

- clusters:

  Integer vector of cluster assignments.

- n_permutations:

  Number of permutations.

## Value

Vector of permuted weighted HSIC statistics.
