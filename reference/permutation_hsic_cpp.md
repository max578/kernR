# Permutation HSIC: compute HSIC for many Y permutations

Efficiently computes HSIC under permutations of the Y kernel matrix.
Only the row/column indices of Ky are permuted (avoiding recomputation).

## Usage

``` r
permutation_hsic_cpp(Kx, Ky, n_perm)
```

## Arguments

- Kx:

  n x n kernel matrix for X.

- Ky:

  n x n kernel matrix for Y.

- n_perm:

  Number of permutations.

## Value

Vector of n_perm HSIC values under permutation.
