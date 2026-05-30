# Bin-Based Permutation for DR Tests

Permutes treatment labels within propensity score bins.

## Usage

``` r
bin_permute_treatment(treatment, propensity_scores, n_bins = 10L)
```

## Arguments

- treatment:

  Binary treatment vector.

- propensity_scores:

  Propensity score vector.

- n_bins:

  Number of bins.

## Value

Permuted treatment vector.
