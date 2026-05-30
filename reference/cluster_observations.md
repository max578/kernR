# Simple K-Means Clustering for Permutation Groups

Clusters observations based on conditional density embeddings using
standard k-means on the density ratio weight space.

## Usage

``` r
cluster_observations(weights, z, n_clusters = "auto")
```

## Arguments

- weights:

  Weight vector (density ratios or propensity scores).

- z:

  Confounder matrix.

- n_clusters:

  Number of clusters. If `"auto"`, selects by silhouette score (2 to 10
  clusters).

## Value

Integer vector of cluster assignments.
