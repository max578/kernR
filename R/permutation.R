#' Cluster-Based Permutation for bd-HSIC
#'
#' Permutes Y indices within clusters of similar conditional densities
#' p(x|z), ensuring valid exchangeability under the null.
#'
#' @param Kx n x n kernel matrix for X.
#' @param Ky n x n kernel matrix for Y.
#' @param weights Density ratio weights.
#' @param clusters Integer vector of cluster assignments.
#' @param n_permutations Number of permutations.
#'
#' @return Vector of permuted weighted HSIC statistics.
#' @keywords internal
cluster_permutation_hsic <- function(Kx, Ky, weights, clusters,
                                     n_permutations) {
  n <- nrow(Kx)
  n_clusters <- max(clusters)
  results <- numeric(n_permutations)

  for (p in seq_len(n_permutations)) {
    # Permute Y indices within each cluster
    perm <- seq_len(n)
    for (k in seq_len(n_clusters)) {
      idx <- which(clusters == k)
      if (length(idx) > 1) {
        perm[idx] <- idx[sample.int(length(idx))]
      }
    }
    Ky_perm <- Ky[perm, perm]
    results[p] <- weighted_hsic_stat_cpp(Kx, Ky_perm, weights)
  }

  results
}

#' Bin-Based Permutation for DR Tests
#'
#' Permutes treatment labels within propensity score bins.
#'
#' @param treatment Binary treatment vector.
#' @param propensity_scores Propensity score vector.
#' @param n_bins Number of bins.
#'
#' @return Permuted treatment vector.
#' @keywords internal
bin_permute_treatment <- function(treatment, propensity_scores, n_bins = 10L) {
  bins <- as.integer(cut(propensity_scores,
    breaks = n_bins,
    labels = FALSE,
    include.lowest = TRUE
  ))

  perm_idx <- stratified_permute_cpp(bins, n_bins)
  # C++ returns 0-indexed
  treatment[perm_idx + 1L]
}

#' Simple K-Means Clustering for Permutation Groups
#'
#' Clusters observations based on conditional density embeddings
#' using standard k-means on the density ratio weight space.
#'
#' @param weights Weight vector (density ratios or propensity scores).
#' @param z Confounder matrix.
#' @param n_clusters Number of clusters. If `"auto"`, selects by
#'   silhouette score (2 to 10 clusters).
#'
#' @return Integer vector of cluster assignments.
#' @keywords internal
cluster_observations <- function(weights, z, n_clusters = "auto") {
  z <- as.matrix(z)
  # Cluster on (weights, z) jointly
  features <- cbind(scale(weights), scale(z))

  if (identical(n_clusters, "auto")) {
    # Try 2 to min(10, n/5) clusters, pick by within-SS
    max_k <- min(10L, as.integer(nrow(features) / 5))
    max_k <- max(max_k, 2L)

    best_k <- 2L
    best_ratio <- Inf
    for (k in 2:max_k) {
      km <- stats::kmeans(features, centers = k, nstart = 5, iter.max = 50)
      ratio <- km$tot.withinss / km$totss
      if (ratio < best_ratio) {
        best_ratio <- ratio
        best_k <- k
      }
    }
    n_clusters <- best_k
  }

  km <- stats::kmeans(features, centers = as.integer(n_clusters),
    nstart = 10, iter.max = 100
  )
  km$cluster
}
