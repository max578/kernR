// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
using namespace Rcpp;
using namespace arma;

// Fisher-Yates permutation of 0..n-1 driven by R's RNG, so callers honour
// set.seed(). Rcpp wraps every exported entry point in an RNGScope
// (rng = true by default), syncing R's RNG state, so R::unif_rand() here
// advances the same stream set.seed() controls. Replaces arma::randperm,
// whose internal Mersenne-Twister ignores set.seed() and silently broke
// the package's documented `seed=` reproducibility contract.
arma::uvec r_randperm(int n) {
  arma::uvec idx(n);
  for (int i = 0; i < n; i++) idx(i) = (arma::uword) i;
  for (int i = n - 1; i > 0; i--) {
    int j = (int) std::floor(R::unif_rand() * (i + 1));
    if (j > i) j = i;  // guard the unif_rand() == 1.0 edge
    arma::uword tmp = idx(i);
    idx(i) = idx(j);
    idx(j) = tmp;
  }
  return idx;
}

//' Permutation HSIC: compute HSIC for many Y permutations
//'
//' Efficiently computes HSIC under permutations of the Y kernel matrix.
//' Only the row/column indices of Ky are permuted (avoiding recomputation).
//'
//' @param Kx n x n kernel matrix for X.
//' @param Ky n x n kernel matrix for Y.
//' @param n_perm Number of permutations.
//' @return Vector of n_perm HSIC values under permutation.
//' @keywords internal
// [[Rcpp::export]]
arma::vec permutation_hsic_cpp(const arma::mat& Kx,
                                const arma::mat& Ky,
                                int n_perm) {
  int n = Kx.n_rows;
  arma::vec results(n_perm);

  for (int p = 0; p < n_perm; p++) {
    arma::uvec perm = r_randperm(n);
    arma::mat Ky_perm = Ky(perm, perm);
    results(p) = 0.0;

    // Centred HSIC inline
    double n2 = (double)n * n;
    arma::vec ones_n = arma::ones<arma::vec>(n);
    arma::vec Kx_rm = Kx * ones_n / n;
    arma::vec Ky_rm = Ky_perm * ones_n / n;
    double Kx_m = arma::accu(Kx) / n2;
    double Ky_m = arma::accu(Ky_perm) / n2;

    double hsic = 0.0;
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        double kxc = Kx(i, j) - Kx_rm(i) - Kx_rm(j) + Kx_m;
        double kyc = Ky_perm(i, j) - Ky_rm(i) - Ky_rm(j) + Ky_m;
        hsic += kxc * kyc;
      }
    }
    results(p) = hsic / n2;
  }

  return results;
}

//' Permutation MMD: compute MMD^2 for many permutations of pooled data
//'
//' Given the full (n+m) x (n+m) kernel matrix of the pooled sample,
//' permute the assignment into two groups and compute MMD^2.
//'
//' @param K_pool (n+m) x (n+m) kernel matrix.
//' @param n Size of first sample.
//' @param m Size of second sample.
//' @param n_perm Number of permutations.
//' @return Vector of n_perm MMD^2 values under permutation.
//' @keywords internal
// [[Rcpp::export]]
arma::vec permutation_mmd_cpp(const arma::mat& K_pool,
                               int n, int m, int n_perm) {
  int N = n + m;
  arma::vec results(n_perm);

  for (int p = 0; p < n_perm; p++) {
    arma::uvec perm = r_randperm(N);
    arma::uvec idx_x = perm.head(n);
    arma::uvec idx_y = perm.tail(m);

    double sum_xx = 0.0, sum_yy = 0.0, sum_xy = 0.0;

    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        if (i != j) sum_xx += K_pool(idx_x(i), idx_x(j));
      }
    }
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < m; j++) {
        if (i != j) sum_yy += K_pool(idx_y(i), idx_y(j));
      }
    }
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < m; j++) {
        sum_xy += K_pool(idx_x(i), idx_y(j));
      }
    }

    results(p) = sum_xx / ((double)n * (n - 1)) +
                 sum_yy / ((double)m * (m - 1)) -
                 2.0 * sum_xy / ((double)n * m);
  }

  return results;
}

//' Stratified permutation of indices within bins
//'
//' Given bin assignments, permute indices within each bin.
//'
//' @param bins Integer vector of bin assignments (1-indexed).
//' @param n_bins Number of bins.
//' @return Permuted index vector (0-indexed for C++).
//' @keywords internal
// [[Rcpp::export]]
arma::uvec stratified_permute_cpp(const arma::ivec& bins, int n_bins) {
  int n = bins.n_elem;
  arma::uvec result(n);

  for (int b = 1; b <= n_bins; b++) {
    arma::uvec idx = arma::find(bins == b);
    arma::uvec perm = r_randperm((int) idx.n_elem);
    for (unsigned int k = 0; k < idx.n_elem; k++) {
      result(idx(k)) = idx(perm(k));
    }
  }

  return result;
}

//' Permutation k-sample MMD: summed pairwise unbiased MMD^2 under joint relabel
//'
//' Given the pooled (N x N) kernel matrix of K stacked groups and their sizes,
//' draws n_perm joint relabelings of the pooled sample into the original group
//' sizes and returns, for each, the summed pairwise unbiased MMD^2 statistic
//' \eqn{\sum_{a < b} \mathrm{MMD}^2_u(a, b)}. The single shared relabeling per
//' replicate (not independent per-pair permutation) is what makes this a valid
//' k-sample null. Relabeling uses r_randperm, so callers honour set.seed().
//'
//' @param K_pool (N x N) kernel matrix of the row-stacked groups.
//' @param sizes Integer vector of the K group sizes (summing to N).
//' @param n_perm Number of permutations.
//' @return Vector of n_perm summed-pairwise MMD^2 values under the null.
//' @keywords internal
// [[Rcpp::export]]
arma::vec permutation_ksample_mmd_cpp(const arma::mat& K_pool,
                                       const arma::ivec& sizes,
                                       int n_perm) {
  int K = sizes.n_elem;
  int N = K_pool.n_rows;

  arma::ivec offset(K);
  int acc = 0;
  for (int k = 0; k < K; k++) {
    offset(k) = acc;
    acc += sizes(k);
  }

  arma::vec results(n_perm);

  for (int p = 0; p < n_perm; p++) {
    arma::uvec perm = r_randperm(N);

    // Within-group sums (i != j), one block per group
    arma::vec within(K, arma::fill::zeros);
    for (int g = 0; g < K; g++) {
      int s = offset(g);
      int ng = sizes(g);
      for (int i = 0; i < ng; i++) {
        for (int j = 0; j < ng; j++) {
          if (i != j) within(g) += K_pool(perm(s + i), perm(s + j));
        }
      }
    }

    // Each within-group term appears in (K - 1) pairs
    double stat = 0.0;
    for (int g = 0; g < K; g++) {
      double ng = (double) sizes(g);
      stat += (K - 1) * within(g) / (ng * (ng - 1));
    }

    // Between-group sums over the upper triangle of group pairs
    for (int a = 0; a < K; a++) {
      for (int b = a + 1; b < K; b++) {
        int sa = offset(a), na = sizes(a);
        int sb = offset(b), nb = sizes(b);
        double between = 0.0;
        for (int i = 0; i < na; i++) {
          for (int j = 0; j < nb; j++) {
            between += K_pool(perm(sa + i), perm(sb + j));
          }
        }
        stat -= 2.0 * between / ((double) na * nb);
      }
    }

    results(p) = stat;
  }

  return results;
}
