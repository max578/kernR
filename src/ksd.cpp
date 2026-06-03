// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
using namespace Rcpp;
using namespace arma;

// Stein-kernel machinery for the kernel Stein discrepancy goodness-of-fit
// test (ksd_test). Each entry of the returned matrix is the Stein kernel
// u_p(x_i, x_j) under the Langevin Stein operator for a target whose score
// (gradient of the log density) is supplied evaluated at the sample. The
// inverse multi-quadric (IMQ) and RBF base kernels both have closed-form
// Stein kernels; the trace term carries the dimension d explicitly so the
// estimator stays correct in the multivariate case. References inline in
// the R documentation of ksd_test().

//' Inverse multi-quadric Stein kernel matrix
//'
//' Builds the n x n Stein-kernel matrix u_p(x_i, x_j) for the IMQ base
//' kernel k(x, y) = (c^2 + ||x - y||^2)^beta under the Langevin Stein
//' operator, given the score (gradient of the log target density) evaluated
//' at each sample point.
//'
//' @param X Numeric matrix (n x d): the sample.
//' @param S Numeric matrix (n x d): the score evaluated row-wise at X.
//' @param beta Negative scalar exponent in (-1, 0).
//' @param c2 Squared offset c^2 (positive scalar).
//' @return Symmetric n x n Stein-kernel matrix.
//' @keywords internal
// [[Rcpp::export]]
arma::mat stein_kernel_imq_cpp(const arma::mat& X,
                               const arma::mat& S,
                               double beta,
                               double c2) {
  int n = X.n_rows;
  int d = X.n_cols;
  arma::mat H(n, n);

  for (int i = 0; i < n; i++) {
    for (int j = i; j < n; j++) {
      arma::rowvec r = X.row(i) - X.row(j);
      double rho = arma::dot(r, r);
      double g   = c2 + rho;
      double gb  = std::pow(g, beta);
      double gb1 = std::pow(g, beta - 1.0);
      double gb2 = std::pow(g, beta - 2.0);

      double si_sj = arma::dot(S.row(i), S.row(j));
      double si_r  = arma::dot(S.row(i), r);
      double sj_r  = arma::dot(S.row(j), r);

      double term_score = si_sj * gb;
      double term_cross = 2.0 * beta * gb1 * (sj_r - si_r);
      double term_trace = -2.0 * beta * d * gb1 -
                          4.0 * beta * (beta - 1.0) * gb2 * rho;

      double val = term_score + term_cross + term_trace;
      H(i, j) = val;
      H(j, i) = val;
    }
  }

  return H;
}

//' RBF Stein kernel matrix
//'
//' Builds the n x n Stein-kernel matrix u_p(x_i, x_j) for the Gaussian
//' base kernel k(x, y) = exp(-||x - y||^2 / (2 h^2)) under the Langevin
//' Stein operator, given the score evaluated at each sample point. The
//' bandwidth convention matches rbf_kernel_matrix_cpp.
//'
//' @param X Numeric matrix (n x d): the sample.
//' @param S Numeric matrix (n x d): the score evaluated row-wise at X.
//' @param h2 Squared bandwidth h^2 (positive scalar).
//' @return Symmetric n x n Stein-kernel matrix.
//' @keywords internal
// [[Rcpp::export]]
arma::mat stein_kernel_rbf_cpp(const arma::mat& X,
                               const arma::mat& S,
                               double h2) {
  int n = X.n_rows;
  int d = X.n_cols;
  arma::mat H(n, n);
  double inv_h2 = 1.0 / h2;

  for (int i = 0; i < n; i++) {
    for (int j = i; j < n; j++) {
      arma::rowvec r = X.row(i) - X.row(j);
      double rho = arma::dot(r, r);
      double k   = std::exp(-0.5 * rho * inv_h2);

      double si_sj = arma::dot(S.row(i), S.row(j));
      double si_r  = arma::dot(S.row(i), r);
      double sj_r  = arma::dot(S.row(j), r);

      double bracket = si_sj +
                       inv_h2 * (si_r - sj_r) +
                       inv_h2 * (d - rho * inv_h2);

      double val = k * bracket;
      H(i, j) = val;
      H(j, i) = val;
    }
  }

  return H;
}

//' Inverse multi-quadric Stein cross-kernel matrix (n x m)
//'
//' Rectangular counterpart of stein_kernel_imq_cpp: builds the n x m block
//' u_p(x_i, z_j) of the IMQ Stein kernel between the full sample X (with score
//' S) and a set of landmark points Z (with score Sm). Used to assemble the
//' Nystrom factorisation of the Stein-kernel matrix for ksd_test_nystrom().
//' The per-pair formula is identical to the symmetric builder; only the index
//' ranges differ (rows over X, columns over Z).
//'
//' @param X Numeric matrix (n x d): the full sample.
//' @param S Numeric matrix (n x d): the score evaluated row-wise at X.
//' @param Z Numeric matrix (m x d): the landmark points.
//' @param Sm Numeric matrix (m x d): the score evaluated row-wise at Z.
//' @param beta Negative scalar exponent in (-1, 0).
//' @param c2 Squared offset c^2 (positive scalar).
//' @return n x m Stein cross-kernel matrix.
//' @keywords internal
// [[Rcpp::export]]
arma::mat stein_kernel_imq_cross_cpp(const arma::mat& X,
                                     const arma::mat& S,
                                     const arma::mat& Z,
                                     const arma::mat& Sm,
                                     double beta,
                                     double c2) {
  int n = X.n_rows;
  int m = Z.n_rows;
  int d = X.n_cols;
  arma::mat C(n, m);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      arma::rowvec r = X.row(i) - Z.row(j);
      double rho = arma::dot(r, r);
      double g   = c2 + rho;
      double gb  = std::pow(g, beta);
      double gb1 = std::pow(g, beta - 1.0);
      double gb2 = std::pow(g, beta - 2.0);

      double si_sj = arma::dot(S.row(i), Sm.row(j));
      double si_r  = arma::dot(S.row(i), r);
      double sj_r  = arma::dot(Sm.row(j), r);

      double term_score = si_sj * gb;
      double term_cross = 2.0 * beta * gb1 * (sj_r - si_r);
      double term_trace = -2.0 * beta * d * gb1 -
                          4.0 * beta * (beta - 1.0) * gb2 * rho;

      C(i, j) = term_score + term_cross + term_trace;
    }
  }

  return C;
}

//' RBF Stein cross-kernel matrix (n x m)
//'
//' Rectangular counterpart of stein_kernel_rbf_cpp: the n x m block
//' u_p(x_i, z_j) of the Gaussian Stein kernel between the full sample X (with
//' score S) and landmark points Z (with score Sm). Bandwidth convention
//' matches stein_kernel_rbf_cpp.
//'
//' @param X Numeric matrix (n x d): the full sample.
//' @param S Numeric matrix (n x d): the score evaluated row-wise at X.
//' @param Z Numeric matrix (m x d): the landmark points.
//' @param Sm Numeric matrix (m x d): the score evaluated row-wise at Z.
//' @param h2 Squared bandwidth h^2 (positive scalar).
//' @return n x m Stein cross-kernel matrix.
//' @keywords internal
// [[Rcpp::export]]
arma::mat stein_kernel_rbf_cross_cpp(const arma::mat& X,
                                     const arma::mat& S,
                                     const arma::mat& Z,
                                     const arma::mat& Sm,
                                     double h2) {
  int n = X.n_rows;
  int m = Z.n_rows;
  int d = X.n_cols;
  arma::mat C(n, m);
  double inv_h2 = 1.0 / h2;

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      arma::rowvec r = X.row(i) - Z.row(j);
      double rho = arma::dot(r, r);
      double k   = std::exp(-0.5 * rho * inv_h2);

      double si_sj = arma::dot(S.row(i), Sm.row(j));
      double si_r  = arma::dot(S.row(i), r);
      double sj_r  = arma::dot(Sm.row(j), r);

      double bracket = si_sj +
                       inv_h2 * (si_r - sj_r) +
                       inv_h2 * (d - rho * inv_h2);

      C(i, j) = k * bracket;
    }
  }

  return C;
}

//' Factor-space wild bootstrap of the kernel Stein discrepancy null
//'
//' Low-rank counterpart of ksd_wild_bootstrap_cpp. Given a Nystrom factor F
//' (n x m) with \eqn{F F^\top \approx U} approximating the Stein-kernel
//' matrix, and the factor trace \eqn{\mathrm{tr}(F F^\top)}, draws n_boot wild-
//' bootstrap replicates of the degenerate U-statistic null via Rademacher
//' multipliers: each replicate is
//' \eqn{(\lVert F^\top w \rVert^2 - \mathrm{tr}) / (n (n - 1))}, computed in
//' O(n m) rather than O(n^2). Multipliers are drawn through R's RNG, so
//' callers honour set.seed().
//'
//' @param F Numeric n x m Nystrom factor of the Stein-kernel matrix.
//' @param tr Trace of F F^t (sum of squared factor entries).
//' @param n_boot Number of bootstrap replicates.
//' @return Vector of n_boot bootstrap KSD statistics.
//' @keywords internal
// [[Rcpp::export]]
arma::vec ksd_wild_bootstrap_factor_cpp(const arma::mat& F,
                                        double tr,
                                        int n_boot) {
  int n = F.n_rows;
  arma::vec results(n_boot);
  double denom = (double) n * (n - 1);

  for (int b = 0; b < n_boot; b++) {
    arma::vec w(n);
    for (int i = 0; i < n; i++) {
      w(i) = (R::unif_rand() < 0.5) ? -1.0 : 1.0;
    }
    arma::rowvec v = w.t() * F;            // 1 x m  ==  (F^t w)^t
    double quad = arma::dot(v, v);         // || F^t w ||^2  ==  w^t F F^t w
    results(b) = (quad - tr) / denom;
  }

  return results;
}

//' Wild bootstrap of the kernel Stein discrepancy null
//'
//' Given the Stein-kernel matrix H, draws n_boot wild-bootstrap replicates
//' of the degenerate U-statistic null via independent Rademacher multipliers
//' \eqn{W_i \in \{-1, +1\}}: each replicate is
//' \eqn{(1 / (n (n - 1))) \sum_{i \ne j} W_i W_j H_{ij}}. The diagonal is
//' excluded to match the unbiased U-statistic. Multipliers are drawn through
//' R's RNG (R::unif_rand within Rcpp's RNGScope), so callers honour
//' set.seed().
//'
//' @param H Symmetric n x n Stein-kernel matrix.
//' @param n_boot Number of bootstrap replicates.
//' @return Vector of n_boot bootstrap KSD statistics.
//' @keywords internal
// [[Rcpp::export]]
arma::vec ksd_wild_bootstrap_cpp(const arma::mat& H, int n_boot) {
  int n = H.n_rows;
  arma::vec results(n_boot);
  double denom = (double) n * (n - 1);
  double tr = arma::trace(H);

  for (int b = 0; b < n_boot; b++) {
    arma::vec w(n);
    for (int i = 0; i < n; i++) {
      w(i) = (R::unif_rand() < 0.5) ? -1.0 : 1.0;
    }
    double quad = arma::as_scalar(w.t() * H * w);
    results(b) = (quad - tr) / denom;
  }

  return results;
}
