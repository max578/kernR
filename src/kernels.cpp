// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//' Compute the median pairwise distance (bandwidth heuristic)
//'
//' @param x Numeric matrix (n x d).
//' @return Scalar: sqrt(median of squared pairwise distances).
//' @keywords internal
// [[Rcpp::export]]
double median_bandwidth_cpp(const arma::mat& x) {
  int n = x.n_rows;

  // For large n, subsample to keep computation manageable

int n_sub = (n > 2500) ? 2500 : n;
  arma::uvec idx;
  arma::mat xs;
  if (n_sub < n) {
    idx = arma::randperm(n, n_sub);
    xs = x.rows(idx);
  } else {
    xs = x;
  }

  int ns = xs.n_rows;
  // Compute squared pairwise distances
  arma::vec sq_norms = arma::sum(arma::square(xs), 1);
  arma::mat dist2 = arma::repmat(sq_norms, 1, ns) +
                     arma::repmat(sq_norms.t(), ns, 1) -
                     2.0 * xs * xs.t();

  // Extract upper triangle
  std::vector<double> upper_tri;
  upper_tri.reserve(ns * (ns - 1) / 2);
  for (int i = 0; i < ns; i++) {
    for (int j = i + 1; j < ns; j++) {
      double d = dist2(i, j);
      if (d > 0) upper_tri.push_back(d);
    }
  }

  if (upper_tri.empty()) return 1.0;

  // Compute median
  size_t mid = upper_tri.size() / 2;
  std::nth_element(upper_tri.begin(), upper_tri.begin() + mid, upper_tri.end());
  double med = upper_tri[mid];

  return std::sqrt(med);
}

//' RBF (Gaussian) kernel matrix
//'
//' K(x, y) = exp(-||x - y||^2 / (2 * bandwidth^2))
//'
//' @param x Numeric matrix (n x d).
//' @param y Numeric matrix (m x d).
//' @param bandwidth Positive scalar.
//' @return n x m kernel matrix.
//' @keywords internal
// [[Rcpp::export]]
arma::mat rbf_kernel_matrix_cpp(const arma::mat& x,
                                 const arma::mat& y,
                                 double bandwidth) {
  double gamma = -1.0 / (2.0 * bandwidth * bandwidth);

  arma::vec sq_x = arma::sum(arma::square(x), 1);
  arma::vec sq_y = arma::sum(arma::square(y), 1);

  arma::mat dist2 = arma::repmat(sq_x, 1, y.n_rows) +
                     arma::repmat(sq_y.t(), x.n_rows, 1) -
                     2.0 * x * y.t();

  // Clamp negative values from numerical error
  dist2.clamp(0.0, arma::datum::inf);

  return arma::exp(gamma * dist2);
}

//' Matern kernel matrix
//'
//' Supports nu = 0.5, 1.5, 2.5, and Inf (RBF).
//'
//' @param x Numeric matrix (n x d).
//' @param y Numeric matrix (m x d).
//' @param bandwidth Positive scalar.
//' @param nu Smoothness parameter.
//' @return n x m kernel matrix.
//' @keywords internal
// [[Rcpp::export]]
arma::mat matern_kernel_matrix_cpp(const arma::mat& x,
                                    const arma::mat& y,
                                    double bandwidth,
                                    double nu) {
  arma::vec sq_x = arma::sum(arma::square(x), 1);
  arma::vec sq_y = arma::sum(arma::square(y), 1);

  arma::mat dist2 = arma::repmat(sq_x, 1, y.n_rows) +
                     arma::repmat(sq_y.t(), x.n_rows, 1) -
                     2.0 * x * y.t();
  dist2.clamp(0.0, arma::datum::inf);
  arma::mat dist = arma::sqrt(dist2);

  double l = bandwidth;

  if (nu == 0.5) {
    // Exponential kernel
    return arma::exp(-dist / l);
  } else if (std::abs(nu - 1.5) < 1e-10) {
    arma::mat s3 = std::sqrt(3.0) * dist / l;
    return (1.0 + s3) % arma::exp(-s3);
  } else if (std::abs(nu - 2.5) < 1e-10) {
    arma::mat s5 = std::sqrt(5.0) * dist / l;
    return (1.0 + s5 + (5.0 / 3.0) * dist2 / (l * l)) % arma::exp(-s5);
  } else {
    // Fall back to RBF for nu = Inf or unsupported values
    double gamma = -1.0 / (2.0 * l * l);
    return arma::exp(gamma * dist2);
  }
}

//' Linear kernel matrix
//'
//' K(x, y) = x %*% t(y)
//'
//' @param x Numeric matrix (n x d).
//' @param y Numeric matrix (m x d).
//' @return n x m kernel matrix.
//' @keywords internal
// [[Rcpp::export]]
arma::mat linear_kernel_matrix_cpp(const arma::mat& x,
                                    const arma::mat& y) {
  return x * y.t();
}

//' Polynomial kernel matrix
//'
//' K(x, y) = (x %*% t(y) + offset)^degree
//'
//' @param x Numeric matrix (n x d).
//' @param y Numeric matrix (m x d).
//' @param degree Integer degree.
//' @param offset Scalar offset.
//' @return n x m kernel matrix.
//' @keywords internal
// [[Rcpp::export]]
arma::mat polynomial_kernel_matrix_cpp(const arma::mat& x,
                                        const arma::mat& y,
                                        int degree,
                                        double offset) {
  arma::mat K = x * y.t() + offset;
  return arma::pow(K, degree);
}
