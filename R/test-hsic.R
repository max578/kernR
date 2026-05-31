#' HSIC Independence Test
#'
#' Tests whether two variables are independent using the Hilbert-Schmidt
#' Independence Criterion (HSIC). Uses a permutation test for inference.
#'
#' @param x Numeric vector, matrix, or data.frame. First variable.
#' @param y Numeric vector, matrix, or data.frame. Second variable.
#' @param kernel_x Kernel specification for `x`. Default is RBF with
#'   median heuristic.
#' @param kernel_y Kernel specification for `y`. Default is RBF with
#'   median heuristic.
#' @param n_permutations Integer. Number of permutations for the null
#'   distribution. Default is 500.
#' @param alpha Numeric. Significance level. Default is 0.05.
#' @param seed Integer or `NULL`. Random seed for reproducibility.
#'
#' @return An object of class `"kernel_test_result"` with components:
#'   \describe{
#'     \item{statistic}{The observed HSIC test statistic.}
#'     \item{p_value}{Permutation p-value.}
#'     \item{method}{`"HSIC"`.}
#'     \item{n}{Sample size.}
#'     \item{n_permutations}{Number of permutations used.}
#'     \item{null_distribution}{Vector of permuted HSIC values.}
#'     \item{kernel_x, kernel_y}{Kernel specifications used (with
#'       resolved bandwidths).}
#'     \item{call}{The matched call.}
#'   }
#'
#' @references
#' Gretton, A., Fukumizu, K., Teo, C. H., Song, L., Scholkopf, B., &
#' Smola, A. J. (2008). A kernel statistical test of independence.
#' *NeurIPS*, 20.
#'
#' @examples
#' set.seed(42)
#' n <- 200
#' x <- rnorm(n)
#'
#' # Dependent case
#' y_dep <- x^2 + rnorm(n, sd = 0.5)
#' result <- hsic_test(x, y_dep)
#' print(result)
#'
#' # Independent case
#' y_ind <- rnorm(n)
#' result <- hsic_test(x, y_ind)
#' print(result)
#'
#' @family independence and two-sample tests
#' @export
hsic_test <- function(x, y,
                      kernel_x = kernel_spec(),
                      kernel_y = kernel_spec(),
                      n_permutations = 500L,
                      alpha = 0.05,
                      seed = NULL) {
  cl <- match.call()

  # Input validation
  x <- as.matrix(x)
  y <- as.matrix(y)
  n <- nrow(x)

  if (nrow(y) != n) {
    stop("`x` and `y` must have the same number of observations.", call. = FALSE)
  }
  if (n < 10) {
    stop("At least 10 observations are required.", call. = FALSE)
  }
  n_permutations <- as.integer(n_permutations)
  if (n_permutations < 1) {
    stop("`n_permutations` must be a positive integer.", call. = FALSE)
  }

  if (!is.null(seed)) set.seed(seed)

  # Resolve kernels and compute kernel matrices
  kernel_x <- resolve_bandwidth(kernel_x, x)
  kernel_y <- resolve_bandwidth(kernel_y, y)

  Kx <- kernel_matrix(x, kernel = kernel_x)
  Ky <- kernel_matrix(y, kernel = kernel_y)

  # Observed statistic
  stat_obs <- hsic_stat_cpp(Kx, Ky)

  # Permutation null distribution
  null_dist <- permutation_hsic_cpp(Kx, Ky, n_permutations)

  # One-sided (upper-tail) permutation p-value with +1 correction
  p_value <- (1 + sum(null_dist >= stat_obs)) / (1 + n_permutations)

  structure(
    list(
      statistic = stat_obs,
      p_value = p_value,
      method = "HSIC",
      n = n,
      n_permutations = n_permutations,
      null_distribution = as.numeric(null_dist),
      ess = NA_real_,
      weights = NULL,
      kernel_x = kernel_x,
      kernel_y = kernel_y,
      call = cl
    ),
    class = "kernel_test_result"
  )
}
