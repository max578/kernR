#' MMD Two-Sample Test
#'
#' Tests whether two samples come from the same distribution using the
#' Maximum Mean Discrepancy (MMD). Uses a permutation test for inference.
#'
#' @param x Numeric vector, matrix, or data.frame. First sample.
#' @param y Numeric vector, matrix, or data.frame. Second sample.
#' @param kernel Kernel specification. Default is RBF with median heuristic.
#' @param n_permutations Integer. Number of permutations. Default is 500.
#' @param alpha Numeric. Significance level. Default is 0.05.
#' @param seed Integer or `NULL`. Random seed for reproducibility.
#'
#' @return An object of class `"kernel_test_result"` with components:
#'   \describe{
#'     \item{statistic}{The observed MMD^2 statistic (unbiased).}
#'     \item{p_value}{Permutation p-value.}
#'     \item{method}{`"MMD"`.}
#'     \item{n}{Total sample size (n_x + n_y).}
#'     \item{n_permutations}{Number of permutations used.}
#'     \item{null_distribution}{Vector of permuted MMD^2 values.}
#'     \item{kernel_x}{Kernel specification used.}
#'     \item{call}{The matched call.}
#'   }
#'
#' @references
#' Gretton, A., Borgwardt, K. M., Rasch, M. J., Scholkopf, B., & Smola,
#' A. (2012). A kernel two-sample test. *JMLR*, 13, 723-773.
#'
#' @examples
#' set.seed(42)
#'
#' # Same distribution
#' x <- matrix(rnorm(200), 100, 2)
#' y <- matrix(rnorm(200), 100, 2)
#' result <- mmd_test(x, y)
#' print(result)
#'
#' # Different distributions
#' y_shifted <- matrix(rnorm(200, mean = 1), 100, 2)
#' result <- mmd_test(x, y_shifted)
#' print(result)
#'
#' @export
mmd_test <- function(x, y,
                     kernel = kernel_spec(),
                     n_permutations = 500L,
                     alpha = 0.05,
                     seed = NULL) {
  cl <- match.call()

  x <- as.matrix(x)
  y <- as.matrix(y)
  nx <- nrow(x)
  ny <- nrow(y)

  if (ncol(x) != ncol(y)) {
    stop("`x` and `y` must have the same number of columns.", call. = FALSE)
  }
  if (nx < 5 || ny < 5) {
    stop("Each sample must have at least 5 observations.", call. = FALSE)
  }
  n_permutations <- as.integer(n_permutations)

  if (!is.null(seed)) set.seed(seed)

  # Resolve bandwidth on pooled data
  pooled <- rbind(x, y)
  kernel <- resolve_bandwidth(kernel, pooled)

  # Compute kernel matrices
  K_pool <- kernel_matrix(pooled, kernel = kernel)
  Kxx <- K_pool[1:nx, 1:nx, drop = FALSE]
  Kyy <- K_pool[(nx + 1):(nx + ny), (nx + 1):(nx + ny), drop = FALSE]
  Kxy <- K_pool[1:nx, (nx + 1):(nx + ny), drop = FALSE]

  # Observed statistic
  stat_obs <- mmd2_unbiased_cpp(Kxx, Kyy, Kxy)

  # Permutation null distribution
  null_dist <- permutation_mmd_cpp(K_pool, nx, ny, n_permutations)

  # One-sided (upper-tail) permutation p-value with +1 correction
  p_value <- (1 + sum(null_dist >= stat_obs)) / (1 + n_permutations)

  structure(
    list(
      statistic = stat_obs,
      p_value = p_value,
      method = "MMD",
      n = nx + ny,
      n_permutations = n_permutations,
      null_distribution = as.numeric(null_dist),
      ess = NA_real_,
      weights = NULL,
      kernel_x = kernel,
      kernel_y = NULL,
      call = cl
    ),
    class = "kernel_test_result"
  )
}
