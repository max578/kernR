#' Compute a Kernel Matrix
#'
#' Computes the kernel (Gram) matrix between two sets of observations.
#'
#' @param x Numeric matrix (n x d) or vector.
#' @param y Numeric matrix (m x d) or vector. If `NULL` (default),
#'   computes the kernel matrix of `x` with itself.
#' @param kernel A `kernel_spec` object. Default is RBF with median
#'   heuristic bandwidth.
#'
#' @return An n x m numeric matrix.
#'
#' @examples
#' x <- matrix(rnorm(100), 50, 2)
#' K <- kernel_matrix(x)
#' dim(K)  # 50 x 50
#'
#' @family kernel primitives
#' @export
kernel_matrix <- function(x, y = NULL, kernel = kernel_spec()) {
  x <- as.matrix(x)
  symmetric <- is.null(y)
  if (symmetric) y <- x else y <- as.matrix(y)

  if (ncol(x) != ncol(y)) {
    stop("`x` and `y` must have the same number of columns.", call. = FALSE)
  }

  # Resolve bandwidth if needed
  if (kernel$type %in% c("rbf", "matern") && identical(kernel$bandwidth, "median")) {
    kernel <- resolve_bandwidth(kernel, if (symmetric) x else rbind(x, y))
  }

  switch(
    kernel$type,
    rbf = rbf_kernel_matrix_cpp(x, y, kernel$bandwidth),
    matern = matern_kernel_matrix_cpp(x, y, kernel$bandwidth, kernel$nu),
    linear = linear_kernel_matrix_cpp(x, y),
    polynomial = polynomial_kernel_matrix_cpp(x, y, kernel$degree, kernel$offset),
    stop("Unknown kernel type: ", kernel$type, call. = FALSE)
  )
}

#' Centre a Kernel Matrix
#'
#' Centres a kernel matrix in feature space (double centring).
#'
#' @param K Square numeric matrix.
#'
#' @return Centred kernel matrix.
#' @keywords internal
centre_kernel_matrix <- function(K) {
  n <- nrow(K)
  H <- diag(n) - 1 / n
  H %*% K %*% H
}
