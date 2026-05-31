#' Create a Kernel Specification
#'
#' Constructs a kernel specification object used throughout `kernR` for
#' computing kernel matrices. Supports RBF (Gaussian), Matern, linear,
#' and polynomial kernels.
#'
#' @param type Character. Kernel type: `"rbf"` (default), `"matern"`,
#'   `"linear"`, or `"polynomial"`.
#' @param bandwidth Numeric or `"median"`. Lengthscale parameter for RBF
#'   and Matern kernels. If `"median"` (default), the median heuristic is
#'   used to select bandwidth automatically from the data.
#' @param nu Numeric. Smoothness parameter for the Matern kernel. Common
#'   choices: 0.5 (Laplace), 1.5, 2.5, Inf (RBF). Default is 2.5.
#' @param degree Integer. Degree for polynomial kernel. Default is 2.
#' @param offset Numeric. Offset for polynomial kernel. Default is 1.
#'
#' @return An object of class `"kernel_spec"`.
#'
#' @examples
#' # Default RBF kernel with median heuristic bandwidth
#' k <- kernel_spec()
#'
#' # RBF with fixed bandwidth
#' k <- kernel_spec("rbf", bandwidth = 1.0)
#'
#' # Matern kernel
#' k <- kernel_spec("matern", nu = 1.5)
#'
#' # Linear kernel (no bandwidth needed)
#' k <- kernel_spec("linear")
#'
#' @family kernel primitives
#' @export
kernel_spec <- function(type = c("rbf", "matern", "linear", "polynomial"),
                        bandwidth = "median",
                        nu = 2.5,
                        degree = 2L,
                        offset = 1.0) {
  type <- match.arg(type)

  if (type %in% c("rbf", "matern")) {
    if (!identical(bandwidth, "median") && !is.numeric(bandwidth)) {
      stop("`bandwidth` must be numeric or \"median\".", call. = FALSE)
    }
    if (is.numeric(bandwidth) && bandwidth <= 0) {
      stop("`bandwidth` must be positive.", call. = FALSE)
    }
  }

  if (type == "matern" && (!is.numeric(nu) || nu <= 0)) {
    stop("`nu` must be a positive number.", call. = FALSE)
  }

  structure(
    list(
      type = type,
      bandwidth = bandwidth,
      nu = nu,
      degree = as.integer(degree),
      offset = offset
    ),
    class = "kernel_spec"
  )
}

#' @export
print.kernel_spec <- function(x, ...) {
  cat("Kernel specification:\n")
  cat("  Type:", x$type, "\n")
  if (x$type %in% c("rbf", "matern")) {
    bw <- if (identical(x$bandwidth, "median")) "median heuristic" else x$bandwidth
    cat("  Bandwidth:", bw, "\n")
  }
  if (x$type == "matern") cat("  nu:", x$nu, "\n")
  if (x$type == "polynomial") {
    cat("  Degree:", x$degree, "\n")
    cat("  Offset:", x$offset, "\n")
  }
  invisible(x)
}

#' Resolve Kernel Bandwidth
#'
#' If `kernel$bandwidth` is `"median"`, compute the median heuristic from
#' the data. Otherwise return the fixed bandwidth.
#'
#' @param kernel A `kernel_spec` object.
#' @param x Numeric matrix (n x d).
#'
#' @return A `kernel_spec` with resolved numeric bandwidth.
#' @keywords internal
resolve_bandwidth <- function(kernel, x) {
  if (!identical(kernel$bandwidth, "median")) return(kernel)

  x <- as.matrix(x)
  bw <- median_bandwidth_cpp(x)
  kernel$bandwidth <- bw
  kernel
}
