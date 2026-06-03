#' Weighted HSIC Statistic
#'
#' Computes the weighted Hilbert--Schmidt Independence Criterion (HSIC) between
#' two pre-computed kernel matrices,
#' \deqn{\frac{1}{W^2} \sum_{i,j} w_i w_j (K_x^c)_{ij} (K_y^c)_{ij},}
#' where `K_x^c` and `K_y^c` are weight-centred and \eqn{W = \sum_i w_i}. This is
#' the exact statistic [bd_hsic_test()] accumulates over its permutation null,
#' exposed so that downstream methods -- for example theory-anchored causal
#' inference -- can build a bd-HSIC-compatible statistic against a custom
#' reference distribution without re-implementing the engine.
#'
#' With uniform weights this reduces to the (biased) HSIC estimator; supplying
#' density-ratio weights \eqn{w(t, z)} yields the backdoor-adjusted statistic
#' used for causal testing.
#'
#' @param Kx Numeric `n` by `n` kernel matrix for the first variable.
#' @param Ky Numeric `n` by `n` kernel matrix for the second variable, with the
#'   same dimensions as `Kx`.
#' @param w Non-negative numeric weight vector of length `n`. Defaults to
#'   uniform weights, giving the unweighted HSIC.
#'
#' @return A single numeric value: the weighted HSIC statistic.
#'
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(40), ncol = 2)
#' y <- matrix(rnorm(40), ncol = 2)
#' Kx <- kernel_matrix(x, kernel = resolve_bandwidth(kernel_spec(), x))
#' Ky <- kernel_matrix(y, kernel = resolve_bandwidth(kernel_spec(), y))
#' weighted_hsic_stat(Kx, Ky)
#'
#' @seealso [bd_hsic_test()] for the full causal test that uses this statistic;
#'   [resolve_bandwidth()] and [kernel_matrix()] for building the kernel inputs.
#' @family kernel primitives
#' @export
weighted_hsic_stat <- function(Kx, Ky, w = rep(1, nrow(Kx))) {
  Kx <- as.matrix(Kx)
  Ky <- as.matrix(Ky)
  n <- nrow(Kx)

  if (ncol(Kx) != n) {
    stop("`Kx` must be a square matrix.", call. = FALSE)
  }
  if (!identical(dim(Ky), dim(Kx))) {
    stop("`Ky` must have the same dimensions as `Kx`.", call. = FALSE)
  }

  w <- as.numeric(w)
  if (length(w) != n) {
    stop("`w` must have length `nrow(Kx)`.", call. = FALSE)
  }
  if (anyNA(w) || any(w < 0)) {
    stop("`w` must be non-negative and free of missing values.", call. = FALSE)
  }
  if (sum(w) <= 0) {
    stop("`w` must have a positive sum.", call. = FALSE)
  }

  weighted_hsic_stat_cpp(Kx, Ky, w)
}
