# test-ksd.R -- Kernel Stein discrepancy goodness-of-fit / calibration test.
#
# Implements ksd_test(), a one-sample kernel goodness-of-fit test that asks
# whether a sample is consistent with a target distribution known only
# through its score (the gradient of its log density). This is the
# score-based complement to the sample-based two-sample mmd_test(): no
# reference sample is needed, only the target's score, so the target may be
# unnormalised. The wild-bootstrap calibration of the degenerate U-statistic
# null follows Chwialkowski et al. (2016); the inverse multi-quadric default
# kernel follows Gorham & Mackey (2017). gaussian_score() is a convenience
# score factory for the common multivariate-normal target.

# Score factory --------------------------------------------------------------

#' Score function for a multivariate normal target
#'
#' Builds a score function -- the gradient of the log density,
#' \eqn{\nabla_x \log p(x) = -\Sigma^{-1}(x - \mu)} -- for a multivariate
#' normal target, suitable for the `score` argument of [ksd_test()].
#'
#' The returned closure accepts the sample matrix and returns the score
#' evaluated row-wise. Leaving `mean` or `sigma` at `NULL` defaults them to
#' the zero vector and the identity matrix of the dimension seen at call
#' time, so `gaussian_score()` with no arguments is the standard-normal
#' score in any dimension.
#'
#' @param mean Numeric vector of length `d`, the target mean. `NULL`
#'   (default) uses the zero vector of the sample dimension.
#' @param sigma Numeric `d x d` covariance matrix. `NULL` (default) uses the
#'   identity matrix of the sample dimension. Must be symmetric and
#'   invertible.
#'
#' @returns A function of one argument (an `n x d` numeric matrix) returning
#'   the `n x d` matrix of scores. The mean and covariance are captured by
#'   the closure.
#'
#' @references
#' Liu, Q., Lee, J. D., & Jordan, M. I. (2016). A kernelized Stein
#' discrepancy for goodness-of-fit tests. *Proceedings of the 33rd
#' International Conference on Machine Learning*, PMLR 48, 276-284.
#'
#' @examples
#' set.seed(1)
#' x <- matrix(stats::rnorm(200L), ncol = 2L)
#'
#' # Standard-normal target in two dimensions
#' s0 <- gaussian_score()
#' str(s0(x))
#'
#' # Correlated-normal target
#' sig <- matrix(c(1, 0.5, 0.5, 1), nrow = 2L)
#' s1 <- gaussian_score(mean = c(0, 0), sigma = sig)
#'
#' @seealso [ksd_test()]
#' @family goodness-of-fit tests
#' @author Max Moldovan, \email{max.moldovan@@adelaide.edu.au}
#' @export
gaussian_score <- function(mean = NULL, sigma = NULL) {
  force(mean)
  force(sigma)

  function(x) {
    x <- as.matrix(x)
    d <- ncol(x)

    mu <- if (is.null(mean)) rep(0, d) else as.numeric(mean)
    if (length(mu) != d) {
      stop("`mean` must have length ", d, " to match the sample dimension.",
           call. = FALSE)
    }

    if (is.null(sigma)) {
      prec <- diag(d)
    } else {
      sigma <- as.matrix(sigma)
      if (nrow(sigma) != d || ncol(sigma) != d) {
        stop("`sigma` must be a ", d, " x ", d, " covariance matrix.",
             call. = FALSE)
      }
      if (!isSymmetric(unname(sigma))) {
        stop("`sigma` must be a symmetric covariance matrix.", call. = FALSE)
      }
      prec <- tryCatch(
        solve(sigma),
        error = function(e) {
          stop("`sigma` must be invertible: ", conditionMessage(e),
               call. = FALSE)
        }
      )
    }

    centred <- sweep(x, 2L, mu, "-")
    -(centred %*% prec)
  }
}

# Goodness-of-fit test -------------------------------------------------------

#' Kernel Stein Discrepancy Goodness-of-Fit Test
#'
#' Tests whether a sample is consistent with a target distribution `p`, where
#' `p` is supplied through its score \eqn{\nabla_x \log p(x)} rather than a
#' reference sample. The statistic is the (unbiased) kernel Stein discrepancy
#' (KSD); calibration uses a wild bootstrap of the degenerate U-statistic
#' null. Because only the score enters, the target may be known up to an
#' unknown normalising constant.
#'
#' KSD is the score-based, one-sample complement to [mmd_test()]: where MMD
#' compares two samples, KSD compares a sample against a *density*. The
#' calibration framing is direct -- given posterior or ensemble draws and the
#' score of the distribution they claim to represent, KSD asks whether the
#' draws actually follow that distribution. It is sensitive to mean,
#' variance, and tail mis-specification.
#'
#' The default base kernel is the inverse multi-quadric (IMQ),
#' \eqn{k(x, y) = (c^2 + \lVert x - y \rVert^2)^\beta} with
#' \eqn{\beta \in (-1, 0)}. Gorham & Mackey (2017) show the IMQ Stein
#' discrepancy detects non-convergence in regimes where the Gaussian (RBF)
#' Stein discrepancy is blind, particularly as dimension grows; the RBF base
#' kernel remains available via `kernel = "rbf"`. The offset `c` (IMQ) and
#' bandwidth `h` (RBF) default to the median heuristic over the sample.
#'
#' Reproducibility: the wild-bootstrap multipliers are drawn through R's RNG,
#' so a non-`NULL` `seed` makes the p-value reproducible under the active RNG
#' kind (the R default Mersenne-Twister unless changed by the caller).
#'
#' The current implementation materialises the `n x n` Stein-kernel matrix,
#' so memory scales as `O(n^2)`; for very large samples, thin first or use
#' the two-sample [mmd_test()] against reference draws.
#'
#' @param x Numeric vector, matrix, or data.frame. The `n x d` sample to
#'   test, one observation per row. At least five rows are required.
#' @param score Either `NULL` or a function. When `NULL` (default) the target
#'   is the standard multivariate normal, with score \eqn{-x}. When a
#'   function, it must accept the `n x d` sample matrix and return the
#'   `n x d` matrix of scores \eqn{\nabla_x \log p(x)} evaluated row-wise;
#'   see [gaussian_score()] for the multivariate-normal factory.
#' @param kernel Character. Base kernel for the Stein kernel: `"imq"`
#'   (inverse multi-quadric, default) or `"rbf"` (Gaussian).
#' @param beta Numeric in `(-1, 0)`. Exponent of the IMQ kernel; ignored when
#'   `kernel = "rbf"`. Default `-0.5`.
#' @param bandwidth Numeric or `"median"`. The IMQ offset `c` or the RBF
#'   bandwidth `h`. `"median"` (default) uses the median-heuristic bandwidth
#'   of the sample.
#' @param n_boot Integer. Number of wild-bootstrap replicates for the null.
#'   Default `1000`.
#' @param alpha Numeric in `(0, 1)`. Significance level for the verdict.
#'   Default `0.05`.
#' @param seed Integer or `NULL`. Random seed for reproducibility.
#'
#' @returns An object of class `c("ksd_test", "kernel_test_result")` carrying
#'   the standard `kernel_test_result` fields plus:
#'   \describe{
#'     \item{statistic}{The unbiased KSD U-statistic.}
#'     \item{p_value}{Upper-tail wild-bootstrap p-value (with `+1`
#'       correction).}
#'     \item{stein_kernel}{Base kernel used (`"imq"` or `"rbf"`).}
#'     \item{beta}{IMQ exponent, or `NA` for the RBF kernel.}
#'     \item{bandwidth}{Resolved IMQ offset or RBF bandwidth.}
#'     \item{surprise_bits}{Shannon-information surprise `-log2(p_value)`.}
#'     \item{alpha, reject}{Verdict level and `p_value <= alpha`.}
#'   }
#'
#' @references
#' Liu, Q., Lee, J. D., & Jordan, M. I. (2016). A kernelized Stein
#' discrepancy for goodness-of-fit tests. *Proceedings of the 33rd
#' International Conference on Machine Learning*, PMLR 48, 276-284.
#'
#' Chwialkowski, K., Strathmann, H., & Gretton, A. (2016). A kernel test of
#' goodness of fit. *Proceedings of the 33rd International Conference on
#' Machine Learning*, PMLR 48, 2606-2615.
#'
#' Gorham, J., & Mackey, L. (2017). Measuring sample quality with kernels.
#' *Proceedings of the 34th International Conference on Machine Learning*,
#' PMLR 70, 1292-1301.
#'
#' @examples
#' set.seed(1)
#'
#' # Well-specified: standard-normal sample against standard-normal target
#' x_ok <- matrix(stats::rnorm(400L), ncol = 2L)
#' fit_ok <- ksd_test(x_ok, n_boot = 199L, seed = 1L)
#' fit_ok
#'
#' # Mis-specified: mean-shifted sample against the same target
#' x_bad <- x_ok + 1
#' fit_bad <- ksd_test(x_bad, n_boot = 199L, seed = 1L)
#' fit_bad
#'
#' # Explicit non-standard target via the Gaussian score factory
#' sig <- matrix(c(1, 0.6, 0.6, 1), nrow = 2L)
#' x_cor <- x_ok %*% chol(sig)
#' ksd_test(x_cor, score = gaussian_score(sigma = sig),
#'          n_boot = 199L, seed = 1L)
#'
#' @seealso [mmd_test()], [mmd_ppc()], [gaussian_score()]
#' @family goodness-of-fit tests
#' @author Max Moldovan, \email{max.moldovan@@adelaide.edu.au}
#' @export
ksd_test <- function(x,
                     score = NULL,
                     kernel = c("imq", "rbf"),
                     beta = -0.5,
                     bandwidth = "median",
                     n_boot = 1000L,
                     alpha = 0.05,
                     seed = NULL) {
  cl <- match.call()
  kernel <- match.arg(kernel)

  x <- as.matrix(x)
  storage.mode(x) <- "double"
  .check_ksd_input(x, kernel, beta, n_boot, alpha)
  n_boot <- as.integer(n_boot)

  # Resolve the score and bandwidth --------------------------------------
  score_matrix <- .resolve_score(score, x)
  bw <- .ksd_bandwidth(bandwidth, x)

  n <- nrow(x)
  d <- ncol(x)

  if (!is.null(seed)) set.seed(seed)

  # Stein-kernel matrix and the observed U-statistic ---------------------
  h_mat <- switch(
    kernel,
    imq = stein_kernel_imq_cpp(x, score_matrix, beta, bw * bw),
    rbf = stein_kernel_rbf_cpp(x, score_matrix, bw * bw)
  )

  denom <- as.numeric(n) * (n - 1)
  stat_obs <- (sum(h_mat) - sum(diag(h_mat))) / denom

  # Wild-bootstrap null and the upper-tail p-value -----------------------
  null_dist <- ksd_wild_bootstrap_cpp(h_mat, n_boot)
  p_value <- (1 + sum(null_dist >= stat_obs)) / (1 + n_boot)

  kernel_x <- structure(
    list(type = kernel, bandwidth = bw),
    class = "kernel_spec"
  )

  structure(
    list(
      statistic = stat_obs,
      p_value = p_value,
      method = "KSD",
      n = n,
      n_permutations = n_boot,
      null_distribution = as.numeric(null_dist),
      ess = NA_real_,
      weights = NULL,
      kernel_x = kernel_x,
      kernel_y = NULL,
      stein_kernel = kernel,
      beta = if (identical(kernel, "imq")) beta else NA_real_,
      bandwidth = bw,
      dimension = d,
      surprise_bits = -log2(p_value),
      alpha = alpha,
      reject = p_value <= alpha,
      call = cl
    ),
    class = c("ksd_test", "kernel_test_result")
  )
}

#' @export
print.ksd_test <- function(x, ...) {
  NextMethod()
  cat("Goodness-of-fit verdict\n")
  cat("  Stein kernel: ", x$stein_kernel,
      if (identical(x$stein_kernel, "imq")) {
        paste0(" (beta = ", formatC(x$beta, digits = 3L, format = "g"), ")")
      } else {
        ""
      },
      "\n", sep = "")
  cat("  Bandwidth:    ",
      formatC(x$bandwidth, digits = 4L, format = "g"),
      " (median heuristic)\n", sep = "")
  cat("  Bootstrap:    wild, B = ", x$n_permutations, "\n", sep = "")
  cat("  Surprise:     ",
      formatC(x$surprise_bits, digits = 3L, format = "f"), " bits\n",
      sep = "")
  cat("  Verdict:      ",
      if (isTRUE(x$reject)) {
        "REJECT (sample inconsistent with target)"
      } else {
        "consistent with target"
      },
      "\n\n", sep = "")
  invisible(x)
}

# Internal validation and resolution -----------------------------------------

#' Validate the ksd_test() inputs
#'
#' @param x Coerced numeric sample matrix.
#' @param kernel Resolved kernel string.
#' @param beta IMQ exponent.
#' @param n_boot Requested bootstrap count.
#' @param alpha Verdict level.
#' @returns `invisible(NULL)`; called for its error side effects.
#' @noRd
#' @keywords internal
.check_ksd_input <- function(x, kernel, beta, n_boot, alpha) {
  if (!is.numeric(x) || any(!is.finite(x))) {
    stop("`x` must be numeric and contain only finite values.",
         call. = FALSE)
  }
  if (nrow(x) < 5L) {
    stop("`x` must have at least 5 rows (observations).", call. = FALSE)
  }
  if (identical(kernel, "imq")) {
    if (!is.numeric(beta) || length(beta) != 1L || beta <= -1 || beta >= 0) {
      stop("`beta` must be a single number in (-1, 0).", call. = FALSE)
    }
  }
  if (!is.numeric(n_boot) || length(n_boot) != 1L || n_boot < 1) {
    stop("`n_boot` must be a single positive integer.", call. = FALSE)
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a single number in (0, 1).", call. = FALSE)
  }
  invisible(NULL)
}

#' Evaluate and validate the score at the sample
#'
#' @param score `NULL` (standard-normal score) or a user score function.
#' @param x Coerced numeric sample matrix.
#' @returns The `n x d` numeric score matrix.
#' @noRd
#' @keywords internal
.resolve_score <- function(score, x) {
  if (is.null(score)) {
    return(-x)
  }
  if (!is.function(score)) {
    stop("`score` must be `NULL` or a function returning the score matrix.",
         call. = FALSE)
  }

  s <- as.matrix(score(x))
  storage.mode(s) <- "double"
  if (nrow(s) != nrow(x) || ncol(s) != ncol(x)) {
    stop("`score(x)` must return a matrix with the same dimensions as `x` ",
         "(", nrow(x), " x ", ncol(x), "); got ",
         nrow(s), " x ", ncol(s), ".", call. = FALSE)
  }
  if (any(!is.finite(s))) {
    stop("`score(x)` returned non-finite values.", call. = FALSE)
  }
  s
}

#' Resolve the Stein-kernel bandwidth
#'
#' @param bandwidth `"median"` or a positive scalar.
#' @param x Coerced numeric sample matrix.
#' @returns A positive numeric scalar.
#' @noRd
#' @keywords internal
.ksd_bandwidth <- function(bandwidth, x) {
  if (identical(bandwidth, "median")) {
    return(median_bandwidth_cpp(x))
  }
  if (!is.numeric(bandwidth) || length(bandwidth) != 1L || bandwidth <= 0) {
    stop("`bandwidth` must be \"median\" or a single positive number.",
         call. = FALSE)
  }
  bandwidth
}
