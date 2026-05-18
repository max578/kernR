#' Nystrom Low-Rank Kernel Factorisation
#'
#' Builds a Nystrom approximation \eqn{K \approx F F^\top} of an
#' `n x n` kernel matrix using `m << n` uniformly-sampled landmarks.
#' The returned factor `F` is an `n x m` matrix; downstream kernel
#' computations (HSIC, MMD, etc.) that can be expressed in `F` reduce
#' from `O(n^2)` to `O(n m)` cost.
#'
#' The construction is:
#' 1. Sample `m` landmark indices uniformly without replacement.
#' 2. Compute the landmark Gram \eqn{W = K_{mm}} (`m x m`) and the
#'    cross-Gram \eqn{C = K_{nm}} (`n x m`).
#' 3. Stabilise: \eqn{W_\epsilon = W + \epsilon I}.
#' 4. Cholesky factor \eqn{W_\epsilon = L L^\top}.
#' 5. Return \eqn{F = C L^{-\top}}, so that
#'    \eqn{F F^\top = C W_\epsilon^{-1} C^\top \approx K}.
#'
#' Any [kernel_spec()] is supported. Bandwidth selection (median
#' heuristic) is performed on the full dataset before landmark sampling.
#'
#' @param x Numeric matrix `n x d` (or vector). Data points.
#' @param kernel A [kernel_spec()]. Default RBF with median heuristic.
#' @param m Integer. Number of landmark points. Capped at `nrow(x) - 1`.
#'   Default `100L`.
#' @param regularise Small positive numeric. Ridge added to `W` before
#'   Cholesky for numerical stability. Default `1e-6`.
#' @param seed Integer or `NULL`. Random seed for landmark sampling.
#'
#' @return An object of class `"kernel_factor"` with components:
#'   \describe{
#'     \item{F}{The `n x m_eff` factor matrix.}
#'     \item{method}{`"nystrom"`.}
#'     \item{m}{Effective rank (`<= m`).}
#'     \item{kernel}{Resolved kernel spec.}
#'     \item{n}{Number of rows in the input.}
#'   }
#'
#' @references
#' Williams, C. K. I., & Seeger, M. (2001). Using the Nystrom method
#' to speed up kernel machines. *NeurIPS*, 13.
#'
#' @seealso [rff_features()], [hsic_test_nystrom()]
#'
#' @examples
#' set.seed(1)
#' x <- matrix(stats::rnorm(2000L), ncol = 2L)
#' f <- nystrom_factor(x, m = 80L, seed = 1L)
#' dim(f$F)
#'
#' @export
nystrom_factor <- function(x,
                           kernel = kernel_spec(),
                           m = 100L,
                           regularise = 1e-6,
                           seed = NULL) {
  x <- as.matrix(x)
  n <- nrow(x)
  m <- as.integer(m)
  if (length(m) != 1L || is.na(m) || m < 2L) {
    stop("`m` must be an integer >= 2.", call. = FALSE)
  }
  if (m >= n) m <- n - 1L
  if (!is.numeric(regularise) || regularise < 0) {
    stop("`regularise` must be a non-negative numeric.", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)

  kernel <- resolve_bandwidth(kernel, x)
  landmarks <- sort(sample.int(n, m))
  xm <- x[landmarks, , drop = FALSE]

  W <- kernel_matrix(xm, kernel = kernel)
  C <- kernel_matrix(x,  xm, kernel = kernel)

  W_reg <- W + regularise * diag(m)
  L <- tryCatch(
    chol(W_reg),
    error = function(e) {
      # Fall back to larger ridge if numerically indefinite.
      chol(W + max(regularise * 1e3, 1e-3) * diag(m))
    }
  )
  # R's chol() returns upper triangular U with W = U^T U, so
  # W^{-1} = U^{-1} U^{-T}. We want F such that F F^T = C W^{-1} C^T,
  # giving F = C U^{-1}, equivalently F^T = U^{-T} C^T. backsolve()
  # with transpose=TRUE solves U^T y = t(C), so y = U^{-T} C^T = F^T.
  Fmat <- t(backsolve(L, t(C), transpose = TRUE))

  structure(
    list(
      F      = Fmat,
      method = "nystrom",
      m      = m,
      kernel = kernel,
      n      = n
    ),
    class = "kernel_factor"
  )
}


#' Random Fourier Features for the RBF Kernel
#'
#' Builds an `n x D` feature map \eqn{\Phi} such that
#' \eqn{\Phi \Phi^\top \approx K} for the RBF kernel
#' \eqn{K(x, y) = \exp(-\|x-y\|^2 / (2\sigma^2))}, via Rahimi & Recht
#' (2007): draw \eqn{\omega_k \sim N(0, \sigma^{-2} I_d)} and
#' \eqn{b_k \sim U[0, 2\pi]}, then
#' \eqn{\phi_k(x) = \sqrt{2/D} \cos(\omega_k^\top x + b_k)}.
#'
#' Currently RBF-only; other shift-invariant kernels (Matern with
#' specific `nu`) require their own Fourier spectra and are not yet
#' implemented.
#'
#' @param x Numeric matrix `n x d` (or vector).
#' @param kernel A [kernel_spec()] with `type = "rbf"`. Bandwidth may be
#'   `"median"` (resolved against `x`) or a positive numeric.
#' @param D Integer. Number of random features. Larger D -> better
#'   approximation but higher memory / compute. Default `200L`.
#' @param seed Integer or `NULL`. Random seed.
#'
#' @return An object of class `"kernel_factor"` with components:
#'   \describe{
#'     \item{F}{The `n x D` random-feature matrix.}
#'     \item{method}{`"rff"`.}
#'     \item{m}{Equal to `D`.}
#'     \item{kernel}{Resolved kernel spec.}
#'     \item{n}{Number of rows in the input.}
#'     \item{omega, b}{Random draws (kept for reproducible re-encoding).}
#'   }
#'
#' @references
#' Rahimi, A., & Recht, B. (2007). Random features for large-scale
#' kernel machines. *NeurIPS*, 20.
#'
#' @seealso [nystrom_factor()], [hsic_test_nystrom()]
#'
#' @examples
#' set.seed(1)
#' x <- matrix(stats::rnorm(2000L), ncol = 2L)
#' f <- rff_features(x, D = 150L, seed = 1L)
#' dim(f$F)
#'
#' @export
rff_features <- function(x,
                         kernel = kernel_spec("rbf"),
                         D = 200L,
                         seed = NULL) {
  x <- as.matrix(x)
  n <- nrow(x)
  d <- ncol(x)
  D <- as.integer(D)
  if (length(D) != 1L || is.na(D) || D < 2L) {
    stop("`D` must be an integer >= 2.", call. = FALSE)
  }
  if (!inherits(kernel, "kernel_spec") || kernel$type != "rbf") {
    stop("`rff_features()` currently supports only RBF kernels.",
         call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)

  kernel <- resolve_bandwidth(kernel, x)
  sigma  <- kernel$bandwidth

  omega <- matrix(stats::rnorm(D * d, sd = 1 / sigma), nrow = d, ncol = D)
  b     <- stats::runif(D, min = 0, max = 2 * pi)

  proj <- x %*% omega                                # n x D
  Fmat <- sqrt(2 / D) * cos(sweep(proj, 2L, b, "+"))

  structure(
    list(
      F      = Fmat,
      method = "rff",
      m      = D,
      kernel = kernel,
      n      = n,
      omega  = omega,
      b      = b
    ),
    class = "kernel_factor"
  )
}


#' @export
print.kernel_factor <- function(x, ...) {
  cat("\n  Kernel factor (",
      x$method, ")\n\n", sep = "")
  cat("Rows:   ", x$n, "\n")
  cat("Rank:   ", x$m, "\n")
  cat("Kernel: ", x$kernel$type)
  if (x$kernel$type %in% c("rbf", "matern") &&
      is.numeric(x$kernel$bandwidth)) {
    cat(" (bw = ", formatC(x$kernel$bandwidth, digits = 4, format = "g"),
        ")", sep = "")
  }
  cat("\n\n")
  invisible(x)
}


# Internal: HSIC^biased from two centred factors.
# Returns (1/n^2) * ||Fxc^T Fyc||_F^2 where Fxc, Fyc are n x m_x, n x m_y
# kernel factors with column means subtracted.
hsic_from_factors <- function(Fxc, Fyc, n) {
  A <- crossprod(Fxc, Fyc)   # m_x x m_y
  sum(A * A) / (n * n)
}


#' HSIC Independence Test via Low-Rank Factorisation
#'
#' Accelerated HSIC test using either Nystrom (default) or random
#' Fourier features (RFF) factorisations of the input kernel matrices.
#' For large `n` this scales as `O(n m)` per permutation instead of
#' `O(n^2)`, with `m << n` controlling the speed / accuracy trade-off.
#'
#' The test uses the *biased* HSIC estimator
#' \eqn{(1/n^2) \mathrm{tr}(H K_x H K_y)}, which is the form that
#' factorises cleanly through low-rank approximations. The bias is
#' `O(1/n)` and negligible in the large-`n` regime where Nystrom / RFF
#' are useful.
#'
#' The permutation null is built by row-permuting the (centred) `y`
#' factor; per-permutation cost is `O(n m_x m_y)`.
#'
#' @param x,y Numeric matrices (or vectors). Same number of rows.
#' @param kernel_x,kernel_y [kernel_spec()]s for the two factors.
#' @param method Character. `"nystrom"` (default) or `"rff"`.
#' @param m Integer. Rank used for the approximation (number of Nystrom
#'   landmarks or RFF features). Used for both factors unless `m_x` /
#'   `m_y` are supplied.
#' @param m_x,m_y Optional integers overriding `m` per factor.
#' @param n_permutations Integer. Default `500L`.
#' @param alpha Numeric in `(0, 1)`. Default `0.05`.
#' @param seed Integer or `NULL`.
#' @param regularise Ridge for Nystrom Cholesky (ignored under
#'   `method = "rff"`). Default `1e-6`.
#'
#' @return An object of class `"kernel_test_result"` with the standard
#'   fields (`statistic`, `p_value`, `method`, `n`, `n_permutations`,
#'   `null_distribution`, `kernel_x`, `kernel_y`, `call`) plus:
#'   \describe{
#'     \item{approximation}{`"nystrom"` or `"rff"`.}
#'     \item{m_x, m_y}{Effective ranks used.}
#'   }
#'
#' @references
#' Williams, C. K. I., & Seeger, M. (2001). Using the Nystrom method
#' to speed up kernel machines. *NeurIPS*, 13.
#'
#' Rahimi, A., & Recht, B. (2007). Random features for large-scale
#' kernel machines. *NeurIPS*, 20.
#'
#' Gretton, A., Bousquet, O., Smola, A., & Scholkopf, B. (2005).
#' Measuring statistical dependence with Hilbert-Schmidt norms.
#' *ALT*, 63-77.
#'
#' @seealso [hsic_test()], [nystrom_factor()], [rff_features()]
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 1500L
#' x <- stats::rnorm(n)
#' y <- x^2 + stats::rnorm(n, sd = 0.5)
#' fit <- hsic_test_nystrom(x, y, m = 80L,
#'                          n_permutations = 99L, seed = 1L)
#' fit
#' }
#'
#' @export
hsic_test_nystrom <- function(x, y,
                              kernel_x = kernel_spec(),
                              kernel_y = kernel_spec(),
                              method = c("nystrom", "rff"),
                              m = 100L,
                              m_x = NULL,
                              m_y = NULL,
                              n_permutations = 500L,
                              alpha = 0.05,
                              seed = NULL,
                              regularise = 1e-6) {
  cl <- match.call()
  method <- match.arg(method)

  x <- as.matrix(x)
  y <- as.matrix(y)
  n <- nrow(x)
  if (nrow(y) != n) {
    stop("`x` and `y` must have the same number of observations.",
         call. = FALSE)
  }
  if (n < 10L) {
    stop("At least 10 observations are required.", call. = FALSE)
  }
  n_permutations <- as.integer(n_permutations)
  if (length(n_permutations) != 1L || is.na(n_permutations) ||
      n_permutations < 1L) {
    stop("`n_permutations` must be a positive integer.", call. = FALSE)
  }
  if (is.null(m_x)) m_x <- m
  if (is.null(m_y)) m_y <- m

  if (!is.null(seed)) set.seed(seed)

  factorise <- function(z, kernel, m_rank) {
    if (method == "nystrom") {
      nystrom_factor(z, kernel = kernel, m = m_rank,
                     regularise = regularise)
    } else {
      rff_features(z, kernel = kernel, D = m_rank)
    }
  }

  fx <- factorise(x, kernel_x, m_x)
  fy <- factorise(y, kernel_y, m_y)

  # Centre column-wise; centring + row permutation commute.
  Fxc <- sweep(fx$F, 2L, colMeans(fx$F), "-")
  Fyc <- sweep(fy$F, 2L, colMeans(fy$F), "-")

  stat_obs <- hsic_from_factors(Fxc, Fyc, n)

  null_dist <- numeric(n_permutations)
  for (p in seq_len(n_permutations)) {
    perm <- sample.int(n)
    null_dist[p] <- hsic_from_factors(Fxc, Fyc[perm, , drop = FALSE], n)
  }
  p_value <- (1 + sum(null_dist >= stat_obs)) / (1 + n_permutations)

  structure(
    list(
      statistic         = stat_obs,
      p_value           = p_value,
      method            = paste0("HSIC (", method, ")"),
      n                 = n,
      n_permutations    = n_permutations,
      null_distribution = null_dist,
      ess               = NA_real_,
      weights           = NULL,
      kernel_x          = fx$kernel,
      kernel_y          = fy$kernel,
      approximation     = method,
      m_x               = fx$m,
      m_y               = fy$m,
      call              = cl
    ),
    class = "kernel_test_result"
  )
}
