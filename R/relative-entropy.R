# Relative-entropy causal metric + objective causal-influence-range integral.
# The kernR-owned statistic of Assimilative Causal Inference (ACI; Andreou, Chen &
# Bollt, Nat. Commun. 17:1854, 2026): a cause is identified at time t when the
# smoother posterior of the (hidden) cause departs from its filter posterior, i.e.
# when incorporating the future of the observed effect reduces uncertainty about
# the cause. kernR supplies the distributional measure (relative entropy) and the
# threshold-free range integral; the filter/smoother engine lives upstream
# (kalmix / PESTO), keeping responsibilities split (orchestra invariant 8).

.as_cov <- function(s, k) {
  if (is.matrix(s)) return(s)
  if (length(s) == 1L) return(diag(as.numeric(s), k))
  if (length(s) == k) return(diag(as.numeric(s), k))
  stop("covariance must be a k x k matrix, a length-k variance vector, or a scalar")
}

#' Relative Entropy Between Two Gaussians (the ACI causal metric)
#'
#' Computes the Kullback-Leibler divergence
#' \eqn{\mathcal{P}(p, q) = \int p \log(p / q)} between two multivariate Gaussian
#' distributions \eqn{p = \mathcal{N}(\mu_p, \Sigma_p)} and
#' \eqn{q = \mathcal{N}(\mu_q, \Sigma_q)}. This is the operational statistic of
#' Assimilative Causal Inference (ACI): with `p` the **smoother** posterior of a
#' hidden cause (using the future of the observed effect) and `q` the **filter**
#' posterior (past only), a non-zero value identifies the hidden variable as a
#' cause of the observed effect at that instant.
#'
#' The closed form for `k`-dimensional Gaussians is
#' \deqn{\tfrac{1}{2}\left[\operatorname{tr}(\Sigma_q^{-1}\Sigma_p)
#'   + (\mu_q-\mu_p)^\top \Sigma_q^{-1} (\mu_q-\mu_p) - k
#'   + \log\frac{\det \Sigma_q}{\det \Sigma_p}\right].}
#' The measure is non-negative, zero if and only if the two Gaussians coincide,
#' and asymmetric in `p` and `q` (it is a divergence, not a distance). It captures
#' differences in both mean and covariance and is invariant under any smooth
#' invertible reparameterisation of the state.
#'
#' Covariance arguments accept a `k x k` matrix, a length-`k` variance vector
#' (diagonal covariance), or a scalar (isotropic). The reference covariance
#' \eqn{\Sigma_q} must be positive definite.
#'
#' @param mu_p Numeric vector. Mean of `p` (the smoother posterior, for ACI).
#' @param sigma_p Covariance of `p`: matrix, variance vector, or scalar.
#' @param mu_q Numeric vector. Mean of `q` (the filter posterior, for ACI).
#' @param sigma_q Covariance of `q`: matrix, variance vector, or scalar. Must be
#'   positive definite (it is inverted).
#'
#' @return A single non-negative numeric: the relative entropy \eqn{\mathcal{P}(p, q)}.
#'
#' @references Andreou, M., Chen, N. & Bollt, E. (2026). Assimilative causal
#'   inference. *Nature Communications* 17, 1854.
#'
#' @examples
#' relative_entropy(0, 1, 0, 1)            # identical -> 0
#' relative_entropy(1, 1, 0, 1)            # mean shift of 1 sd -> 0.5
#' relative_entropy(c(0, 0), diag(2), c(1, 0), diag(2))
#' @export
relative_entropy <- function(mu_p, sigma_p, mu_q, sigma_q) {
  mu_p <- as.numeric(mu_p); mu_q <- as.numeric(mu_q)
  k <- length(mu_p)
  if (length(mu_q) != k) stop("mu_p and mu_q must have the same length")
  Sp <- .as_cov(sigma_p, k); Sq <- .as_cov(sigma_q, k)
  if (any(dim(Sp) != k) || any(dim(Sq) != k))
    stop("covariance dimensions must match the mean length")
  Sq_chol <- tryCatch(chol(Sq), error = function(e)
    stop("sigma_q must be positive definite"))
  Sq_inv <- chol2inv(Sq_chol)
  d <- mu_q - mu_p
  tr_term <- sum(Sq_inv * Sp)                       # tr(Sq^-1 Sp)
  quad <- as.numeric(crossprod(d, Sq_inv %*% d))
  log_det_q <- 2 * sum(log(diag(Sq_chol)))
  log_det_p <- as.numeric(determinant(Sp, logarithm = TRUE)$modulus)
  val <- 0.5 * (tr_term + quad - k + (log_det_q - log_det_p))
  max(val, 0)                                       # guard floating-point noise
}

#' Relative Entropy Between Two Sample Ensembles (Gaussian approximation)
#'
#' Convenience wrapper around [relative_entropy()] for ensemble inputs (for
#' example a PESTO posterior ensemble): each ensemble is summarised by its sample
#' mean and covariance and the Gaussian relative entropy is returned. This is a
#' moment-matched approximation -- exact when the ensembles are Gaussian and a
#' second-order surrogate otherwise.
#'
#' @param p_draws Numeric matrix `n_p x k` (or vector for `k = 1`): draws from `p`
#'   (the smoother ensemble, for ACI).
#' @param q_draws Numeric matrix `n_q x k` (or vector): draws from `q` (the filter
#'   ensemble). Must have the same number of columns as `p_draws`.
#'
#' @return A single non-negative numeric: the moment-matched relative entropy.
#'
#' @seealso [relative_entropy()]
#' @examples
#' set.seed(1)
#' p <- matrix(rnorm(2000), ncol = 2)
#' q <- matrix(rnorm(2000), ncol = 2) + 1
#' relative_entropy_ensemble(p, q)
#' @export
relative_entropy_ensemble <- function(p_draws, q_draws) {
  if (is.null(dim(p_draws))) p_draws <- matrix(p_draws, ncol = 1L)
  if (is.null(dim(q_draws))) q_draws <- matrix(q_draws, ncol = 1L)
  if (ncol(p_draws) != ncol(q_draws))
    stop("p_draws and q_draws must have the same number of columns")
  k <- ncol(p_draws)
  cov_p <- if (nrow(p_draws) > 1L) stats::cov(p_draws) else diag(0, k)
  cov_q <- stats::cov(q_draws)
  relative_entropy(colMeans(p_draws), cov_p, colMeans(q_draws), cov_q)
}

#' Objective Causal Influence Range (threshold-free)
#'
#' Reduces a divergence-versus-lag profile to a single, threshold-free causal
#' influence range (CIR) -- the effective horizon over which a cause keeps
#' informing the estimate of its effect (Andreou et al. 2026, eqs. 8-9). Given the
#' relative entropy `divergence` between the **lagged** smoother (future of the
#' effect included only up to a lag `L`) and the **complete** smoother (all future
#' included), the divergence falls from \eqn{M} at lag 0 (the filter, no future)
#' towards 0 as `L` grows. The subjective CIR at tolerance \eqn{\varepsilon} is
#' \eqn{\tau_\varepsilon = \inf\{L : D(L) \le \varepsilon\}}; integrating it out,
#' \deqn{\tau = \frac{1}{M}\int_0^M \tau_\varepsilon \, d\varepsilon
#'   = \frac{1}{M}\int_0^{L_{\max}} D(L)\, dL,}
#' the second equality by parts. The result is a decorrelation-time analogue: a
#' lead time in the units of `lag`, free of any cut-off threshold.
#'
#' The producer of the `divergence` profile (the lagged/complete smoother passes)
#' is the upstream filtering engine; this function owns only the range integral.
#'
#' @param lag Numeric vector of non-negative, increasing lags `L` (units of time).
#'   Should start at (or include) 0.
#' @param divergence Numeric vector, same length as `lag`: the relative entropy
#'   \eqn{D(L)} of the lag-`L` smoother from the complete smoother. The value at
#'   the smallest lag is taken as \eqn{M} (the normalising maximum).
#'
#' @return A single non-negative numeric: the objective CIR (a lead time in the
#'   units of `lag`). Zero when there is no recoverable future information
#'   (\eqn{M = 0}).
#'
#' @references Andreou, M., Chen, N. & Bollt, E. (2026). Assimilative causal
#'   inference. *Nature Communications* 17, 1854.
#'
#' @examples
#' lag <- seq(0, 10, by = 0.5)
#' D <- 2 * exp(-0.6 * lag)            # divergence decays with lag
#' cir_objective(lag, D)
#' @export
cir_objective <- function(lag, divergence) {
  lag <- as.numeric(lag); divergence <- as.numeric(divergence)
  if (length(lag) != length(divergence) || length(lag) < 2L)
    stop("lag and divergence must be vectors of equal length >= 2")
  if (is.unsorted(lag)) stop("lag must be increasing")
  M <- divergence[1L]
  if (!is.finite(M) || M <= 0) return(0)
  n <- length(lag)
  area <- sum(diff(lag) * (divergence[-n] + divergence[-1L]) / 2)   # trapezoid
  max(area / M, 0)
}
