#' Aggregate-Likelihood Downscaling
#'
#' Inverts a known aggregation operator `Y = T(X) + eps` to recover the
#' fine-scale latent `X` from coarse / aggregate observations `Y`, using
#' a Gaussian-mixture prior on the latent. Implements the
#' aggregate-likelihood / kernel-downsizing framework (Sejdinovic et
#' al.) as a kernR-side method consuming an optional
#' [proxymix::fit_proxymix()] latent prior.
#'
#' This is the **third** downscaling method in kernR, structurally
#' different from the existing pair:
#'
#' * [kernel_downscale()] (CME, Park-Muandet-Fukumizu-Sejdinovic 2013) --
#'   paired (coarse, fine) training data; supervised regression.
#' * [dist_regression()] (Szabo-Sriperumbudur-Poczos-Gretton 2016) --
#'   distribution-to-distribution mapping from a bag-of-points design.
#' * `aggregate_downscale()` (this function) -- only aggregate
#'   observations + known aggregator + parametric latent prior. Used when
#'   no paired training data exists (only coarse-grid observations) but
#'   the aggregation operator is known (spatial averaging, temporal
#'   averaging, linear or non-linear projection).
#'
#' Two computational paths, selected on `aggregator`'s class:
#'
#' * **Linear-Gaussian closed form** (when `aggregator` is a matrix `A`):
#'   each prior component's posterior is a Kalman update;
#'   `K_k = Sigma_k A^T (A Sigma_k A^T + sigma_y^2 I)^{-1}`;
#'   `mu_k|y = mu_k + K_k (y - A mu_k)`,
#'   `Sigma_k|y = (I - K_k A) Sigma_k`;
#'   posterior mixture weights reweight by per-component evidence
#'   `N(y | A mu_k, A Sigma_k A^T + sigma_y^2 I)`.
#' * **Non-linear importance sampling** (when `aggregator` is a
#'   function): draw `n_samples_per_component` samples from each prior
#'   component, evaluate `T(.)`, weight by Gaussian likelihood
#'   `N(y | T(x), sigma_y^2 I)`, recover posterior moments + reweighted
#'   mixture weights from the importance-weighted samples. Reports
#'   per-component effective sample size; warns when below a stated
#'   floor.
#'
#' @param y Numeric vector of length `dim_y`, or `1 x dim_y` matrix.
#'   The observed aggregate.
#' @param aggregator Either a numeric `dim_y x dim_x` matrix (treated
#'   as a linear aggregator and dispatched to the closed-form path) or
#'   a function `function(x) -> y_matrix` that maps an `n x dim_x`
#'   matrix of latent samples to an `n x dim_y` matrix of aggregates
#'   (non-linear IS path).
#' @param latent_prior Either (a) a list with elements `means`
#'   (`N x dim_x` matrix or list of `dim_x`-vectors), `covariances`
#'   (list of `N` `dim_x x dim_x` matrices), `weights` (length-`N`
#'   numeric, summing to 1) -- or (b) a `proxymix::gmm_fit` (any object
#'   exposing `@means`, `@covariances`, `@weights` slots).
#' @param sigma_y Numeric. Observation noise standard deviation
#'   (scalar; `eps ~ N(0, sigma_y^2 I)`). Default `0.1`.
#' @param n_samples_per_component Integer. Importance-sampling sample
#'   count per prior component (non-linear path only). Default `200L`.
#' @param min_ess_fraction Numeric in `(0, 1]`. ESS-floor reliability
#'   gate for the IS path: when per-component ESS drops below
#'   `min_ess_fraction * n_samples_per_component`, a warning is
#'   emitted. Default `0.1`. Set `0` to disable.
#' @param seed Integer or `NULL`. Random seed (non-linear path).
#'
#' @return An object of class `"aggregate_downscale"` with components:
#'   * `posterior_mean` -- length-`dim_x` numeric, `E[X | y]`.
#'   * `posterior_cov` -- `dim_x x dim_x` matrix, `Cov[X | y]`
#'     (law-of-total-covariance over mixture components).
#'   * `posterior_weights` -- length-`N` numeric, posterior mixture
#'     weights (sum to 1).
#'   * `posterior_components_means` -- list of `N` length-`dim_x`
#'     posterior component means.
#'   * `posterior_components_covariances` -- list of `N` posterior
#'     component covariances.
#'   * `aggregator_type` -- `"linear"` or `"nonlinear"`.
#'   * `method` -- `"linear_closed_form"` or `"nonlinear_is"`.
#'   * `ess_per_component` -- length-`N` per-component ESS (IS path
#'     only; `NA` for closed form).
#'   * `ess_warning` -- `TRUE` if any per-component ESS fell below the
#'     floor (IS path only).
#'   * `n_components`, `sigma_y`, `n_samples_per_component`, `call`.
#'
#' @references
#' Sejdinovic, D. (talk; 2025). *Kernel downsizing and aggregate
#' likelihoods.* Companion talk to Hoek-Elliott (2024) for the
#' aggregate-likelihood / GMM-proxy direction. The companion proxymix
#' Tier-2 stub `proxymix::from_aggregate_likelihood()` targets the
#' same problem from the prior-fitting side; this function targets it
#' from the consumption side (inversion given a fitted prior).
#'
#' @seealso [kernel_downscale()], [dist_regression()].
#' @examples
#' set.seed(1L)
#' # Linear-Gaussian: spatial averaging of two adjacent cells.
#' A <- matrix(c(0.5, 0.5), nrow = 1L)
#' prior <- list(
#'   means = list(c(0, 0), c(2, 2)),
#'   covariances = list(diag(2L), diag(2L)),
#'   weights = c(0.5, 0.5)
#' )
#' fit <- aggregate_downscale(y = 1.0, aggregator = A,
#'                            latent_prior = prior, sigma_y = 0.2)
#' fit$posterior_mean
#'
#' # Non-linear: aggregator is sin of the sum.
#' agg_fn <- function(x) matrix(sin(rowSums(x)), ncol = 1L)
#' fit2 <- aggregate_downscale(y = 0.5, aggregator = agg_fn,
#'                             latent_prior = prior, sigma_y = 0.1,
#'                             n_samples_per_component = 300L, seed = 1L)
#' fit2$posterior_mean
#' @export
aggregate_downscale <- function(y, aggregator, latent_prior,
                                sigma_y = 0.1,
                                n_samples_per_component = 200L,
                                min_ess_fraction = 0.1,
                                seed = NULL) {
  cl <- match.call()
  if (length(sigma_y) != 1L || !is.finite(sigma_y) || sigma_y <= 0) {
    stop("`sigma_y` must be a single positive number.", call. = FALSE)
  }
  if (length(min_ess_fraction) != 1L || !is.finite(min_ess_fraction) ||
      min_ess_fraction < 0 || min_ess_fraction > 1) {
    stop("`min_ess_fraction` must be in [0, 1].", call. = FALSE)
  }
  prior <- .normalise_aggregate_prior(latent_prior)
  N <- length(prior$weights)
  dim_x <- length(prior$means[[1L]])

  y_vec <- as.numeric(y)
  dim_y <- length(y_vec)

  if (is.matrix(aggregator) && !is.function(aggregator)) {
    .aggregate_downscale_linear(y_vec, aggregator, prior,
                                sigma_y, dim_x, dim_y, N, cl)
  } else if (is.function(aggregator)) {
    if (!is.null(seed)) set.seed(seed)
    .aggregate_downscale_nonlinear(y_vec, aggregator, prior,
                                   sigma_y, n_samples_per_component,
                                   min_ess_fraction,
                                   dim_x, dim_y, N, cl)
  } else {
    stop("`aggregator` must be a matrix (linear) or a function ",
         "(non-linear).", call. = FALSE)
  }
}

# ---- internals -------------------------------------------------------

.normalise_aggregate_prior <- function(p) {
  # Accept either a bare list or any object exposing means/
  # covariances/weights slots (e.g. proxymix::gmm_fit).
  if (!is.list(p) ||
      !all(c("means", "covariances", "weights") %in% names(p))) {
    p_slot <- tryCatch(
      list(means       = methods::slot(p, "means"),
           covariances = methods::slot(p, "covariances"),
           weights     = methods::slot(p, "weights")),
      error = function(e) NULL
    )
    if (is.null(p_slot)) {
      stop("`latent_prior` must be a list with elements means / ",
           "covariances / weights, or an object with those slots ",
           "(e.g., a proxymix::gmm_fit).", call. = FALSE)
    }
    p <- p_slot
  }
  means <- p$means
  if (is.matrix(means)) {
    means <- lapply(seq_len(nrow(means)), function(k) means[k, ])
  } else if (!is.list(means)) {
    stop("`means` must be a matrix or list of vectors.", call. = FALSE)
  }
  covs <- p$covariances
  if (!is.list(covs)) {
    stop("`covariances` must be a list of matrices.", call. = FALSE)
  }
  if (length(means) != length(covs)) {
    stop("`means` and `covariances` must have the same length.",
         call. = FALSE)
  }
  w <- as.numeric(p$weights)
  if (length(w) != length(means)) {
    stop("`weights` length must match the number of components.",
         call. = FALSE)
  }
  if (any(!is.finite(w)) || any(w < 0) || sum(w) <= 0) {
    stop("`weights` must be non-negative, finite, and sum > 0.",
         call. = FALSE)
  }
  list(means = means, covariances = covs, weights = w / sum(w))
}

.aggregate_downscale_linear <- function(y_vec, A, prior, sigma_y,
                                        dim_x, dim_y, N, cl) {
  if (nrow(A) != dim_y || ncol(A) != dim_x) {
    stop("aggregator matrix has dims ", nrow(A), "x", ncol(A),
         " but expected ", dim_y, "x", dim_x,
         " (rows = dim_y, cols = dim_x).", call. = FALSE)
  }
  Sigma_y_full <- sigma_y^2 * diag(dim_y)
  log_evidence <- numeric(N)
  post_means <- vector("list", N)
  post_covs  <- vector("list", N)
  for (k in seq_len(N)) {
    mu_k    <- prior$means[[k]]
    Sigma_k <- prior$covariances[[k]]
    S_yk    <- A %*% Sigma_k %*% t(A) + Sigma_y_full
    # Solve via Cholesky for numerical stability.
    L_yk <- tryCatch(chol(S_yk),
                     error = function(e) chol(S_yk +
                       1e-10 * diag(dim_y)))
    resid <- y_vec - as.numeric(A %*% mu_k)
    # log N(y | A mu_k, S_yk)
    half_logdet <- sum(log(diag(L_yk)))
    quad <- sum(backsolve(L_yk, resid, transpose = TRUE)^2)
    log_evidence[k] <- -0.5 * dim_y * log(2 * pi) -
                       half_logdet - 0.5 * quad
    # Kalman gain & per-component posterior
    Kgain <- Sigma_k %*% t(A) %*%
      chol2inv(L_yk)
    post_means[[k]] <- as.numeric(mu_k + Kgain %*% resid)
    post_covs[[k]]  <- Sigma_k - Kgain %*% A %*% Sigma_k
    # Symmetrise to defend against ~ULP asymmetry
    post_covs[[k]]  <- 0.5 * (post_covs[[k]] + t(post_covs[[k]]))
  }
  log_w_post <- log(prior$weights) + log_evidence
  log_w_post <- log_w_post - max(log_w_post)
  w_post <- exp(log_w_post)
  w_post <- w_post / sum(w_post)

  post_mean <- Reduce(`+`,
    lapply(seq_len(N), function(k) w_post[k] * post_means[[k]]))
  post_cov <- Reduce(`+`,
    lapply(seq_len(N), function(k) {
      diff <- post_means[[k]] - post_mean
      w_post[k] * (post_covs[[k]] + diff %*% t(diff))
    }))

  structure(
    list(
      posterior_mean                   = as.numeric(post_mean),
      posterior_cov                    = post_cov,
      posterior_weights                = w_post,
      posterior_components_means       = post_means,
      posterior_components_covariances = post_covs,
      aggregator_type                  = "linear",
      method                           = "linear_closed_form",
      ess_per_component                = rep(NA_real_, N),
      ess_warning                      = FALSE,
      n_components                     = N,
      sigma_y                          = sigma_y,
      n_samples_per_component          = NA_integer_,
      call                             = cl
    ),
    class = "aggregate_downscale"
  )
}

.aggregate_downscale_nonlinear <- function(y_vec, agg_fn, prior,
                                           sigma_y,
                                           n_samples_per_component,
                                           min_ess_fraction,
                                           dim_x, dim_y, N, cl) {
  M <- as.integer(n_samples_per_component)
  if (M < 1L) {
    stop("`n_samples_per_component` must be a positive integer.",
         call. = FALSE)
  }
  log_evidence <- numeric(N)
  post_means   <- vector("list", N)
  post_covs    <- vector("list", N)
  ess_k        <- numeric(N)
  any_ess_low  <- FALSE

  for (k in seq_len(N)) {
    mu_k    <- prior$means[[k]]
    Sigma_k <- prior$covariances[[k]]
    L_k <- tryCatch(chol(Sigma_k),
                    error = function(e) chol(Sigma_k +
                      1e-10 * diag(dim_x)))
    Z <- matrix(stats::rnorm(M * dim_x), M, dim_x)
    X_k <- Z %*% L_k + matrix(mu_k, M, dim_x, byrow = TRUE)
    Y_pred <- agg_fn(X_k)
    if (!is.matrix(Y_pred)) Y_pred <- matrix(Y_pred, ncol = dim_y)
    if (nrow(Y_pred) != M || ncol(Y_pred) != dim_y) {
      stop("aggregator function must return an n x dim_y matrix; ",
           "got ", nrow(Y_pred), "x", ncol(Y_pred), ".",
           call. = FALSE)
    }
    resid <- sweep(Y_pred, 2L, y_vec, FUN = "-")
    log_lik <- -0.5 * dim_y * log(2 * pi) - dim_y * log(sigma_y) -
               0.5 * rowSums(resid^2) / sigma_y^2
    # MC estimate of log p(y | k) = log( mean exp(log_lik) )
    lmax <- max(log_lik)
    log_evidence[k] <- lmax + log(mean(exp(log_lik - lmax)))
    # Importance weights within component (normalised)
    iw <- exp(log_lik - lmax)
    iw <- iw / sum(iw)
    ess_k[k] <- 1 / sum(iw^2)
    if (ess_k[k] < min_ess_fraction * M) any_ess_low <- TRUE
    mu_post_k <- colSums(iw * X_k)
    # Weighted sample covariance about the weighted mean
    Xc <- sweep(X_k, 2L, mu_post_k, FUN = "-")
    Sigma_post_k <- crossprod(Xc * sqrt(iw))
    Sigma_post_k <- 0.5 * (Sigma_post_k + t(Sigma_post_k))
    post_means[[k]] <- as.numeric(mu_post_k)
    post_covs[[k]]  <- Sigma_post_k
  }

  log_w_post <- log(prior$weights) + log_evidence
  log_w_post <- log_w_post - max(log_w_post)
  w_post <- exp(log_w_post)
  w_post <- w_post / sum(w_post)

  post_mean <- Reduce(`+`,
    lapply(seq_len(N), function(k) w_post[k] * post_means[[k]]))
  post_cov <- Reduce(`+`,
    lapply(seq_len(N), function(k) {
      diff <- post_means[[k]] - post_mean
      w_post[k] * (post_covs[[k]] + diff %*% t(diff))
    }))

  if (any_ess_low) {
    low_idx <- which(ess_k < min_ess_fraction * M)
    warning(
      "aggregate_downscale(): per-component IS effective sample size ",
      "fell below ",
      formatC(100 * min_ess_fraction, digits = 0, format = "f"),
      "% of n_samples_per_component for component(s) ",
      paste(low_idx, collapse = ", "),
      " (ESS = ",
      paste(formatC(ess_k[low_idx], digits = 1, format = "f"),
            collapse = ", "),
      "). Increase `n_samples_per_component`, sharpen the prior, or ",
      "consider a coarser aggregator if posterior is locally narrow.",
      call. = FALSE
    )
  }

  structure(
    list(
      posterior_mean                   = as.numeric(post_mean),
      posterior_cov                    = post_cov,
      posterior_weights                = w_post,
      posterior_components_means       = post_means,
      posterior_components_covariances = post_covs,
      aggregator_type                  = "nonlinear",
      method                           = "nonlinear_is",
      ess_per_component                = ess_k,
      ess_warning                      = any_ess_low,
      n_components                     = N,
      sigma_y                          = sigma_y,
      n_samples_per_component          = M,
      call                             = cl
    ),
    class = "aggregate_downscale"
  )
}

#' @export
print.aggregate_downscale <- function(x, digits = 3L, ...) {
  cat("Aggregate-likelihood downscaling\n")
  cat("  method:        ", x$method, " (", x$aggregator_type, ")\n",
      sep = "")
  cat("  components:    ", x$n_components, "\n", sep = "")
  cat("  sigma_y:       ", formatC(x$sigma_y, digits = digits,
                                   format = "g"), "\n", sep = "")
  if (identical(x$method, "nonlinear_is")) {
    cat("  IS samples/k:  ", x$n_samples_per_component, "\n",
        sep = "")
    cat("  ESS per k:     ",
        paste(formatC(x$ess_per_component, digits = 1, format = "f"),
              collapse = ", "), "\n", sep = "")
    if (isTRUE(x$ess_warning)) {
      cat("  WARNING:       per-component ESS below floor\n")
    }
  }
  cat("  posterior mean:", paste(formatC(x$posterior_mean,
                                         digits = digits,
                                         format = "g"),
                                 collapse = ", "), "\n", sep = " ")
  cat("  posterior wts: ",
      paste(formatC(x$posterior_weights, digits = digits,
                    format = "f"), collapse = ", "),
      "\n", sep = "")
  invisible(x)
}

#' Sample from the posterior of an aggregate-downscale fit
#'
#' Draws `n` samples from the per-component posterior mixture.
#' Each draw picks a component by `posterior_weights`, then samples
#' from `N(posterior_components_means[[k]], posterior_components_covariances[[k]])`.
#'
#' @param object An `aggregate_downscale` fit.
#' @param n Integer. Number of posterior samples. Default `1000L`.
#' @param seed Integer or `NULL`. Random seed.
#' @return An `n x dim_x` numeric matrix.
#' @export
posterior_sample_aggregate <- function(object, n = 1000L,
                                       seed = NULL) {
  if (!inherits(object, "aggregate_downscale")) {
    stop("`object` must be an `aggregate_downscale` fit.",
         call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)
  n <- as.integer(n)
  N <- object$n_components
  dim_x <- length(object$posterior_mean)
  k_draws <- sample.int(N, n, replace = TRUE,
                        prob = object$posterior_weights)
  out <- matrix(NA_real_, n, dim_x)
  for (k in seq_len(N)) {
    rows <- which(k_draws == k)
    if (!length(rows)) next
    mu_k <- object$posterior_components_means[[k]]
    Sigma_k <- object$posterior_components_covariances[[k]]
    L_k <- tryCatch(chol(Sigma_k),
                    error = function(e) chol(Sigma_k +
                      1e-10 * diag(dim_x)))
    Z <- matrix(stats::rnorm(length(rows) * dim_x), length(rows),
                dim_x)
    out[rows, ] <- Z %*% L_k +
                   matrix(mu_k, length(rows), dim_x, byrow = TRUE)
  }
  out
}
