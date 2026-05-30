#' Fit a Density-Ratio Model
#'
#' Trains a density-ratio estimator for the do-null reweighting
#' `w(x, z) = p*(x) / p(x | z)` used by [bd_hsic_test()]. The fitted
#' model is decoupled from evaluation: [predict_density_ratio()] applies
#' it to held-out rows, so train/test splits are honoured cleanly.
#'
#' Four backends are supported (the `method` argument):
#'
#' * `"logistic"` (default), `"ranger"`, `"xgboost"` -- classifier-based
#'   noise-contrastive estimation. The classifier is trained to
#'   distinguish joint samples `(x, z)` from product-of-marginals
#'   samples `(x_perm, z)`; the density ratio is recovered from the
#'   calibrated class probabilities. Log-ratios are stored internally
#'   for numerical stability.
#' * `"proxymix"` -- Gaussian-mixture density-ratio. Fits one GMM to the
#'   joint sample cloud `(x, z)` and one to a permuted
#'   product-of-marginals cloud via `proxymix::fit_proxymix(regime =
#'   "sample")`; ratios are evaluated in log-space from
#'   `proxymix::dgmm()`. Per-GMM convergence diagnostics (BIC, AIC,
#'   final log-likelihood, iteration count) are surfaced on the
#'   returned fit; query them via `fit$diagnostics`.
#'
#' Introduced in kernR 0.0.0.9014 to close the documented-but-
#' unimplemented sample-split gap in [bd_hsic_test()] (see NEWS).
#' [estimate_density_ratio()] is now a thin backwards-compatible
#' wrapper that fits and predicts on the same data.
#'
#' @param x Numeric vector or matrix. Treatment variable (training).
#' @param z Numeric matrix or data.frame. Confounders (training).
#' @param method Character. Backend: `"logistic"` (default),
#'   `"ranger"`, `"xgboost"`, or `"proxymix"`.
#' @param n_noise Integer. Noise samples per real sample for
#'   classifier backends. Default `1L`.
#' @param proxymix_components Integer. Mixture components per density
#'   when `method = "proxymix"`. Default `2L`.
#' @param seed Integer or `NULL`. Random seed.
#'
#' @return An object of class `density_ratio_fit` (plus
#'   `density_ratio_fit_<method>` as the dispatch class). Carries:
#'   `method`, the backend-specific fit (`model` for classifiers;
#'   `fit_joint` + `fit_marg` for proxymix), `diagnostics`,
#'   `n_train`, `ncol_x`, `ncol_z`, `seed`.
#' @seealso [predict_density_ratio()], [estimate_density_ratio()],
#'   [bd_hsic_test()].
#' @examples
#' set.seed(1L)
#' n <- 200L
#' z <- matrix(rnorm(n * 2L), n, 2L)
#' x <- z[, 1L] + rnorm(n)
#' fit <- fit_density_ratio(x, z, method = "logistic", seed = 1L)
#' fit$diagnostics
#'
#' @export
fit_density_ratio <- function(x, z,
                              method = c("logistic", "ranger",
                                         "xgboost", "proxymix"),
                              n_noise = 1L,
                              proxymix_components = 2L,
                              seed = NULL) {
  method <- match.arg(method)
  x <- as.matrix(x)
  z <- as.matrix(z)
  n <- nrow(x)
  if (nrow(z) != n) {
    stop("`x` and `z` must have the same number of rows.",
         call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)

  if (method == "proxymix") {
    return(.fit_density_ratio_proxymix(x, z, proxymix_components, seed))
  }

  n_neg <- n * n_noise
  perm_idx <- sample.int(n, n_neg, replace = TRUE)
  x_perm <- x[perm_idx, , drop = FALSE]
  z_rep <- z[rep(seq_len(n), n_noise), , drop = FALSE]
  features_pos <- cbind(x, z)
  features_neg <- cbind(x_perm, z_rep)
  features <- rbind(features_pos, features_neg)
  labels <- c(rep(1L, n), rep(0L, n_neg))
  colnames(features) <- c(paste0("x", seq_len(ncol(x))),
                          paste0("z", seq_len(ncol(z))))

  model <- switch(method,
    logistic = {
      df <- data.frame(y = labels, features)
      stats::glm(y ~ ., data = df, family = stats::binomial())
    },
    ranger = {
      if (!requireNamespace("ranger", quietly = TRUE)) {
        stop("Package 'ranger' required for method = 'ranger'.",
             call. = FALSE)
      }
      df <- data.frame(y = factor(labels), features)
      ranger::ranger(y ~ ., data = df, probability = TRUE,
                     num.trees = 500L)
    },
    xgboost = {
      if (!requireNamespace("xgboost", quietly = TRUE)) {
        stop("Package 'xgboost' required for method = 'xgboost'.",
             call. = FALSE)
      }
      dtrain <- xgboost::xgb.DMatrix(data = features, label = labels)
      xgboost::xgb.train(
        params = list(objective = "binary:logistic",
                      max_depth = 4L, eta = 0.1),
        data = dtrain, nrounds = 100L, verbose = 0L
      )
    }
  )

  structure(
    list(
      method        = method,
      model         = model,
      n_train       = n,
      n_noise       = as.integer(n_noise),
      ncol_x        = ncol(x),
      ncol_z        = ncol(z),
      feature_names = colnames(features),
      diagnostics   = list(method = method, n_train = n,
                           n_noise = as.integer(n_noise)),
      seed          = seed
    ),
    class = c(paste0("density_ratio_fit_", method),
              "density_ratio_fit")
  )
}

.fit_density_ratio_proxymix <- function(x, z, n_components, seed) {
  if (!requireNamespace("proxymix", quietly = TRUE)) {
    stop("method = \"proxymix\" requires the `proxymix` package ",
         "(>= 0.3.0).", call. = FALSE)
  }
  n <- nrow(x)
  joint <- cbind(x, z)
  perm  <- sample.int(n)
  marg  <- cbind(x[perm, , drop = FALSE], z)
  target_joint <- proxymix::gmm_target_from_samples(joint)
  target_marg  <- proxymix::gmm_target_from_samples(marg)
  fit_joint <- proxymix::fit_proxymix(target_joint,
                                      N = as.integer(n_components),
                                      regime = "sample")
  fit_marg  <- proxymix::fit_proxymix(target_marg,
                                      N = as.integer(n_components),
                                      regime = "sample")
  diagnostics <- list(
    method            = "proxymix",
    n_components      = as.integer(n_components),
    joint_converged   = fit_joint@converged,
    marg_converged    = fit_marg@converged,
    joint_loglik      = fit_joint@diagnostics$loglik_final,
    marg_loglik       = fit_marg@diagnostics$loglik_final,
    joint_bic         = fit_joint@diagnostics$bic,
    marg_bic          = fit_marg@diagnostics$bic,
    joint_aic         = fit_joint@diagnostics$aic,
    marg_aic          = fit_marg@diagnostics$aic,
    joint_iterations  = fit_joint@iterations,
    marg_iterations   = fit_marg@iterations
  )
  structure(
    list(
      method       = "proxymix",
      fit_joint    = fit_joint,
      fit_marg     = fit_marg,
      n_train      = n,
      n_components = as.integer(n_components),
      ncol_x       = ncol(x),
      ncol_z       = ncol(z),
      diagnostics  = diagnostics,
      seed         = seed
    ),
    class = c("density_ratio_fit_proxymix", "density_ratio_fit")
  )
}

#' Predict from a Fitted Density-Ratio Model
#'
#' Applies a `density_ratio_fit` object (from [fit_density_ratio()])
#' to new `(x, z)` rows. All four backends compute ratios in log-space
#' internally for numerical stability; `type` controls the returned
#' representation.
#'
#' @param object A `density_ratio_fit`.
#' @param new_x Numeric vector or matrix. Treatment values to evaluate.
#' @param new_z Numeric matrix or data.frame. Confounders to evaluate.
#' @param type Character. Return type:
#'   `"log_ratio"` -- natural-log density ratio (default; preferred for
#'   downstream calculation);
#'   `"ratio"` -- raw density ratio (`exp(log_ratio)`);
#'   `"weight"` -- IES-compatible normalised weights (positive,
#'   sum-to-`n_new`).
#'
#' @return Numeric vector of length `nrow(new_x)`.
#' @seealso [fit_density_ratio()].
#' @examples
#' set.seed(1L)
#' n <- 200L
#' z <- matrix(rnorm(n * 2L), n, 2L)
#' x <- z[, 1L] + rnorm(n)
#' fit <- fit_density_ratio(x, z, method = "logistic", seed = 1L)
#' weights <- predict_density_ratio(fit, new_x = x, new_z = z,
#'                                  type = "weight")
#' summary(weights)
#'
#' @export
predict_density_ratio <- function(object, new_x, new_z,
                                  type = c("log_ratio", "weight",
                                           "ratio")) {
  if (!inherits(object, "density_ratio_fit")) {
    stop("`object` must be a `density_ratio_fit` from ",
         "fit_density_ratio().", call. = FALSE)
  }
  type <- match.arg(type)
  new_x <- as.matrix(new_x)
  new_z <- as.matrix(new_z)
  n_new <- nrow(new_x)
  if (nrow(new_z) != n_new) {
    stop("`new_x` and `new_z` must have the same number of rows.",
         call. = FALSE)
  }
  if (ncol(new_x) != object$ncol_x ||
      ncol(new_z) != object$ncol_z) {
    stop("`new_x` / `new_z` column counts do not match the training ",
         "data (training had ncol_x = ", object$ncol_x,
         ", ncol_z = ", object$ncol_z, ").", call. = FALSE)
  }

  if (object$method == "proxymix") {
    joint_new <- cbind(new_x, new_z)
    p_joint <- proxymix::dgmm(joint_new, object$fit_joint)
    p_marg  <- proxymix::dgmm(joint_new, object$fit_marg)
    # Floor in log-space to avoid -Inf on extreme tails.
    log_ratio <- log(pmax(p_joint, 1e-300)) -
                 log(pmax(p_marg,  1e-300))
  } else {
    features_new <- cbind(new_x, new_z)
    colnames(features_new) <- object$feature_names
    probs <- switch(object$method,
      logistic = stats::predict(object$model,
                                newdata = data.frame(features_new),
                                type = "response"),
      ranger   = stats::predict(object$model,
                                data = data.frame(features_new))$
                   predictions[, "1"],
      xgboost  = stats::predict(object$model, features_new)
    )
    probs <- pmax(pmin(as.numeric(probs), 1 - 1e-6), 1e-6)
    # log[p / (1 - p)] + log(n_noise) is the calibrated log-ratio
    # under noise-contrastive estimation.
    log_ratio <- log(probs) - log(1 - probs) +
                 log(object$n_noise)
  }

  switch(type,
    log_ratio = log_ratio,
    ratio     = exp(log_ratio),
    weight    = {
      ratios <- exp(log_ratio)
      ratios <- pmax(ratios, 1e-8)
      ratios / sum(ratios) * n_new
    }
  )
}

#' @export
print.density_ratio_fit <- function(x, ...) {
  cat("Density-ratio fit (", x$method, ")\n", sep = "")
  cat("  n_train: ", x$n_train, "\n", sep = "")
  cat("  ncol_x:  ", x$ncol_x,  "\n", sep = "")
  cat("  ncol_z:  ", x$ncol_z,  "\n", sep = "")
  if (identical(x$method, "proxymix")) {
    cat("  components: ", x$n_components,
        "  joint_converged: ", x$diagnostics$joint_converged,
        "  marg_converged: ",  x$diagnostics$marg_converged, "\n",
        sep = "")
    cat("  joint BIC: ",
        formatC(x$diagnostics$joint_bic, digits = 3, format = "g"),
        "  marg BIC: ",
        formatC(x$diagnostics$marg_bic,  digits = 3, format = "g"),
        "\n", sep = "")
  }
  invisible(x)
}

#' Estimate Density Ratios (backwards-compatible wrapper)
#'
#' Wraps [fit_density_ratio()] + [predict_density_ratio()] on the same
#' data. Preserved for backwards compatibility with kernR 0.0.0.901x
#' callers; new code should prefer the explicit fit/predict pair so
#' train/test splits are honoured.
#'
#' The return shape (`weights`, `ratios`, `ess`, `method`, `n`) is
#' unchanged from previous versions. Internally, ratios are now
#' computed in log-space (which fixes pathological tail behaviour
#' that the classifier-based 0.0.0.9012 implementation occasionally
#' showed under extreme imbalance).
#'
#' @inheritParams fit_density_ratio
#' @return A list of class `density_ratio_fit_estimate` with components
#'   `weights`, `ratios`, `ess`, `method`, `n`, and `fit` (the
#'   underlying `density_ratio_fit` for callers that want diagnostics).
#' @seealso [fit_density_ratio()] for the fit/predict surface.
#' @examples
#' set.seed(42)
#' n <- 200
#' z <- matrix(rnorm(n * 2), n, 2)
#' x <- z[, 1] + rnorm(n)
#' dr <- estimate_density_ratio(x, z)
#' dr$ess
#' @export
estimate_density_ratio <- function(x, z,
                                   method = c("logistic", "ranger",
                                              "xgboost", "proxymix"),
                                   n_noise = 1L,
                                   proxymix_components = 2L,
                                   seed = NULL) {
  method <- match.arg(method)
  fit <- fit_density_ratio(x, z, method = method, n_noise = n_noise,
                           proxymix_components = proxymix_components,
                           seed = seed)
  weights <- predict_density_ratio(fit, x, z, type = "weight")
  ratios  <- predict_density_ratio(fit, x, z, type = "ratio")
  ess <- effective_sample_size(weights)
  structure(
    list(
      weights = as.numeric(weights),
      ratios  = as.numeric(ratios),
      ess     = ess,
      method  = method,
      n       = nrow(as.matrix(x)),
      fit     = fit
    ),
    # Inherits "density_ratio_fit" so legacy call-sites that do
    # `inherits(dr, "density_ratio_fit")` keep working after the
    # 0.0.0.9014 fit/predict refactor.
    class = c("density_ratio_fit_estimate", "density_ratio_fit_legacy",
              "density_ratio_fit")
  )
}

#' @export
print.density_ratio_fit_legacy <- function(x, ...) {
  cat("Density ratio estimation (", x$method, ")\n", sep = "")
  cat("  N:   ", x$n, "\n")
  cat("  ESS: ", formatC(x$ess, digits = 1, format = "f"), "\n")
  cat("  Weight range: [",
    formatC(min(x$weights), digits = 3, format = "g"), ", ",
    formatC(max(x$weights), digits = 3, format = "g"), "]\n",
    sep = ""
  )
  invisible(x)
}


#' Estimate Density Ratios via RuLSIF
#'
#' Relative unconstrained Least-Squares Importance Fitting.
#' Kernel-based closed-form density ratio estimation.
#'
#' @param x_num Numeric matrix. Numerator samples.
#' @param x_den Numeric matrix. Denominator samples.
#' @param kernel Kernel specification. Default is RBF with median heuristic.
#' @param lambda Regularisation parameter. Default is 0.1.
#' @param alpha Relative parameter (0 = LSIF, 0.5 = symmetric). Default is 0.
#'
#' @return Named list with `weights` and `ess`.
#'
#' @keywords internal
estimate_rulsif <- function(x_num, x_den,
                            kernel = kernel_spec(),
                            lambda = 0.1,
                            alpha = 0) {
  x_num <- as.matrix(x_num)
  x_den <- as.matrix(x_den)

  n_num <- nrow(x_num)
  n_den <- nrow(x_den)

  kernel <- resolve_bandwidth(kernel, rbind(x_num, x_den))

  K_num_num <- kernel_matrix(x_num, kernel = kernel)
  K_den_num <- kernel_matrix(x_den, x_num, kernel = kernel)

  H <- (1 - alpha) / n_den * crossprod(K_den_num) +
    alpha / n_num * K_num_num +
    lambda * diag(n_num)
  h <- colMeans(K_num_num)

  theta <- solve(H, h)
  theta <- pmax(theta, 0)

  weights <- as.numeric(K_den_num %*% theta)
  weights <- pmax(weights, 1e-8)
  weights <- weights / sum(weights) * n_den

  ess <- effective_sample_size(weights)

  list(weights = weights, ess = ess, method = "rulsif")
}
