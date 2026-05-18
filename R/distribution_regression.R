# ---------------------------------------------------------------------------
# Internal helpers (not exported)
# ---------------------------------------------------------------------------

# Validate a list-of-matrices bag input.
.validate_bags <- function(bags, name = "bags", min_bags = 5L) {
  if (!is.list(bags) || length(bags) < min_bags) {
    stop(sprintf("`%s` must be a list of at least %d matrices.",
                 name, min_bags),
         call. = FALSE)
  }
  bags <- lapply(bags, as.matrix)
  d <- ncol(bags[[1L]])
  for (i in seq_along(bags)) {
    if (ncol(bags[[i]]) != d) {
      stop("All bags must have the same number of columns.",
           call. = FALSE)
    }
    if (nrow(bags[[i]]) < 1L) {
      stop(sprintf("Bag %d is empty.", i), call. = FALSE)
    }
  }
  bags
}

# Compute the bag-inner-mean Gram block G[i, j] = mean(k(bag_i_pts, bag_j_pts))
# for a pair of bag lists. When `bags_b` is NULL, computes the symmetric
# training Gram. Resolves median bandwidth on pooled `bags_a` data.
.bag_gram <- function(bags_a, bags_b = NULL,
                      inner_kernel = kernel_spec()) {
  symmetric <- is.null(bags_b)
  if (symmetric) bags_b <- bags_a
  M <- length(bags_a); N <- length(bags_b)

  if (inner_kernel$type %in% c("rbf", "matern") &&
      identical(inner_kernel$bandwidth, "median")) {
    pooled <- do.call(rbind, bags_a)
    inner_kernel <- resolve_bandwidth(inner_kernel, pooled)
  }

  G <- matrix(NA_real_, M, N)
  for (i in seq_len(M)) {
    j_start <- if (symmetric) i else 1L
    for (j in j_start:N) {
      Kij <- kernel_matrix(bags_a[[i]], bags_b[[j]],
                           kernel = inner_kernel)
      G[i, j] <- mean(Kij)
      if (symmetric && i != j) G[j, i] <- G[i, j]
    }
  }
  list(G = G, inner_kernel = inner_kernel)
}

# LOO-CV ridge for kernel ridge regression with scalar / matrix y.
.cv_ridge_lambda_y <- function(K, y) {
  n <- nrow(K)
  y <- as.matrix(y)
  lambdas <- 10^seq(-6, 1, length.out = 15L)
  errors <- vapply(lambdas, function(lam) {
    H <- K %*% solve(K + n * lam * diag(n))
    resid <- y - H %*% y
    diag_h <- pmax(1 - diag(H), 1e-10)
    mean((resid / diag_h)^2)
  }, numeric(1L))
  lambdas[which.min(errors)]
}


# ---------------------------------------------------------------------------
# User-facing API
# ---------------------------------------------------------------------------

#' Kernel Distribution Regression
#'
#' Regression where each input is a *bag* of points (a sample from a
#' distribution) rather than a single feature vector. Each bag is
#' implicitly mapped to its empirical mean embedding in the RKHS of an
#' inner kernel; the outer (between-bag) kernel acts on those
#' embeddings, and kernel ridge regression predicts a scalar or
#' multivariate output. This is the Szabó-Sriperumbudur-Póczos-Gretton
#' (2016) "learning theory for distribution regression" setup;
#' kernel-mean-embedding background is from Muandet et al. (2017).
#'
#' **Outer kernels supported.**
#' - `"linear"`: \eqn{K(P_i, P_j) = \langle \hat{\mu}_{P_i},
#'   \hat{\mu}_{P_j} \rangle =
#'   \frac{1}{n_i n_j} \sum_{k,l} k(x_i^{(k)}, x_j^{(l)})}.
#' - `"rbf"`: \eqn{K(P_i, P_j) = \exp(-\| \hat{\mu}_{P_i} -
#'   \hat{\mu}_{P_j} \|^2 / (2\sigma^2))}, where the embedding-space
#'   distance is recovered from the inner Gram via
#'   \eqn{\|\hat\mu_i - \hat\mu_j\|^2 = G_{ii} - 2 G_{ij} + G_{jj}}.
#'
#' Typical ag-systems use: each paddock contributes a *bag* of soil-core
#' measurements (variable depth, multiple cores per paddock); the
#' regression predicts paddock-level yield from the distributional
#' shape of the soil profile. Distinct from [kernel_downscale()] in
#' that the input is itself a distribution, not a fixed-length vector.
#'
#' @param bags A list of length `M` of numeric matrices. Each
#'   element is a bag (`n_i x d` points). All bags must have the same
#'   number of columns; bag sizes `n_i` may differ.
#' @param y Numeric vector of length `M`, or numeric matrix `M x d_y`,
#'   of training targets.
#' @param inner_kernel A [kernel_spec()] applied between points within
#'   and across bags. Median bandwidth heuristic resolved on the
#'   pooled training points. Default RBF with median heuristic.
#' @param outer Character. Outer-kernel form: `"linear"` (default) or
#'   `"rbf"`.
#' @param outer_bandwidth Outer-kernel bandwidth (RBF only). `"median"`
#'   (default) resolves to the median embedding-space pairwise distance
#'   on the training Gram; a positive numeric overrides.
#' @param lambda Ridge regularisation for kernel ridge regression. If
#'   `"cv"` (default), selected by leave-one-out CV over
#'   `10^seq(-6, 1, length.out = 15)`.
#'
#' @return An object of class `"dist_regression"` with components:
#'   \describe{
#'     \item{alpha}{Ridge weights (length `M` or `M x d_y`).}
#'     \item{G_train}{Training bag-inner-mean Gram `(M x M)`.}
#'     \item{K_outer_train}{Outer kernel matrix `(M x M)` actually used.}
#'     \item{bags_train}{Bags used (kept for prediction).}
#'     \item{y_train}{Targets used.}
#'     \item{inner_kernel}{Resolved inner kernel.}
#'     \item{outer, outer_bandwidth}{Outer kernel choice and resolved
#'       bandwidth (`NA` for linear outer).}
#'     \item{lambda}{Ridge parameter used.}
#'     \item{M, d_y, call}{Metadata.}
#'   }
#'
#' @references
#' Szabó, Z., Sriperumbudur, B. K., Póczos, B., & Gretton, A. (2016).
#' Learning theory for distribution regression.
#' *Journal of Machine Learning Research*, 17(152), 1-40.
#'
#' Muandet, K., Fukumizu, K., Sriperumbudur, B., & Scholkopf, B.
#' (2017). *Kernel mean embedding of distributions: A review and
#' beyond.* Foundations and Trends in Machine Learning, 10(1-2).
#'
#' @seealso [predict.dist_regression()], [kernel_downscale()]
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' # 30 bags, each a sample from N(mu_i, 1); predict mu_i
#' M <- 30L
#' mu <- stats::rnorm(M)
#' bags <- lapply(mu, function(m)
#'   matrix(stats::rnorm(40, mean = m), ncol = 1L))
#' fit <- dist_regression(bags, y = mu, outer = "linear")
#' fit
#'
#' # Predict at new bags
#' new_mu <- stats::rnorm(5L)
#' new_bags <- lapply(new_mu, function(m)
#'   matrix(stats::rnorm(40, mean = m), ncol = 1L))
#' predict(fit, new_bags)
#' }
#'
#' @export
dist_regression <- function(bags, y,
                            inner_kernel = kernel_spec(),
                            outer = c("linear", "rbf"),
                            outer_bandwidth = "median",
                            lambda = "cv") {
  cl <- match.call()
  outer <- match.arg(outer)

  bags <- .validate_bags(bags)
  M <- length(bags)
  if (is.numeric(y) && is.null(dim(y))) {
    if (length(y) != M) {
      stop("`y` length must equal the number of bags.", call. = FALSE)
    }
    y_mat <- matrix(y, ncol = 1L)
  } else {
    y_mat <- as.matrix(y)
    if (nrow(y_mat) != M) {
      stop("`y` must have one row per bag.", call. = FALSE)
    }
  }

  # Inner Gram on training bags (resolves bandwidth on pooled points).
  g <- .bag_gram(bags, inner_kernel = inner_kernel)
  G <- g$G
  inner_kernel <- g$inner_kernel

  if (outer == "linear") {
    K_outer <- G
    sigma_outer <- NA_real_
  } else {  # rbf
    D2 <- outer(diag(G), diag(G), "+") - 2 * G
    D2[D2 < 0] <- 0  # numerical floor
    if (identical(outer_bandwidth, "median")) {
      D <- sqrt(D2)
      lower <- D[lower.tri(D)]
      lower <- lower[lower > 0]
      sigma_outer <- if (length(lower) > 0L) stats::median(lower) else 1
      if (sigma_outer <= 0) sigma_outer <- 1
    } else {
      if (!is.numeric(outer_bandwidth) || outer_bandwidth <= 0) {
        stop("`outer_bandwidth` must be \"median\" or a positive number.",
             call. = FALSE)
      }
      sigma_outer <- outer_bandwidth
    }
    K_outer <- exp(-D2 / (2 * sigma_outer^2))
  }

  if (identical(lambda, "cv")) {
    lambda <- .cv_ridge_lambda_y(K_outer, y_mat)
  }
  if (!is.numeric(lambda) || lambda < 0) {
    stop("`lambda` must be \"cv\" or a non-negative number.",
         call. = FALSE)
  }

  alpha <- solve(K_outer + M * lambda * diag(M), y_mat)
  if (ncol(y_mat) == 1L && is.null(dim(y))) alpha <- drop(alpha)

  out <- list(
    alpha           = alpha,
    G_train         = G,
    K_outer_train   = K_outer,
    bags_train      = bags,
    y_train         = y_mat,
    inner_kernel    = inner_kernel,
    outer           = outer,
    outer_bandwidth = sigma_outer,
    lambda          = lambda,
    M               = M,
    d_y             = ncol(y_mat),
    call            = cl
  )
  structure(out, class = "dist_regression")
}


#' Predict from a Fitted Distribution Regression Model
#'
#' @param object A `dist_regression` fit.
#' @param newdata A list of new bags (matrices with the same number of
#'   columns as the training bags). Bag sizes may differ from training.
#' @param ... Currently ignored.
#'
#' @return Numeric vector (length `length(newdata)`) when the model was
#'   fitted with a scalar `y`; numeric matrix (`length(newdata) x d_y`)
#'   for multivariate `y`.
#' @export
predict.dist_regression <- function(object, newdata, ...) {
  new_bags <- .validate_bags(newdata, name = "newdata", min_bags = 1L)

  if (ncol(new_bags[[1L]]) != ncol(object$bags_train[[1L]])) {
    stop("`newdata` bags must have the same number of columns as training bags.",
         call. = FALSE)
  }

  # Cross Gram: rows = new bags, cols = training bags.
  M_new   <- length(new_bags)
  M_train <- object$M
  G_cross <- matrix(NA_real_, M_new, M_train)
  for (i in seq_len(M_new)) {
    for (j in seq_len(M_train)) {
      Kij <- kernel_matrix(new_bags[[i]], object$bags_train[[j]],
                           kernel = object$inner_kernel)
      G_cross[i, j] <- mean(Kij)
    }
  }

  if (object$outer == "linear") {
    K_cross <- G_cross
  } else {
    # Diagonals of new bags
    G_new_diag <- vapply(new_bags, function(b) {
      mean(kernel_matrix(b, b, kernel = object$inner_kernel))
    }, numeric(1L))
    G_train_diag <- diag(object$G_train)
    D2 <- outer(G_new_diag, G_train_diag, "+") - 2 * G_cross
    D2[D2 < 0] <- 0
    K_cross <- exp(-D2 / (2 * object$outer_bandwidth^2))
  }

  pred <- K_cross %*% as.matrix(object$alpha)
  if (object$d_y == 1L && is.null(dim(object$alpha))) {
    pred <- drop(pred)
  } else if (object$d_y > 1L) {
    colnames(pred) <- colnames(object$y_train)
  }
  pred
}


#' @export
print.dist_regression <- function(x, ...) {
  cat("\n  Kernel Distribution Regression\n\n")
  cat("Training bags:    ", x$M, "\n")
  bag_sizes <- vapply(x$bags_train, nrow, integer(1L))
  cat("Bag sizes:        ",
      min(bag_sizes), "-", max(bag_sizes),
      " (median ", stats::median(bag_sizes), ")\n", sep = "")
  cat("Point dim:        ", ncol(x$bags_train[[1L]]), "\n")
  cat("Output dim:       ", x$d_y, "\n")
  cat("Inner kernel:     ", x$inner_kernel$type)
  if (is.numeric(x$inner_kernel$bandwidth)) {
    cat(" (bw = ",
        formatC(x$inner_kernel$bandwidth, digits = 4L, format = "g"),
        ")", sep = "")
  }
  cat("\n")
  cat("Outer kernel:     ", x$outer)
  if (x$outer == "rbf") {
    cat(" (bw = ",
        formatC(x$outer_bandwidth, digits = 4L, format = "g"),
        ")", sep = "")
  }
  cat("\n")
  cat("Ridge lambda:     ",
      formatC(x$lambda, digits = 4L, format = "g"), "\n\n")
  invisible(x)
}


#' @export
as.data.frame.dist_regression <- function(x, row.names = NULL,
                                          optional = FALSE, ...) {
  y_hat <- as.numeric(x$K_outer_train %*% as.matrix(x$alpha))
  data.frame(
    bag_id        = seq_len(x$M),
    bag_size      = vapply(x$bags_train, nrow, integer(1L)),
    y_train       = as.numeric(x$y_train[, 1L]),
    y_fit         = if (x$d_y == 1L) y_hat else NA_real_,
    stringsAsFactors = FALSE
  )
}
