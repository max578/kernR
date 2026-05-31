#' Estimate Conditional Mean Embedding via Kernel Ridge Regression
#'
#' Estimates the conditional mean embedding \eqn{\mu_{Y|X=x}} in the
#' RKHS using kernel ridge regression (Park, Muandet, Fukumizu &
#' Sejdinovic, 2013; Muandet et al., 2017).
#'
#' This is the lower-level building block used by [kernel_downscale()].
#' Most users should call `kernel_downscale()`; use `fit_cme()` directly
#' when you need access to the trained operator (the weight matrix `W`)
#' for custom downstream computations -- e.g. constructing a kernel
#' Bayes' rule update, plugging into a manuscript figure pipeline, or
#' composing with other RKHS operators.
#'
#' @param x Numeric matrix of conditioning variables (n x d_x).
#' @param y Numeric matrix of target variables (n x d_y).
#' @param kernel_x Kernel specification for `x`.
#' @param kernel_y Kernel specification for `y`.
#' @param lambda Ridge regularisation parameter. If `"cv"`, selected
#'   by leave-one-out cross-validation. Default is `"cv"`.
#'
#' @return A list of class `"cme_fit"` with components:
#'   \describe{
#'     \item{W}{Operator matrix `(K_x + n lambda I)^{-1}` (n x n).}
#'     \item{Ky}{Kernel matrix of `y`.}
#'     \item{x_train}{Training `x` data.}
#'     \item{kernel_x, kernel_y}{Resolved kernel specifications.}
#'     \item{lambda}{Regularisation parameter used.}
#'   }
#'
#' @references
#' Park, J., Muandet, K., Fukumizu, K., & Sejdinovic, D. (2013).
#' *Kernel embeddings of conditional distributions: A unified kernel
#' framework for nonparametric inference in graphical models.*
#' IEEE Signal Processing Magazine.
#'
#' Muandet, K., Fukumizu, K., Sriperumbudur, B., & Scholkopf, B.
#' (2017). *Kernel mean embedding of distributions: A review and
#' beyond.* Foundations and Trends in Machine Learning, 10(1-2).
#'
#' @seealso [predict.cme_fit()], [kernel_downscale()]
#' @examples
#' set.seed(1L)
#' x <- matrix(rnorm(60L), ncol = 2L)
#' y <- matrix(x[, 1L] + rnorm(30L, sd = 0.2), ncol = 1L)
#' fit <- fit_cme(x, y, lambda = 1e-2)
#' dim(fit$W)
#'
#' @family downscaling and embeddings
#' @export
fit_cme <- function(x, y,
                    kernel_x = kernel_spec(),
                    kernel_y = kernel_spec(),
                    lambda = "cv") {
  x <- as.matrix(x)
  y <- as.matrix(y)
  n <- nrow(x)

  kernel_x <- resolve_bandwidth(kernel_x, x)
  kernel_y <- resolve_bandwidth(kernel_y, y)

  Kx <- kernel_matrix(x, kernel = kernel_x)
  Ky <- kernel_matrix(y, kernel = kernel_y)

  # Select lambda by LOO-CV if requested
  if (identical(lambda, "cv")) {
    lambda <- cv_ridge_lambda(Kx, Ky)
  }

  # Solve (Kx + n*lambda*I)^{-1}
  # W = (Kx + n*lambda*I)^{-1}, alpha = W %*% Ky
  reg <- Kx + n * lambda * diag(n)
  W <- solve(reg)

  structure(
    list(
      W = W,
      Ky = Ky,
      x_train = x,
      kernel_x = kernel_x,
      kernel_y = kernel_y,
      lambda = lambda
    ),
    class = "cme_fit"
  )
}

#' Predict Conditional Mean Embedding Weights at New Points
#'
#' Returns the row weights \eqn{\alpha(x^*) = k(x^*, X_{\mathrm{train}})
#' W} from a fitted [fit_cme()] object. Each row of the result is the
#' weight vector for combining training-`y` quantities (kernel values
#' or values themselves) to produce a CME prediction at `x_new`.
#'
#' For a typical "predict Y at new X" workflow use [kernel_downscale()],
#' which combines this with the training Y matrix to return predictions
#' directly.
#'
#' @param object A `cme_fit` object.
#' @param x_new Numeric matrix of new conditioning points.
#' @param ... Currently ignored.
#'
#' @return An n_new x n_train matrix of embedding weights.
#' @export
predict.cme_fit <- function(object, x_new, ...) {
  x_new <- as.matrix(x_new)
  Kxs <- kernel_matrix(x_new, object$x_train, kernel = object$kernel_x)
  Kxs %*% object$W
}

#' Print a Conditional Mean Embedding Fit
#'
#' Prints a compact summary of a fitted conditional mean embedding: the
#' training-sample size, the input and output dimensions, the resolved
#' kernels, and the ridge regularisation parameter that was used.
#'
#' @param x A `cme_fit` object returned by [fit_cme()].
#' @param ... Unused; present for S3 generic compatibility.
#'
#' @return The `cme_fit` object `x`, invisibly.
#'
#' @seealso [fit_cme()]
#' @examples
#' set.seed(1L)
#' x <- matrix(rnorm(60L), ncol = 2L)
#' y <- matrix(x[, 1L] + rnorm(30L, sd = 0.2), ncol = 1L)
#' print(fit_cme(x, y, lambda = 1e-2))
#'
#' @export
print.cme_fit <- function(x, ...) {
  cat("Conditional mean embedding (kernel ridge regression)\n")
  cat("  Training points: ", nrow(x$x_train), "\n", sep = "")
  cat("  Input dim:       ", ncol(x$x_train), "\n", sep = "")
  cat("  Output dim:      ", ncol(x$Ky), "\n", sep = "")
  cat("  Kernel (x):      ", x$kernel_x$type, "\n", sep = "")
  cat("  Kernel (y):      ", x$kernel_y$type, "\n", sep = "")
  cat("  Ridge lambda:    ",
    formatC(x$lambda, digits = 4, format = "g"), "\n",
    sep = ""
  )
  invisible(x)
}

#' Cross-Validate Ridge Parameter for KRR
#'
#' Simple LOO-CV for the ridge parameter in kernel ridge regression.
#' Tests a grid of lambda values and picks the one minimising LOO error.
#'
#' @param Kx n x n kernel matrix.
#' @param Ky n x n kernel matrix.
#'
#' @return Optimal lambda (numeric scalar).
#' @keywords internal
cv_ridge_lambda <- function(Kx, Ky) {
  n <- nrow(Kx)
  lambdas <- 10^seq(-6, 1, length.out = 15)

  # LOO error via the hat matrix: H = K(K + n*lam*I)^{-1}
  # LOO residual for obs i: (Ky - H %*% Ky)[i,] / (1 - H[i,i])
  # Exact LOO-CV mean squared error via the hat-matrix shortcut
  errors <- vapply(lambdas, function(lam) {
    reg <- Kx + n * lam * diag(n)
    H <- Kx %*% solve(reg)
    resid <- Ky - H %*% Ky
    diag_h <- pmax(1 - diag(H), 1e-10)
    # Mean squared LOO error (summing over kernel dimensions)
    mean((resid / diag_h)^2)
  }, numeric(1))

  lambdas[which.min(errors)]
}
