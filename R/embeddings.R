#' Estimate Conditional Mean Embedding via Kernel Ridge Regression
#'
#' Estimates the conditional mean embedding mu(Y|X=x) in the RKHS
#' using kernel ridge regression (Muandet et al., 2017).
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
#'     \item{alpha}{Weight matrix (n x n) for the embedding.}
#'     \item{Ky}{Kernel matrix of `y`.}
#'     \item{x_train}{Training `x` data.}
#'     \item{kernel_x}{Resolved kernel specification.}
#'     \item{kernel_y}{Resolved kernel specification.}
#'     \item{lambda}{Regularisation parameter used.}
#'   }
#'
#' @keywords internal
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

#' Predict Conditional Mean Embedding at New Points
#'
#' @param object A `cme_fit` object.
#' @param x_new Numeric matrix of new conditioning points.
#' @param ... Currently ignored.
#'
#' @return An n_new x n_train matrix of embedding weights. Each row
#'   gives the weights to combine training Y kernel values.
#' @keywords internal
predict.cme_fit <- function(object, x_new, ...) {
  x_new <- as.matrix(x_new)
  # k(x_new, x_train) %*% W
  Kxs <- kernel_matrix(x_new, object$x_train, kernel = object$kernel_x)
  Kxs %*% object$W
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
  # We approximate with trace of squared residuals
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
