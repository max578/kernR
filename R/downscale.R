#' Kernel-Based Statistical Downscaling
#'
#' Predicts fine-resolution outputs at new coarse-resolution inputs via
#' conditional mean embedding (CME) regression in an RKHS. Given paired
#' coarse-fine training data `(coarse, fine)`, fits the operator
#' \eqn{E[Y \mid X = x]} in closed form via kernel ridge regression
#' and returns predictions at `new_coarse`. This is the
#' Park-Muandet-Fukumizu-Sejdinovic conditional-mean-embedding scheme,
#' specialised to the regression form needed for spatial / temporal
#' downscaling.
#'
#' Typical ag-systems use: coarse climate-grid inputs (e.g. monthly
#' temperature, rainfall on a 25 km grid) -> fine-resolution outputs
#' (paddock yield, biomass) at the same time index. Train on years
#' where both coarse and fine are observed; predict fine outputs at
#' new coarse inputs.
#'
#' Compared to a linear regression baseline, the kernel approach
#' captures non-linear coarse-fine relationships without specifying
#' the functional form. Compared to deep-learning downscalers, it has
#' a closed-form solution, uses orders-of-magnitude less data, and
#' carries an interpretable kernel-bandwidth degrees-of-freedom knob.
#'
#' @param coarse Numeric matrix `n x d_coarse` of training
#'   coarse-resolution inputs. Vectors are coerced via [as.matrix()].
#' @param fine Numeric matrix `n x d_fine` of training fine-resolution
#'   outputs. Multivariate outputs (`d_fine > 1`) are supported and
#'   predicted jointly.
#' @param new_coarse Numeric matrix `n_new x d_coarse` of coarse inputs
#'   at which to predict the fine outputs.
#' @param kernel_coarse,kernel_fine [kernel_spec()] for the coarse and
#'   fine spaces. Defaults to RBF with median heuristic. `kernel_fine`
#'   is used only for the bandwidth-CV step; predictions are returned
#'   in the original `fine` units.
#' @param lambda Ridge regularisation parameter for the CME ridge
#'   regression. If `"cv"` (default), selected by leave-one-out
#'   cross-validation over `10^seq(-6, 1, length.out = 15)`.
#' @param return_weights Logical. If `TRUE`, the result carries the
#'   `n_new x n_train` weight matrix used to combine training fine
#'   values. Default `FALSE` (saves memory for large designs).
#'
#' @return An object of class `"kernel_downscale"` with components:
#'   \describe{
#'     \item{prediction}{`n_new x d_fine` matrix of predicted fine
#'       outputs at `new_coarse`.}
#'     \item{n_train}{Number of training pairs used.}
#'     \item{n_new}{Number of prediction points.}
#'     \item{lambda}{Regularisation used (CV-selected when
#'       `lambda = "cv"`).}
#'     \item{kernel_coarse, kernel_fine}{Resolved kernel specs.}
#'     \item{weights}{Optional `n_new x n_train` weight matrix.}
#'     \item{call}{The matched call.}
#'   }
#'
#' @references
#' Park, J., Muandet, K., Fukumizu, K., & Sejdinovic, D. (2013).
#' *Kernel embeddings of conditional distributions.* IEEE Signal
#' Processing Magazine.
#'
#' Muandet, K., Fukumizu, K., Sriperumbudur, B., & Scholkopf, B.
#' (2017). *Kernel mean embedding of distributions: A review and
#' beyond.* Foundations and Trends in Machine Learning, 10(1-2).
#'
#' @seealso [fit_cme()], [dist_regression()] (for downscaling when
#'   each coarse "input" is a bag of points rather than a single
#'   vector).
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 80L
#' coarse <- matrix(stats::rnorm(n * 2L), n, 2L,
#'                  dimnames = list(NULL, c("temp", "rainfall")))
#' fine <- cbind(
#'   yield   = 2 * coarse[, "rainfall"] -
#'             0.5 * coarse[, "temp"]^2 +
#'             stats::rnorm(n, sd = 0.2),
#'   biomass = coarse[, "temp"] + coarse[, "rainfall"]^2 +
#'             stats::rnorm(n, sd = 0.2)
#' )
#' # Predict at a held-out grid
#' new_coarse <- matrix(stats::rnorm(20L * 2L), 20L, 2L,
#'                      dimnames = list(NULL, c("temp", "rainfall")))
#' fit <- kernel_downscale(coarse, fine, new_coarse)
#' print(fit)
#' head(fit$prediction)
#' }
#'
#' @family downscaling and embeddings
#' @export
kernel_downscale <- function(coarse, fine, new_coarse,
                             kernel_coarse = kernel_spec(),
                             kernel_fine   = kernel_spec(),
                             lambda = "cv",
                             return_weights = FALSE) {
  cl <- match.call()

  coarse     <- as.matrix(coarse)
  fine       <- as.matrix(fine)
  new_coarse <- as.matrix(new_coarse)

  n <- nrow(coarse)
  if (nrow(fine) != n) {
    stop("`coarse` and `fine` must have the same number of rows.",
         call. = FALSE)
  }
  if (n < 10L) {
    stop("At least 10 training pairs are required.", call. = FALSE)
  }
  if (ncol(new_coarse) != ncol(coarse)) {
    stop("`new_coarse` must have the same number of columns as `coarse`.",
         call. = FALSE)
  }

  cme <- fit_cme(coarse, fine,
                 kernel_x = kernel_coarse,
                 kernel_y = kernel_fine,
                 lambda   = lambda)

  weights <- predict(cme, new_coarse)        # n_new x n_train
  prediction <- weights %*% fine             # n_new x d_fine
  colnames(prediction) <- colnames(fine)

  out <- list(
    prediction     = prediction,
    n_train        = n,
    n_new          = nrow(new_coarse),
    lambda         = cme$lambda,
    kernel_coarse  = cme$kernel_x,
    kernel_fine    = cme$kernel_y,
    weights        = if (isTRUE(return_weights)) weights else NULL,
    call           = cl
  )
  structure(out, class = "kernel_downscale")
}


#' @export
print.kernel_downscale <- function(x, ...) {
  cat("\n  Kernel Downscaling (CME)\n\n")
  cat("Training pairs:  ", x$n_train, "\n")
  cat("Prediction points:", x$n_new, "\n")
  cat("Output dims:     ", ncol(x$prediction), "\n")
  cat("Kernel (coarse): ", x$kernel_coarse$type)
  if (is.numeric(x$kernel_coarse$bandwidth)) {
    cat(" (bw = ",
        formatC(x$kernel_coarse$bandwidth, digits = 4L, format = "g"),
        ")", sep = "")
  }
  cat("\n")
  cat("Lambda (ridge):  ",
      formatC(x$lambda, digits = 4L, format = "g"), "\n\n")
  invisible(x)
}


#' @export
as.data.frame.kernel_downscale <- function(x, row.names = NULL,
                                           optional = FALSE, ...) {
  pred <- x$prediction
  out <- as.data.frame(pred, optional = TRUE)
  if (!is.null(row.names)) rownames(out) <- row.names
  out
}
