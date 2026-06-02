# score-helpers.R -- score-function adapters for ksd_test().
#
# ksd_test() needs the score of the target (the gradient of its log density).
# gaussian_score() (in test-ksd.R) covers the multivariate-normal target in
# closed form. numeric_score() covers everything else: it turns any log-density
# function into a score by central finite differences, so a target need only be
# expressible as a (possibly unnormalised) log density, not hand-differentiated.
# This is also the kernR-side adapter for the constellation C4 edge: once the
# leader exports a log-posterior generic, wrapping it in numeric_score() feeds
# ksd_test() directly -- no sibling dependency is added (the user supplies a
# plain function), so the edge stays dormant-but-correct until C4 ships.

#' Finite-difference score from a log-density
#'
#' Builds a score function -- the gradient of the log density -- by central
#' finite differences, for use as the `score` argument of [ksd_test()]. The
#' target need only be expressible as a function returning a (possibly
#' unnormalised) log density; any additive normalising constant cancels in the
#' gradient, so the constant may be omitted.
#'
#' This makes [ksd_test()] usable against any target a user can write down as a
#' log density, not only targets with a hand-derived score. It is also the
#' kernR-side adapter for a log-posterior contract: wrapping an exported
#' log-posterior evaluator in `numeric_score()` yields the score
#' [ksd_test()] needs to check whether posterior draws are calibrated against
#' that posterior, with no dependency added on the producer -- the evaluator is
#' passed in as an ordinary function.
#'
#' Central differences cost `2 * d` evaluations of `log_density` per call, where
#' `d` is the dimension. For a cheap closed-form target prefer [gaussian_score()]
#' or a hand-written score; `numeric_score()` is the general fallback.
#'
#' @param log_density A function accepting an `n x d` numeric matrix and
#'   returning a numeric vector of length `n`: the log density (up to an
#'   additive constant) evaluated row-wise.
#' @param h Positive numeric. Finite-difference step. Default `1e-4`. Too large
#'   biases the gradient; too small amplifies floating-point noise.
#'
#' @returns A function of one argument (an `n x d` numeric matrix) returning the
#'   `n x d` matrix of finite-difference scores, suitable as the `score`
#'   argument of [ksd_test()].
#'
#' @examples
#' set.seed(1)
#' x <- matrix(stats::rnorm(200L), ncol = 2L)
#'
#' # Standard-normal log density (unnormalised): -||x||^2 / 2
#' ld <- function(z) -0.5 * rowSums(z^2)
#' s <- numeric_score(ld)
#'
#' # Matches the closed-form standard-normal score -x to finite-difference order
#' max(abs(s(x) - (-x)))
#'
#' # Use directly in a goodness-of-fit test
#' ksd_test(x, score = numeric_score(ld), n_boot = 199L, seed = 1L)
#'
#' @seealso [ksd_test()], [gaussian_score()]
#' @family goodness-of-fit tests
#' @author Max Moldovan, \email{max.moldovan@@adelaide.edu.au}
#' @export
numeric_score <- function(log_density, h = 1e-4) {
  if (!is.function(log_density)) {
    stop("`log_density` must be a function of an `n x d` matrix.",
         call. = FALSE)
  }
  if (!is.numeric(h) || length(h) != 1L || h <= 0) {
    stop("`h` must be a single positive number.", call. = FALSE)
  }
  force(log_density)
  force(h)

  function(x) {
    x <- as.matrix(x)
    storage.mode(x) <- "double"
    n <- nrow(x)
    d <- ncol(x)

    base_ld <- log_density(x)
    if (!is.numeric(base_ld) || length(base_ld) != n) {
      stop("`log_density(x)` must return a numeric vector of length nrow(x) (",
           n, "); got length ", length(base_ld), ".", call. = FALSE)
    }

    s <- matrix(0, nrow = n, ncol = d)
    for (j in seq_len(d)) {
      x_plus <- x
      x_minus <- x
      x_plus[, j] <- x_plus[, j] + h
      x_minus[, j] <- x_minus[, j] - h
      s[, j] <- (log_density(x_plus) - log_density(x_minus)) / (2 * h)
    }

    if (any(!is.finite(s))) {
      stop("Finite-difference score contains non-finite values; check ",
           "`log_density` or adjust `h`.", call. = FALSE)
    }
    s
  }
}
