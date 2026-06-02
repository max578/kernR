# test-concordance.R -- k-sample kernel concordance test over ensembles.
#
# Implements concordance_test(), a kernel k-sample test asking whether several
# samples -- typically posterior draws from different inference engines or
# scenario ensembles from different simulators -- are mutually consistent at
# the distributional level. The statistic is the summed pairwise unbiased MMD,
# calibrated by a joint label permutation of the pooled sample; the returned
# pairwise discrepancy matrix shows not just whether the sources disagree but
# which pair disagrees and by how much. This is the constellation validation
# role made concrete: convergence across engines is evidence, divergence
# localises the assumption-load.

#' Kernel k-sample Concordance Test
#'
#' Tests whether two or more samples come from a common distribution, using the
#' summed pairwise Maximum Mean Discrepancy with a joint-permutation null. The
#' samples are typically posterior draws from different inference engines, or
#' scenario ensembles from different simulators; the test asks whether they are
#' mutually concordant. Unlike repeated two-sample [mmd_test()] calls, the null
#' is a single shared relabeling of the pooled sample, so the family-wise error
#' is controlled and the overall verdict is one calibrated p-value.
#'
#' The returned object carries the full pairwise MMD discrepancy matrix, so a
#' rejection can be read down to the offending pair: convergence across sources
#' is corroborating evidence, and divergence localises which source departs and
#' on which margin. This is the cross-engine concordance role -- a sample-based
#' complement to the score-based [ksd_test()].
#'
#' @param x A list of two or more samples to compare. Each element is a numeric
#'   vector, matrix, or data.frame with `n_k` observations (rows) over a shared
#'   `d` columns. A named list labels the sources in the output; an unnamed list
#'   is labelled `Source 1`, `Source 2`, and so on. Each sample needs at least
#'   five rows.
#' @param kernel Kernel specification. Default is RBF with the median heuristic
#'   over the pooled sample.
#' @param n_permutations Integer. Number of joint permutations for the null.
#'   Default `500`.
#' @param alpha Numeric in `(0, 1)`. Significance level for the verdict.
#'   Default `0.05`.
#' @param seed Integer or `NULL`. Random seed for reproducibility.
#'
#' @returns An object of class `c("concordance_test", "kernel_test_result")`
#'   carrying the standard `kernel_test_result` fields plus:
#'   \describe{
#'     \item{statistic}{Summed pairwise unbiased MMD-squared.}
#'     \item{p_value}{Upper-tail joint-permutation p-value (with `+1`
#'       correction).}
#'     \item{n_groups}{Number of sources compared.}
#'     \item{group_sizes}{Named integer vector of per-source sample sizes.}
#'     \item{pairwise}{Symmetric `K x K` matrix of pairwise unbiased
#'       MMD-squared, row/column-named by source.}
#'     \item{alpha, reject}{Verdict level and `p_value <= alpha`.}
#'   }
#'
#' @references
#' Gretton, A., Borgwardt, K. M., Rasch, M. J., Scholkopf, B., & Smola, A.
#' (2012). A kernel two-sample test. *Journal of Machine Learning Research*,
#' 13, 723-773.
#'
#' @examples
#' set.seed(1)
#'
#' # Three concordant sources (same distribution): not rejected
#' draws <- list(
#'   engine_a = matrix(stats::rnorm(400L), ncol = 2L),
#'   engine_b = matrix(stats::rnorm(400L), ncol = 2L),
#'   engine_c = matrix(stats::rnorm(400L), ncol = 2L)
#' )
#' fit_ok <- concordance_test(draws, n_permutations = 199L, seed = 1L)
#' fit_ok
#'
#' # One source departs (mean-shifted): rejected, and the pairwise matrix
#' # localises engine_c
#' draws$engine_c <- draws$engine_c + 1
#' fit_bad <- concordance_test(draws, n_permutations = 199L, seed = 1L)
#' fit_bad$pairwise
#'
#' @seealso [mmd_test()], [ksd_test()], [mmd_ppc()]
#' @family goodness-of-fit tests
#' @author Max Moldovan, \email{max.moldovan@@adelaide.edu.au}
#' @export
concordance_test <- function(x,
                             kernel = kernel_spec(),
                             n_permutations = 500L,
                             alpha = 0.05,
                             seed = NULL) {
  cl <- match.call()
  groups <- .check_concordance_input(x, alpha)
  n_permutations <- as.integer(n_permutations)

  group_names <- names(groups)
  sizes <- vapply(groups, nrow, integer(1L))
  k <- length(groups)

  if (!is.null(seed)) set.seed(seed)

  # Pooled kernel matrix and per-source index blocks ---------------------
  pooled <- do.call(rbind, groups)
  kernel <- resolve_bandwidth(kernel, pooled)
  k_pool <- kernel_matrix(pooled, kernel = kernel)

  ends <- cumsum(sizes)
  starts <- c(1L, ends[-length(ends)] + 1L)

  # Observed pairwise unbiased MMD^2 + summed statistic ------------------
  pairwise <- matrix(0, nrow = k, ncol = k,
                     dimnames = list(group_names, group_names))
  for (a in seq_len(k - 1L)) {
    idx_a <- starts[a]:ends[a]
    for (b in (a + 1L):k) {
      idx_b <- starts[b]:ends[b]
      kxx <- k_pool[idx_a, idx_a, drop = FALSE]
      kyy <- k_pool[idx_b, idx_b, drop = FALSE]
      kxy <- k_pool[idx_a, idx_b, drop = FALSE]
      mmd2 <- mmd2_unbiased_cpp(kxx, kyy, kxy)
      pairwise[a, b] <- mmd2
      pairwise[b, a] <- mmd2
    }
  }
  stat_obs <- sum(pairwise[upper.tri(pairwise)])

  # Joint-permutation null and the upper-tail p-value -------------------
  null_dist <- permutation_ksample_mmd_cpp(k_pool, as.integer(sizes),
                                           n_permutations)
  p_value <- (1 + sum(null_dist >= stat_obs)) / (1 + n_permutations)

  structure(
    list(
      statistic = stat_obs,
      p_value = p_value,
      method = "Concordance",
      n = sum(sizes),
      n_permutations = n_permutations,
      null_distribution = as.numeric(null_dist),
      ess = NA_real_,
      weights = NULL,
      kernel_x = kernel,
      kernel_y = NULL,
      n_groups = k,
      group_sizes = sizes,
      pairwise = pairwise,
      surprise_bits = -log2(p_value),
      alpha = alpha,
      reject = p_value <= alpha,
      call = cl
    ),
    class = c("concordance_test", "kernel_test_result")
  )
}

#' @export
print.concordance_test <- function(x, ...) {
  NextMethod()
  cat("Concordance verdict\n")
  cat("  Sources:    ", x$n_groups, " (",
      paste(names(x$group_sizes), collapse = ", "), ")\n", sep = "")
  cat("  Verdict:    ",
      if (isTRUE(x$reject)) {
        "REJECT (sources are not mutually concordant)"
      } else {
        "concordant"
      },
      "\n", sep = "")
  cat("  Pairwise MMD-squared:\n")
  pw <- formatC(x$pairwise, digits = 3L, format = "g")
  print(noquote(pw))
  cat("\n")
  invisible(x)
}

# Internal validation --------------------------------------------------------

#' Validate and normalise the concordance_test() input list
#'
#' @param x The user-supplied list of samples.
#' @param alpha Verdict level.
#' @returns A named list of coerced numeric matrices with a shared column count.
#' @noRd
#' @keywords internal
.check_concordance_input <- function(x, alpha) {
  if (!is.list(x) || length(x) < 2L) {
    stop("`x` must be a list of at least two samples to compare.",
         call. = FALSE)
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a single number in (0, 1).", call. = FALSE)
  }

  groups <- lapply(x, function(g) {
    g <- as.matrix(g)
    storage.mode(g) <- "double"
    g
  })

  ncols <- vapply(groups, ncol, integer(1L))
  if (length(unique(ncols)) != 1L) {
    stop("All samples in `x` must have the same number of columns; got ",
         paste(ncols, collapse = ", "), ".", call. = FALSE)
  }
  nrows <- vapply(groups, nrow, integer(1L))
  if (any(nrows < 5L)) {
    stop("Each sample in `x` must have at least 5 rows (observations).",
         call. = FALSE)
  }
  if (any(vapply(groups, function(g) any(!is.finite(g)), logical(1L)))) {
    stop("All samples in `x` must contain only finite values.",
         call. = FALSE)
  }

  if (is.null(names(groups)) || any(!nzchar(names(groups)))) {
    names(groups) <- paste("Source", seq_along(groups))
  }
  groups
}
