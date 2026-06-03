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
#' The exact test materialises the pooled `n x n` kernel matrix (`O(n^2)`). To
#' keep large ensembles tractable without a silent loss of exactness, a pooled
#' sample with more than `n_exact_max` rows is delegated to
#' [concordance_test_nystrom()] -- a low-rank approximation that is *announced*
#' by a message and *recorded* in the returned object's `approximation` and `m`
#' fields. Set `n_exact_max = Inf` to force the exact test, or call
#' [concordance_test_nystrom()] directly to control the approximation rank `m`.
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
#' @param n_exact_max Integer or `Inf`. Pooled-sample-size ceiling for the
#'   exact `O(n^2)` test. Above it, the call is delegated to
#'   [concordance_test_nystrom()] (with a message; the verdict object records
#'   `approximation = "nystrom"`). `Inf` forces the exact test at any size.
#'   Default `5000L`.
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
                             seed = NULL,
                             n_exact_max = 5000L) {
  cl <- match.call()
  groups <- .check_concordance_input(x, alpha)
  n_permutations <- as.integer(n_permutations)

  group_names <- names(groups)
  sizes <- vapply(groups, nrow, integer(1L))
  k <- length(groups)

  # Auto-switch to the low-rank test above n_exact_max: announced, recorded
  # in the object, escapable via Inf. See the Details section.
  if (!is.numeric(n_exact_max) || length(n_exact_max) != 1L ||
      is.na(n_exact_max) || n_exact_max < 1) {
    stop("`n_exact_max` must be a single positive number or `Inf`.",
         call. = FALSE)
  }
  n_total <- sum(sizes)
  if (n_total > n_exact_max) {
    m_auto <- 100L
    message("concordance_test(): pooled n = ", n_total,
            " > n_exact_max = ", n_exact_max,
            "; delegating to concordance_test_nystrom(m = ", m_auto,
            "). Call concordance_test_nystrom() to control the rank, or set ",
            "n_exact_max = Inf to force the exact O(n^2) test.")
    return(concordance_test_nystrom(x, kernel = kernel, m = m_auto,
                                    n_permutations = n_permutations,
                                    alpha = alpha, seed = seed))
  }

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

# Accelerated (low-rank) concordance ----------------------------------------

#' Internal: per-source column sums and squared-norm (diagonal) sums
#'
#' Given an `n x m` kernel factor `F` (with `F F^\top \approx K`), the row
#' grouping `g` (integer labels `1..k`), and the diagonal `rn = rowSums(F^2)`,
#' returns the `k x m` matrix of per-source column sums `S` and the length-`k`
#' vector of per-source diagonal sums `D`. Both are computed in `O(n m)` via
#' [base::rowsum()], which orders its output by sorted group label.
#'
#' @param Fmat The `n x m` factor matrix.
#' @param rn Length-`n` vector of row squared norms, `rowSums(Fmat^2)`.
#' @param g Integer length-`n` group labels in `1..k`.
#' @param k Number of sources.
#' @returns A list with `S` (`k x m`) and `D` (length `k`).
#' @noRd
#' @keywords internal
.concordance_block_moments <- function(Fmat, rn, g, k) {
  gf <- factor(g, levels = seq_len(k))
  S <- rowsum(Fmat, gf, reorder = FALSE)
  D <- as.numeric(rowsum(rn, gf, reorder = FALSE))
  list(S = S, D = D)
}

#' Internal: pairwise unbiased MMD^2 from block moments
#'
#' Reconstructs the `k x k` pairwise unbiased MMD-squared matrix from the
#' per-source column sums `S` and diagonal sums `D` using the low-rank
#' identity verified against [mmd2_unbiased_cpp()]: for sources `a`, `b`,
#' \eqn{\mathrm{MMD}^2 = (\lVert S_a \rVert^2 - D_a) / (n_a (n_a - 1)) +
#' (\lVert S_b \rVert^2 - D_b) / (n_b (n_b - 1)) -
#' 2 \langle S_a, S_b \rangle / (n_a n_b)}.
#'
#' @param S `k x m` per-source column sums.
#' @param D Length-`k` per-source diagonal sums.
#' @param sizes Length-`k` integer source sizes.
#' @returns Symmetric `k x k` matrix of pairwise unbiased MMD-squared.
#' @noRd
#' @keywords internal
.concordance_pairwise_lowrank <- function(S, D, sizes) {
  k <- length(sizes)
  within <- (rowSums(S * S) - D) / (sizes * (sizes - 1))
  pw <- matrix(0, nrow = k, ncol = k)
  for (a in seq_len(k - 1L)) {
    for (b in (a + 1L):k) {
      cross <- 2 * sum(S[a, ] * S[b, ]) / (sizes[a] * sizes[b])
      val <- within[a] + within[b] - cross
      pw[a, b] <- val
      pw[b, a] <- val
    }
  }
  pw
}

#' Accelerated Kernel k-sample Concordance Test (Nystrom / RFF)
#'
#' Low-rank counterpart to [concordance_test()] for large ensembles. The
#' pooled sample is factorised once -- by the Nystrom method (default) or
#' random Fourier features -- into an `n x m` factor `F` with
#' \eqn{F F^\top \approx K}; the summed pairwise unbiased MMD-squared and its
#' joint-permutation null are then computed from `F` in `O(n m)` per
#' permutation rather than `O(n^2)`, with `m << n` controlling the
#' speed / accuracy trade-off. The verdict object and its interpretation are
#' identical to [concordance_test()]; only the cost scales differently.
#'
#' The factorisation of the pooled sample preserves the per-source mean
#' embeddings (per-source column sums of `F`), so the pairwise discrepancy
#' matrix still localises which source departs. The joint-permutation null is
#' built by relabelling the rows of `F` -- the low-rank analogue of permuting
#' the pooled-sample labels in [concordance_test()].
#'
#' Use [concordance_test()] for exact results at moderate `n`; reach for this
#' function when the pooled sample is large enough that the `O(n^2)` kernel
#' matrix is the bottleneck. RFF (`method = "rff"`) requires an RBF kernel;
#' Nystrom supports any [kernel_spec()].
#'
#' @inheritParams concordance_test
#' @param method Character. `"nystrom"` (default) or `"rff"`. RFF requires an
#'   RBF `kernel`.
#' @param m Integer. Rank of the approximation: the number of Nystrom
#'   landmarks or RFF features. Larger `m` improves accuracy at higher cost.
#'   Default `100L`.
#' @param regularise Small positive numeric. Ridge added before the Nystrom
#'   Cholesky for numerical stability; ignored under `method = "rff"`.
#'   Default `1e-6`.
#'
#' @returns An object of class `c("concordance_test", "kernel_test_result")`
#'   carrying the same fields as [concordance_test()] plus:
#'   \describe{
#'     \item{approximation}{`"nystrom"` or `"rff"`.}
#'     \item{m}{Effective rank used for the factorisation.}
#'   }
#'
#' @references
#' Gretton, A., Borgwardt, K. M., Rasch, M. J., Scholkopf, B., & Smola, A.
#' (2012). A kernel two-sample test. *Journal of Machine Learning Research*,
#' 13, 723-773.
#'
#' Williams, C. K. I., & Seeger, M. (2001). Using the Nystrom method to speed
#' up kernel machines. *NeurIPS*, 13.
#'
#' Rahimi, A., & Recht, B. (2007). Random features for large-scale kernel
#' machines. *NeurIPS*, 20.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' big <- list(
#'   engine_a = matrix(stats::rnorm(4000L), ncol = 2L),
#'   engine_b = matrix(stats::rnorm(4000L), ncol = 2L),
#'   engine_c = matrix(stats::rnorm(4000L), ncol = 2L) + 0.4
#' )
#' fit <- concordance_test_nystrom(big, m = 80L,
#'                                 n_permutations = 199L, seed = 1L)
#' fit
#' fit$pairwise
#' }
#'
#' @seealso [concordance_test()], [nystrom_factor()], [rff_features()],
#'   [hsic_test_nystrom()]
#' @family goodness-of-fit tests
#' @family low-rank acceleration
#' @author Max Moldovan, \email{max.moldovan@@adelaide.edu.au}
#' @export
concordance_test_nystrom <- function(x,
                                     kernel = kernel_spec(),
                                     method = c("nystrom", "rff"),
                                     m = 100L,
                                     n_permutations = 500L,
                                     alpha = 0.05,
                                     seed = NULL,
                                     regularise = 1e-6) {
  cl <- match.call()
  method <- match.arg(method)
  groups <- .check_concordance_input(x, alpha)
  n_permutations <- as.integer(n_permutations)

  group_names <- names(groups)
  sizes <- vapply(groups, nrow, integer(1L))
  k <- length(groups)
  n_total <- sum(sizes)

  if (!is.null(seed)) set.seed(seed)

  # Factorise the pooled sample once (bandwidth resolved on the pool) ----
  pooled <- do.call(rbind, groups)
  fac <- if (method == "nystrom") {
    nystrom_factor(pooled, kernel = kernel, m = m, regularise = regularise)
  } else {
    rff_features(pooled, kernel = kernel, D = m)
  }
  Fmat <- fac$F
  rn <- rowSums(Fmat * Fmat)

  # Fixed row -> source labels for the observed statistic ----------------
  g_obs <- rep.int(seq_len(k), times = sizes)

  observed <- .concordance_block_moments(Fmat, rn, g_obs, k)
  pairwise <- .concordance_pairwise_lowrank(observed$S, observed$D, sizes)
  dimnames(pairwise) <- list(group_names, group_names)
  stat_obs <- sum(pairwise[upper.tri(pairwise)])

  # Joint-permutation null: relabel the rows of the factor ---------------
  null_dist <- numeric(n_permutations)
  for (p in seq_len(n_permutations)) {
    g_perm <- g_obs[sample.int(n_total)]
    mom <- .concordance_block_moments(Fmat, rn, g_perm, k)
    pw <- .concordance_pairwise_lowrank(mom$S, mom$D, sizes)
    null_dist[p] <- sum(pw[upper.tri(pw)])
  }
  p_value <- (1 + sum(null_dist >= stat_obs)) / (1 + n_permutations)

  structure(
    list(
      statistic = stat_obs,
      p_value = p_value,
      method = paste0("Concordance (", method, ")"),
      n = n_total,
      n_permutations = n_permutations,
      null_distribution = null_dist,
      ess = NA_real_,
      weights = NULL,
      kernel_x = fac$kernel,
      kernel_y = NULL,
      n_groups = k,
      group_sizes = sizes,
      pairwise = pairwise,
      approximation = method,
      m = fac$m,
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
