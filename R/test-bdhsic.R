#' Backdoor-HSIC Test for Causal Association
#'
#' Tests the do-null hypothesis H_0: p(y | do(x)) = p*(y) using a
#' kernel-based test with backdoor adjustment via density ratio
#' estimation. Detects causal associations including non-linear effects
#' that standard linear methods miss.
#'
#' @param x Numeric vector or matrix. Treatment variable.
#' @param y Numeric vector or matrix. Outcome variable.
#' @param z Numeric matrix, data.frame, or data.table. Confounders.
#' @param kernel_x Kernel specification for treatment space. Default is
#'   RBF with median heuristic.
#' @param kernel_y Kernel specification for outcome space. Default is
#'   RBF with median heuristic.
#' @param density_ratio Character. Method for density ratio estimation:
#'   `"logistic"` (default), `"ranger"`, `"xgboost"`, `"proxymix"`, or
#'   `"rulsif"`. The `"proxymix"` backend fits Gaussian-mixture proxies
#'   to the joint and product-of-marginals sample clouds via
#'   classical EM (Hoek & Elliott, 2024), giving a parametric
#'   alternative to NCE-based classifiers; useful for multimodal
#'   densities or when classifier calibration is unreliable. Requires
#'   the `proxymix` package (`>= 0.3.0`).
#' @param n_permutations Integer. Number of permutations for the null
#'   distribution. Default is 500.
#' @param n_clusters Integer or `"auto"`. Number of *propensity* clusters
#'   for valid permutation when `cluster_id = NULL`. Default is `"auto"`.
#' @param split_ratio Numeric in (0, 1). Proportion of data for training
#'   the density ratio estimator. Default is 0.5.
#' @param alpha Numeric. Significance level. Default is 0.05.
#' @param seed Integer or `NULL`. Random seed for reproducibility.
#' @param verbose Logical. Print progress. Default is `FALSE`.
#' @param cluster_id Optional vector of length `nrow(x)` identifying
#'   external clusters (e.g. site, season, paddock, farm). When supplied,
#'   the permutation null is built by within-cluster reshuffling of `y`,
#'   which preserves cluster-level effects in the null. Coerced to
#'   factor; the test split inherits the cluster assignment. The result
#'   then carries a per-cluster stratified bd-HSIC alongside the pooled
#'   statistic.
#' @param permutation Character. Permutation scheme:
#'   * `"auto"` (default) -- when `cluster_id` is supplied, equivalent to
#'     `"within_cluster"`; otherwise falls back to k-means clustering on
#'     propensity weights (the original Hu/Sejdinovic/Evans scheme).
#'   * `"within_cluster"` -- requires `cluster_id`. Permutes `y` indices
#'     only within clusters; preserves cluster-level effects.
#'   * `"naive"` -- unrestricted permutation across all observations.
#'     Use only when independence within clusters is plausible (rarely
#'     true in ag-systems data).
#' @param min_ess_fraction Numeric in `(0, 1)` or `0` / non-finite to
#'   disable. ESS-floor reliability gate (added 0.0.0.9013): if the
#'   weighted-HSIC effective sample size is below
#'   `min_ess_fraction * n_test`, a warning is emitted and
#'   `result$ess_warning` is `TRUE`. The default `0.1` (10%) is a
#'   conservative floor; tighten it for studies with strict reliability
#'   requirements.
#'
#' @section Train/test split (0.0.0.9014):
#' The density-ratio estimator is now **fit on the train split and
#' predicted on the held-out test split** via [fit_density_ratio()] +
#' [predict_density_ratio()]. The documented `split_ratio` is honoured
#' end-to-end for all four classifier / proxymix backends. The 0.0.0.9013
#' sample-split leak warning is therefore retired. RuLSIF, the
#' kernel-based closed-form backend, still uses [estimate_rulsif()] on
#' the train/test split natively.
#'
#' The fitted density-ratio model is preserved on
#' `result$density_ratio_fit` for callers that want backend
#' diagnostics (see `?fit_density_ratio` Value; proxymix exposes BIC,
#' AIC, log-likelihood, convergence per GMM).
#'
#' @return An object of class `"kernel_test_result"`. When `cluster_id`
#'   is supplied, the result additionally carries:
#'   \describe{
#'     \item{permutation_scheme}{Character: which scheme was used.}
#'     \item{cluster_id}{Integer cluster assignment on the test split.}
#'     \item{cluster_levels}{Character cluster labels.}
#'     \item{per_cluster_statistic}{Per-cluster weighted HSIC (stratified
#'       contributions); `NA` for clusters with `< 2` test observations.}
#'   }
#'
#' @details
#' The bd-HSIC test (Hu, Sejdinovic & Evans, 2024) tests whether
#' treatment X has a causal effect on outcome Y after adjusting for
#' confounders Z via the backdoor criterion.
#'
#' The test works by:
#' 1. Estimating density ratios w(x, z) = p*(x) / p(x|z) to reweight
#'    observational samples to the interventional distribution.
#' 2. Computing a weighted HSIC statistic between X and Y.
#' 3. Obtaining p-values via permutation of Y within exchangeability
#'    clusters -- propensity-similarity clusters by default, or external
#'    design clusters (site / season / paddock) when `cluster_id` is
#'    supplied.
#'
#' **Hierarchical extension.** When the design is naturally clustered
#' (multi-site agricultural trials, paddock x season factorial designs,
#' patient x hospital data), supplying `cluster_id` activates
#' within-cluster permutation: indices of `y` are reshuffled only within
#' each cluster, preserving cluster-level effects in the null. This is
#' the safer default for clustered data; naive permutation across
#' clusters can inflate Type I error when cluster effects exist.
#'
#' Unlike PDS or Double ML, bd-HSIC can detect non-linear causal effects
#' (e.g., U-shaped relationships) where the treatment affects higher
#' moments of the outcome but not necessarily the mean.
#'
#' @references
#' Hu, R., Sejdinovic, D., & Evans, R. J. (2024). A kernel test for
#' causal association via noise contrastive backdoor adjustment. *JMLR*,
#' 25(160), 1-56.
#'
#' @examples
#' set.seed(42)
#' n <- 300
#' z <- matrix(rnorm(n * 2), n, 2)
#' x <- z[, 1] + rnorm(n)
#' y <- 0.5 * x + z[, 2] + rnorm(n, sd = 0.5)
#'
#' result <- bd_hsic_test(x, y, z, n_permutations = 200, seed = 1)
#' print(result)
#'
#' @export
bd_hsic_test <- function(x, y, z,
                         kernel_x = kernel_spec(),
                         kernel_y = kernel_spec(),
                         density_ratio = c("logistic", "ranger", "xgboost",
                                           "proxymix", "rulsif"),
                         n_permutations = 500L,
                         n_clusters = "auto",
                         split_ratio = 0.5,
                         alpha = 0.05,
                         seed = NULL,
                         verbose = FALSE,
                         cluster_id = NULL,
                         permutation = c("auto", "within_cluster", "naive"),
                         min_ess_fraction = 0.1) {
  cl <- match.call()
  density_ratio <- match.arg(density_ratio)
  permutation <- match.arg(permutation)
  n_permutations <- as.integer(n_permutations)

  x <- validate_input(x, "x", min_n = 20)
  y <- validate_input(y, "y", min_n = 20)
  z <- validate_input(z, "z", min_n = 20)
  n <- nrow(x)

  if (nrow(y) != n || nrow(z) != n) {
    stop("`x`, `y`, and `z` must have the same number of observations.",
      call. = FALSE
    )
  }
  if (!is.null(cluster_id)) {
    if (length(cluster_id) != n) {
      stop("`cluster_id` length must match `nrow(x)`.", call. = FALSE)
    }
    cluster_id <- as.factor(cluster_id)
    if (nlevels(cluster_id) < 2L) {
      stop("`cluster_id` must define at least 2 distinct clusters.",
           call. = FALSE)
    }
  }
  if (permutation == "within_cluster" && is.null(cluster_id)) {
    stop("permutation = \"within_cluster\" requires `cluster_id`.",
         call. = FALSE)
  }

  if (!is.null(seed)) set.seed(seed)

  # Split data: train density ratio on one half, test on the other
  n_train <- floor(n * split_ratio)
  idx_train <- sample.int(n, n_train)
  idx_test <- setdiff(seq_len(n), idx_train)
  n_test <- length(idx_test)

  if (verbose) message("Training density ratio estimator on ", n_train, " obs...")

  # Honest train/test split (kernR 0.0.0.9014): fit the density-ratio
  # estimator on the training split, predict on the held-out test
  # split. Closes the P0 #2 sample-split leak surfaced by the
  # 2026-05-16 critical review.
  if (density_ratio == "rulsif") {
    dr <- estimate_rulsif(x[idx_train, , drop = FALSE],
      x[idx_test, , drop = FALSE],
      kernel = kernel_x
    )
    weights <- dr$weights
    dr_fit  <- NULL
  } else {
    dr_fit <- fit_density_ratio(
      x = x[idx_train, , drop = FALSE],
      z = z[idx_train, , drop = FALSE],
      method = density_ratio
    )
    weights <- predict_density_ratio(
      dr_fit,
      new_x = x[idx_test, , drop = FALSE],
      new_z = z[idx_test, , drop = FALSE],
      type  = "weight"
    )
  }

  ess <- effective_sample_size(weights)

  if (verbose) message("ESS: ", round(ess, 1), " / ", n_test)

  # ESS-floor reliability gate (0.0.0.9013). When ESS collapses
  # the weighted statistic is effectively driven by a handful of
  # points; surface this loudly rather than reporting a finite-but-
  # noisy p-value as if it were trustworthy.
  if (is.finite(min_ess_fraction) && min_ess_fraction > 0 &&
      ess < min_ess_fraction * n_test) {
    warning(
      "bd_hsic_test(): ESS (", round(ess, 1L), ") is below ",
      formatC(100 * min_ess_fraction, digits = 0, format = "f"),
      "% of n_test (", n_test, "). The weighted test statistic ",
      "is dominated by a small number of high-weight observations; ",
      "the resulting p-value is not a reliable verdict. Increase ",
      "n, switch density_ratio backend, or tighten the design.",
      call. = FALSE
    )
    ess_warning <- TRUE
  } else {
    ess_warning <- FALSE
  }

  # Compute kernel matrices on test data
  kernel_x <- resolve_bandwidth(kernel_x, x[idx_test, , drop = FALSE])
  kernel_y <- resolve_bandwidth(kernel_y, y[idx_test, , drop = FALSE])

  Kx <- kernel_matrix(x[idx_test, , drop = FALSE], kernel = kernel_x)
  Ky <- kernel_matrix(y[idx_test, , drop = FALSE], kernel = kernel_y)

  # Observed weighted HSIC
  stat_obs <- weighted_hsic_stat_cpp(Kx, Ky, weights)

  if (verbose) message("Computing ", n_permutations, " permutations...")

  # Choose permutation scheme
  if (!is.null(cluster_id)) {
    cluster_test_int <- as.integer(cluster_id[idx_test])
    cluster_test_int <- match(cluster_test_int, sort(unique(cluster_test_int)))
    if (length(unique(cluster_test_int)) < 2L) {
      stop("`cluster_id` must define at least 2 distinct clusters in the test split.",
           call. = FALSE)
    }
    if (permutation == "naive") {
      clusters <- rep(1L, n_test)
      permutation_scheme <- "naive"
    } else {
      clusters <- cluster_test_int
      permutation_scheme <- "within_cluster"
    }
    K <- max(cluster_test_int)
    per_cluster <- vapply(seq_len(K), function(k) {
      idx <- which(cluster_test_int == k)
      if (length(idx) < 2L) return(NA_real_)
      weighted_hsic_stat_cpp(
        Kx[idx, idx, drop = FALSE],
        Ky[idx, idx, drop = FALSE],
        weights[idx]
      )
    }, numeric(1L))
    cluster_levels_test <- levels(cluster_id)[
      sort(unique(as.integer(cluster_id[idx_test])))
    ]
    names(per_cluster) <- cluster_levels_test
  } else {
    if (permutation == "naive") {
      clusters <- rep(1L, n_test)
      permutation_scheme <- "naive"
    } else {
      clusters <- cluster_observations(
        weights,
        z[idx_test, , drop = FALSE],
        n_clusters = n_clusters
      )
      permutation_scheme <- "propensity"
    }
    cluster_test_int <- NULL
    cluster_levels_test <- NULL
    per_cluster <- NULL
  }

  # Permutation null distribution
  null_dist <- cluster_permutation_hsic(
    Kx, Ky, weights, clusters, n_permutations
  )

  # p-value
  p_value <- (1 + sum(null_dist >= stat_obs)) / (1 + n_permutations)

  structure(
    list(
      statistic = stat_obs,
      p_value = p_value,
      method = "bd-HSIC",
      n = n_test,
      n_permutations = n_permutations,
      null_distribution = null_dist,
      ess = ess,
      ess_fraction = ess / n_test,
      ess_warning = ess_warning,
      min_ess_fraction = min_ess_fraction,
      density_ratio_fit = dr_fit,
      weights = weights,
      kernel_x = kernel_x,
      kernel_y = kernel_y,
      permutation_scheme = permutation_scheme,
      cluster_id = cluster_test_int,
      cluster_levels = cluster_levels_test,
      per_cluster_statistic = per_cluster,
      call = cl
    ),
    class = "kernel_test_result"
  )
}
