#' Doubly Robust Distributional Average Treatment Effect Test (DR-DATE)
#'
#' Tests whether the distributions of potential outcomes Y(1) and Y(0)
#' differ using a doubly robust kernel MMD statistic. Detects
#' distributional effects (variance, shape) that mean-based tests miss.
#'
#' @param y Numeric vector or matrix. Outcome variable.
#' @param treatment Binary vector (0/1). Treatment indicator.
#' @param covariates Numeric matrix, data.frame, or data.table. Confounders.
#' @param kernel_y Kernel specification for outcome space. Default is RBF.
#' @param propensity_model Character. Method for propensity estimation:
#'   `"logistic"` (default), `"ranger"`, or `"xgboost"`.
#' @param outcome_model Character. `"krr"` (kernel ridge regression,
#'   default) fits a conditional mean embedding for each arm and forms
#'   the doubly robust statistic; `"zero"` drops the outcome model and
#'   returns the inverse-probability-weighted (singly robust) statistic.
#' @param cross_fit Logical. If `TRUE` (default), both nuisances -- the
#'   propensity score and the conditional mean embedding -- are estimated
#'   by `n_folds`-fold cross-fitting and evaluated out-of-fold, as the
#'   doubly robust theory requires under flexible nuisance estimators. If
#'   `FALSE`, both are fit in-sample (faster, but the test can be
#'   anti-conservative).
#' @param n_folds Integer. Number of cross-fitting folds. Default is 5.
#' @param n_permutations Integer. Number of permutations. Default is 500.
#' @param n_bins Integer. Propensity score bins for permutation.
#'   Default is 10.
#' @param regularisation Numeric or `"cv"`. Ridge parameter for the CME.
#'   Default is `"cv"`.
#' @param min_ess_fraction Numeric in (0, 1). If the effective sample
#'   size of either arm's inverse-probability weights falls below this
#'   fraction of `n`, a reliability `warning()` is emitted. Default 0.1.
#' @param alpha Numeric. Significance level. Default is 0.05.
#' @param seed Integer or `NULL`. Random seed. Permutations are drawn
#'   through R's RNG, so a fixed `seed` makes the test fully reproducible.
#' @param verbose Logical. Print progress. Default is `FALSE`.
#'
#' @return An object of class `"kernel_test_result"`. The `ess` element
#'   holds the smaller of the two per-arm effective sample sizes and
#'   `ess_warning` records whether the reliability floor was hit.
#'
#' @details
#' The DR-DATE test (Fawkes, Hu, Evans & Sejdinovic, 2024) constructs
#' doubly robust (augmented inverse-probability-weighted) estimators for
#' the counterfactual mean embeddings of Y(1) and Y(0) in a reproducing
#' kernel Hilbert space. For arm \eqn{a} the augmented embedding is
#' \deqn{\hat\mu_a = \frac{1}{n}\sum_i \tilde w_{a,i}
#'   \bigl(k(y_i,\cdot) - \hat m_a(x_i)\bigr) + \hat m_a(x_i),}
#' where \eqn{\hat m_a} is the conditional mean embedding fitted on arm
#' \eqn{a} and \eqn{\tilde w_{a,i}} are stabilised inverse-probability
#' weights. The statistic is \eqn{\|\hat\mu_1 - \hat\mu_0\|^2} in the
#' RKHS. Setting `outcome_model = "zero"` sets \eqn{\hat m_a \equiv 0}
#' and recovers the inverse-probability-weighted statistic.
#'
#' **Double robustness**: the test is consistent if *either* the
#' propensity model *or* the outcome (CME) model is correctly specified.
#' Cross-fitting (`cross_fit = TRUE`) makes this hold under flexible
#' machine-learning nuisances by removing own-observation overfitting
#' bias (Chernozhukov et al., 2018).
#'
#' **Permutation null**: the reference distribution permutes treatment
#' labels *within* propensity-score bins, holding the fitted nuisances
#' fixed. This is valid under within-bin exchangeability of treatment
#' given the (binned) propensity score; calibration degrades as the bins
#' coarsen relative to the propensity variation inside them.
#'
#' **Key advantage**: unlike DML or TMLE which test only for mean
#' shifts, DR-DATE detects *any* distributional difference including
#' changes in variance, skewness, or shape.
#'
#' @references
#' Fawkes, J., Hu, R., Evans, R. J., & Sejdinovic, D. (2024). Doubly
#' robust kernel statistics for testing distributional treatment effects.
#' *Transactions on Machine Learning Research*.
#'
#' Chernozhukov, V., Chetverikov, D., Demirer, M., Duflo, E., Hansen, C.,
#' Newey, W., & Robins, J. (2018). Double/debiased machine learning for
#' treatment and structural parameters. *The Econometrics Journal*,
#' 21(1), C1-C68.
#'
#' @examples
#' set.seed(42)
#' n <- 300
#' x <- matrix(rnorm(n * 2), n, 2)
#' logit_p <- 0.5 * x[, 1]
#' t <- rbinom(n, 1, plogis(logit_p))
#' y <- t * 1.0 + x[, 1] + rnorm(n, sd = 0.5)
#'
#' result <- dr_date_test(y, t, x, n_permutations = 200, seed = 1)
#' print(result)
#'
#' @family distributional treatment effects
#' @export
dr_date_test <- function(y, treatment, covariates,
                         kernel_y = kernel_spec(),
                         propensity_model = c("logistic", "ranger", "xgboost"),
                         outcome_model = c("krr", "zero"),
                         cross_fit = TRUE,
                         n_folds = 5L,
                         n_permutations = 500L,
                         n_bins = 10L,
                         regularisation = "cv",
                         min_ess_fraction = 0.1,
                         alpha = 0.05,
                         seed = NULL,
                         verbose = FALSE) {
  cl <- match.call()
  propensity_model <- match.arg(propensity_model)
  outcome_model <- match.arg(outcome_model)
  n_permutations <- as.integer(n_permutations)
  n_folds <- as.integer(n_folds)

  y <- validate_input(y, "y", min_n = 30)
  covariates <- validate_input(covariates, "covariates", min_n = 30)
  treatment <- as.integer(treatment)
  n <- nrow(y)

  if (length(treatment) != n || nrow(covariates) != n) {
    stop("All inputs must have the same number of observations.", call. = FALSE)
  }
  if (!all(treatment %in% c(0L, 1L))) {
    stop("`treatment` must be binary (0/1).", call. = FALSE)
  }

  if (!is.null(seed)) set.seed(seed)

  # Estimate propensity scores (cross-fitted by default)
  if (verbose) message("Estimating propensity scores...")
  ps <- estimate_propensity(treatment, covariates,
    method = propensity_model,
    cross_fit = cross_fit, n_folds = n_folds
  )
  e_hat <- ps$scores

  # Resolve kernel and compute outcome kernel matrix
  kernel_y <- resolve_bandwidth(kernel_y, y)
  Ky <- kernel_matrix(y, kernel = kernel_y)

  # Conditional mean embedding coefficient matrices (out-of-fold when
  # cross_fit = TRUE). Zero matrices recover the IPW-only statistic.
  warn_reliability <- n_permutations > 0L
  if (outcome_model == "krr") {
    if (verbose) message("Fitting conditional mean embeddings...")
    lambda <- if (identical(regularisation, "cv")) {
      "cv"
    } else {
      as.numeric(regularisation)
    }
    fold_id <- if (cross_fit) {
      sample(rep(seq_len(n_folds), length.out = n))
    } else {
      rep(1L, n)
    }
    C1 <- .cme_coef_matrix(y, covariates, which(treatment == 1L), n,
      lambda, cross_fit, n_folds, fold_id, kernel_y, warn_reliability)
    C0 <- .cme_coef_matrix(y, covariates, which(treatment == 0L), n,
      lambda, cross_fit, n_folds, fold_id, kernel_y, warn_reliability)
  } else {
    C1 <- matrix(0, n, n)
    C0 <- matrix(0, n, n)
  }

  # Effective-sample-size reliability gate on the IPW weights
  ess1 <- effective_sample_size((treatment / e_hat)[treatment == 1L])
  ess0 <- effective_sample_size(((1 - treatment) / (1 - e_hat))[treatment == 0L])
  ess <- min(ess1, ess0)
  ess_warning <- ess < min_ess_fraction * n
  if (warn_reliability && ess_warning) {
    warning(sprintf(
      paste0(
        "Effective sample size %.1f (%.1f%% of n) is below the %.0f%% ",
        "floor; the propensity weights are concentrated and the test ",
        "may be unreliable."
      ),
      ess, 100 * ess / n, 100 * min_ess_fraction
    ), call. = FALSE)
  }

  # Observed statistic
  stat_obs <- compute_dr_date_stat(Ky, treatment, e_hat, C1, C0)

  if (verbose) message("Computing ", n_permutations, " permutations...")
  null_dist <- numeric(n_permutations)
  for (p in seq_len(n_permutations)) {
    t_perm <- bin_permute_treatment(treatment, e_hat, n_bins)
    null_dist[p] <- compute_dr_date_stat(Ky, t_perm, e_hat, C1, C0)
  }

  p_value <- (1 + sum(null_dist >= stat_obs)) / (1 + n_permutations)

  structure(
    list(
      statistic = stat_obs,
      p_value = p_value,
      method = "DR-DATE",
      n = n,
      n_permutations = n_permutations,
      null_distribution = null_dist,
      ess = ess,
      ess_warning = ess_warning,
      weights = e_hat,
      kernel_x = NULL,
      kernel_y = kernel_y,
      call = cl
    ),
    class = "kernel_test_result"
  )
}

#' Compute DR-DATE Test Statistic (augmented IPW)
#'
#' @param Ky n x n outcome kernel matrix.
#' @param treatment Binary treatment vector.
#' @param e_hat Propensity scores.
#' @param C1,C0 n x n conditional-mean-embedding coefficient matrices for
#'   the treated and control arms. Row `i` holds the coefficients of
#'   \eqn{\hat m_a(x_i)} over the n outcome embeddings; zero matrices give
#'   the inverse-probability-weighted statistic.
#'
#' @return Scalar test statistic.
#' @keywords internal
compute_dr_date_stat <- function(Ky, treatment, e_hat, C1, C0) {
  n <- length(treatment)
  t <- treatment

  # Stabilised (Hajek) inverse-probability weights
  w1 <- t / e_hat
  w0 <- (1 - t) / (1 - e_hat)
  w1 <- w1 / mean(w1)
  w0 <- w0 / mean(w0)

  # Augmented (doubly robust) coefficient vectors over the n outcome
  # embeddings: alpha_a = (1/n) [ w_a + C_a^T (1 - w_a) ]. With C_a = 0
  # this reduces to the IPW coefficients w_a / n.
  alpha1 <- (w1 + as.numeric(crossprod(C1, 1 - w1))) / n
  alpha0 <- (w0 + as.numeric(crossprod(C0, 1 - w0))) / n

  d <- alpha1 - alpha0
  mmd2 <- as.numeric(crossprod(d, Ky %*% d))
  max(mmd2, 0)
}

#' Out-of-Fold Conditional-Mean-Embedding Coefficient Matrix
#'
#' Builds the n x n matrix `C` whose row `i` holds the coefficients of
#' the fitted conditional mean embedding \eqn{\hat m(x_i)} over the n
#' outcome embeddings \eqn{k(y_l, \cdot)}. The CME smoother weights are a
#' function of the conditioning (covariate) kernel only, so they are
#' valid coordinates over the global outcome basis. Columns are non-zero
#' only for the arm's training rows and, under cross-fitting, only for
#' rows outside `i`'s fold. An arm (or a fold's training arm) with fewer
#' than 10 units leaves a zero block (IPW-only for those rows); when
#' `warn = TRUE` this is surfaced via a single `warning()`.
#'
#' @keywords internal
.cme_coef_matrix <- function(y, covariates, arm_idx, n, lambda,
                             cross_fit, n_folds, fold_id, kernel_y,
                             warn = TRUE) {
  C <- matrix(0, n, n)
  if (length(arm_idx) < 10L) {
    if (warn) {
      warning("A treatment arm has fewer than 10 units; using IPW only (zero outcome model) for that arm.",
        call. = FALSE
      )
    }
    return(C)
  }
  if (!cross_fit) {
    cme <- fit_cme(covariates[arm_idx, , drop = FALSE],
      y[arm_idx, , drop = FALSE],
      kernel_y = kernel_y, lambda = lambda
    )
    C[, arm_idx] <- predict(cme, covariates)
    return(C)
  }
  any_fold_short <- FALSE
  for (k in seq_len(n_folds)) {
    test_i <- which(fold_id == k)
    if (length(test_i) == 0L) next
    train_arm <- arm_idx[fold_id[arm_idx] != k]
    if (length(train_arm) < 10L) {
      any_fold_short <- TRUE
      next # rows test_i stay zero -> IPW-only for those rows
    }
    cme <- fit_cme(covariates[train_arm, , drop = FALSE],
      y[train_arm, , drop = FALSE],
      kernel_y = kernel_y, lambda = lambda
    )
    C[test_i, train_arm] <- predict(cme, covariates[test_i, , drop = FALSE])
  }
  if (warn && any_fold_short) {
    warning("Some cross-fitting folds had fewer than 10 arm units; those rows fall back to IPW only.",
      call. = FALSE
    )
  }
  C
}
