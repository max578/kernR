#' Doubly Robust Distributional Effect on the Treated Test (DR-DETT)
#'
#' Tests whether the distribution of the treated potential outcome Y(1)
#' differs from the control potential outcome Y(0) among the treated
#' subpopulation. Requires only one-sided overlap.
#'
#' @inheritParams dr_date_test
#'
#' @return An object of class `"kernel_test_result"`. The `ess` element
#'   holds the effective sample size of the control reconstruction
#'   weights and `ess_warning` records whether the reliability floor was
#'   hit.
#'
#' @details
#' DR-DETT is analogous to DR-DATE but focuses on the **effect on the
#' treated** (ETT) rather than the average treatment effect. The treated
#' counterfactual Y(1) | T = 1 is observed directly, so only the
#' *control* arm needs an outcome model. The control counterfactual
#' Y(0) | T = 1 is reconstructed by augmented inverse-probability
#' weighting, reweighting controls by the treatment odds
#' \eqn{e(x) / (1 - e(x))} so that their covariate distribution matches
#' the treated. With `outcome_model = "krr"` the control conditional mean
#' embedding \eqn{\hat m_0} supplies the doubly robust augmentation;
#' `outcome_model = "zero"` returns the inverse-probability-weighted
#' statistic. It requires only one-sided overlap: P(T = 1 | X) bounded
#' away from 0 (not necessarily from 1), so it applies where positivity
#' fails for the control group.
#'
#' @references
#' Fawkes, J., Hu, R., Evans, R. J., & Sejdinovic, D. (2024). Doubly
#' robust kernel statistics for testing distributional treatment effects.
#' *Transactions on Machine Learning Research*.
#'
#' @examples
#' set.seed(42)
#' n <- 300
#' x <- matrix(rnorm(n * 2), n, 2)
#' t <- rbinom(n, 1, plogis(0.5 * x[, 1]))
#' y <- t * rnorm(n, sd = 2) + (1 - t) * rnorm(n, sd = 1) + x[, 1]
#'
#' result <- dr_dett_test(y, t, x, n_permutations = 200, seed = 1)
#' print(result)
#'
#' @export
dr_dett_test <- function(y, treatment, covariates,
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
  ps <- estimate_propensity(treatment, covariates,
    method = propensity_model,
    cross_fit = cross_fit, n_folds = n_folds
  )
  e_hat <- ps$scores

  kernel_y <- resolve_bandwidth(kernel_y, y)
  Ky <- kernel_matrix(y, kernel = kernel_y)

  # Only the control arm needs an outcome model for the ETT estimand
  # (treated outcomes are observed). Cross-fitted when cross_fit = TRUE.
  warn_reliability <- n_permutations > 0L
  if (outcome_model == "krr") {
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
    C0 <- .cme_coef_matrix(y, covariates, which(treatment == 0L), n,
      lambda, cross_fit, n_folds, fold_id, kernel_y, warn_reliability)
  } else {
    C0 <- matrix(0, n, n)
  }

  # ESS reliability gate on the treatment-odds control weights, relative
  # to the treated count (the ETT compares within the treated).
  n1 <- sum(treatment)
  ess <- effective_sample_size(((1 - treatment) * e_hat / (1 - e_hat))[treatment == 0L])
  ess_warning <- ess < min_ess_fraction * n1
  if (warn_reliability && ess_warning) {
    warning(sprintf(
      paste0(
        "Control-reconstruction effective sample size %.1f (%.1f%% of ",
        "n_treated) is below the %.0f%% floor; the treatment-odds ",
        "weights are concentrated and the test may be unreliable."
      ),
      ess, 100 * ess / n1, 100 * min_ess_fraction
    ), call. = FALSE)
  }

  # Compute DR-DETT statistic
  stat_obs <- compute_dr_dett_stat(Ky, treatment, e_hat, C0)

  # Permutation null
  null_dist <- numeric(n_permutations)
  for (p in seq_len(n_permutations)) {
    t_perm <- bin_permute_treatment(treatment, e_hat, n_bins)
    null_dist[p] <- compute_dr_dett_stat(Ky, t_perm, e_hat, C0)
  }

  p_value <- (1 + sum(null_dist >= stat_obs)) / (1 + n_permutations)

  structure(
    list(
      statistic = stat_obs,
      p_value = p_value,
      method = "DR-DETT",
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

#' Compute DR-DETT Test Statistic (augmented IPW, effect on the treated)
#'
#' @param Ky n x n outcome kernel matrix.
#' @param treatment Binary treatment vector.
#' @param e_hat Propensity scores.
#' @param C0 n x n control-arm conditional-mean-embedding coefficient
#'   matrix (zero matrix gives the inverse-probability-weighted statistic).
#'
#' @return Scalar test statistic.
#' @keywords internal
compute_dr_dett_stat <- function(Ky, treatment, e_hat, C0) {
  n <- length(treatment)
  t <- treatment
  n1 <- sum(t)

  if (n1 < 2) {
    return(0)
  }

  # Treated arm Y(1) | T = 1 is observed: uniform over the treated.
  beta1 <- t / n1

  # Control counterfactual Y(0) | T = 1 by augmented IPW, reweighting
  # controls by the treatment odds e / (1 - e) toward the treated
  # covariate distribution:
  #   beta_0 = (1/n1) [ u + C0^T (t - u) ],  u_i = (1 - t_i) e_i / (1 - e_i).
  # With C0 = 0 this is the inverse-probability-weighted ETT statistic.
  u <- (1 - t) * e_hat / (1 - e_hat)
  beta0 <- (u + as.numeric(crossprod(C0, t - u))) / n1

  d <- beta1 - beta0
  mmd2 <- as.numeric(crossprod(d, Ky %*% d))
  max(mmd2, 0)
}
