# test-end-to-end.R -- whole-package analytical-correctness integration test.
#
# One coherent confounded data-generating process is run through the entire
# causal-inference pipeline, and every analytical corner is checked for a
# meaningful, correct verdict -- not merely that the call returns. The thresholds
# are pinned to a fixed seed and were calibrated empirically: each assertion is
# a statistical claim the method is supposed to satisfy (Type I control under a
# true null, power under a genuine effect, double robustness, low-rank
# approximation fidelity, and full-pipeline reproducibility).

# Shared data-generating process -----------------------------------------------
# The confounder z (two columns) drives both treatment assignment (the backdoor
# path) and the outcome. `effect` toggles a genuine distributional treatment
# effect: a mean shift plus a variance inflation on the treated arm. When
# effect == 0 the only treatment-outcome association is the spurious one routed
# through z, which the backdoor-adjusted methods must remove.

make_confounded_dgp <- function(n, effect, seed) {
  set.seed(seed)
  z1 <- stats::rnorm(n)
  z2 <- stats::rnorm(n)
  e_true <- stats::plogis(0.8 * z1 - 0.6 * z2)
  treatment <- stats::rbinom(n, 1L, e_true)
  sd_treated <- ifelse(treatment == 1L, 1 + 0.8 * effect, 1)
  y <- z1 - 0.5 * z2 + effect * treatment + stats::rnorm(n, sd = sd_treated)
  list(
    y = y,
    treatment = treatment,
    z = cbind(z1, z2),
    e_true = e_true
  )
}

n_perm <- 299L

# Corner 1: marginal independence testing -------------------------------------

test_that("HSIC detects the marginal treatment-outcome association and respects independence", {
  d <- make_confounded_dgp(300L, effect = 1.5, seed = 101L)
  dependent <- hsic_test(d$treatment, d$y, n_permutations = n_perm, seed = 1L)
  independent <- hsic_test(
    stats::rnorm(300L), stats::rnorm(300L),
    n_permutations = n_perm, seed = 1L
  )

  expect_lt(dependent$p_value, 0.05)
  expect_gt(independent$p_value, 0.05)
  # Power separation is large, not borderline.
  expect_gt(independent$p_value, dependent$p_value)
})

# Corner 2: backdoor-adjusted causal association ------------------------------
# The defining claim of bd-HSIC: a purely confounder-induced association is
# removed by the adjustment, while a genuine causal effect survives it. The
# split-sample density-ratio test is low-power, so the genuine regime uses a
# larger n; this is the slow corner and is skipped on CRAN.
#
# Permutation uses explicit confounder-decile clusters with within-cluster
# permutation rather than the default `n_clusters = "auto"`. The "auto" path
# k-means-clusters the propensity scores, and its cluster assignment -- hence
# the permutation null and the p-value -- depends on the platform BLAS;
# deterministic `cluster_id` makes the verdict reproducible across platforms.

# Deterministic exchangeability blocks: deciles of the true propensity.
confounder_blocks <- function(e_true) {
  cut(
    e_true,
    breaks = stats::quantile(e_true, probs = seq(0, 1, length.out = 11L)),
    include.lowest = TRUE, labels = FALSE
  )
}

test_that("bd-HSIC removes a confounded association and detects a genuine causal effect", {
  skip_on_cran()
  genuine <- make_confounded_dgp(800L, effect = 2.5, seed = 101L)
  confounded_only <- make_confounded_dgp(800L, effect = 0, seed = 101L)

  res_genuine <- bd_hsic_test(
    genuine$treatment, genuine$y, genuine$z,
    cluster_id = confounder_blocks(genuine$e_true),
    permutation = "within_cluster", n_permutations = n_perm, seed = 7L
  )
  res_confounded <- bd_hsic_test(
    confounded_only$treatment, confounded_only$y, confounded_only$z,
    cluster_id = confounder_blocks(confounded_only$e_true),
    permutation = "within_cluster", n_permutations = n_perm, seed = 7L
  )

  # Genuine causal effect detected after backdoor adjustment (p at the floor).
  expect_lt(res_genuine$p_value, 0.05)
  # Spurious confounder-induced association removed -- no false rejection.
  expect_gt(res_confounded$p_value, 0.05)
  # The adjustment moves the confounded p-value far above the genuine one.
  expect_gt(res_confounded$p_value, res_genuine$p_value)
})

# Corner 3: propensity recovery -----------------------------------------------

test_that("estimate_propensity recovers the true propensity surface", {
  d <- make_confounded_dgp(800L, effect = 2, seed = 101L)
  fit <- estimate_propensity(d$treatment, d$z, seed = 3L)

  expect_s3_class(fit, "propensity_fit")
  expect_gt(stats::cor(fit$scores, d$e_true), 0.9)
  expect_lt(sqrt(mean((fit$scores - d$e_true)^2)), 0.1)
  # Scores are valid probabilities (and trimmed away from 0/1).
  expect_true(all(fit$scores > 0 & fit$scores < 1))
})

# Corner 4: doubly robust distributional treatment effect (DR-DATE) -----------
# Four cells: power and Type I control, each under the full AIPW estimator and
# under the IPW-only fallback (outcome_model = "zero"). Correct verdicts in all
# four are the operational signature of double robustness.

test_that("DR-DATE has power, controls Type I, and is doubly robust", {
  effect <- make_confounded_dgp(400L, effect = 1.2, seed = 202L)
  null <- make_confounded_dgp(400L, effect = 0, seed = 202L)

  aipw_effect <- dr_date_test(
    effect$y, effect$treatment, effect$z,
    n_permutations = n_perm, seed = 11L
  )
  aipw_null <- dr_date_test(
    null$y, null$treatment, null$z,
    n_permutations = n_perm, seed = 11L
  )
  ipw_effect <- dr_date_test(
    effect$y, effect$treatment, effect$z,
    outcome_model = "zero", n_permutations = n_perm, seed = 11L
  )
  ipw_null <- dr_date_test(
    null$y, null$treatment, null$z,
    outcome_model = "zero", n_permutations = n_perm, seed = 11L
  )

  # AIPW: power under the alternative, Type I control under the null.
  expect_lt(aipw_effect$p_value, 0.05)
  expect_gt(aipw_null$p_value, 0.05)
  # IPW-only fallback: still valid with the outcome model switched off.
  expect_lt(ipw_effect$p_value, 0.05)
  expect_gt(ipw_null$p_value, 0.05)
})

# Corner 5: effect on the treated (DR-DETT) -----------------------------------

test_that("DR-DETT detects an effect on the treated and controls Type I", {
  effect <- make_confounded_dgp(400L, effect = 1.2, seed = 202L)
  null <- make_confounded_dgp(400L, effect = 0, seed = 202L)

  res_effect <- dr_dett_test(
    effect$y, effect$treatment, effect$z,
    n_permutations = n_perm, seed = 13L
  )
  res_null <- dr_dett_test(
    null$y, null$treatment, null$z,
    n_permutations = n_perm, seed = 13L
  )

  expect_lt(res_effect$p_value, 0.05)
  expect_gt(res_null$p_value, 0.05)
})

# Corner 6: two-sample MMD ----------------------------------------------------

test_that("MMD distinguishes different distributions and respects equality", {
  set.seed(55L)
  a <- stats::rnorm(150L)
  b_same <- stats::rnorm(150L)
  b_diff <- stats::rnorm(150L, mean = 1)

  different <- mmd_test(a, b_diff, n_permutations = n_perm, seed = 5L)
  same <- mmd_test(a, b_same, n_permutations = n_perm, seed = 5L)

  expect_lt(different$p_value, 0.05)
  expect_gt(same$p_value, 0.05)
})

# Corner 7: low-rank approximation fidelity -----------------------------------
# Nystrom and RFF accelerations must reproduce the exact-HSIC verdict and track
# its statistic closely on the same data.

test_that("Nystrom and RFF HSIC agree with the exact statistic", {
  set.seed(77L)
  x <- stats::rnorm(400L)
  y <- 0.7 * x + stats::rnorm(400L)

  exact <- hsic_test(x, y, n_permutations = n_perm, seed = 9L)
  nystrom <- hsic_test_nystrom(
    x, y, method = "nystrom", m = 80L, n_permutations = n_perm, seed = 9L
  )
  rff <- hsic_test_nystrom(
    x, y, method = "rff", m = 80L, n_permutations = n_perm, seed = 9L
  )

  # Same verdict as the exact test.
  expect_lt(exact$p_value, 0.05)
  expect_lt(nystrom$p_value, 0.05)
  expect_lt(rff$p_value, 0.05)
  # Statistic fidelity: both approximations track exact to within 10%.
  expect_lt(abs(nystrom$statistic - exact$statistic) / exact$statistic, 0.1)
  expect_lt(abs(rff$statistic - exact$statistic) / exact$statistic, 0.1)
})

# Corner 8: hierarchical within-cluster permutation ---------------------------
# A clustered null with strong cluster random effects but no treatment effect:
# within-cluster permutation must preserve the cluster structure under the null
# and not manufacture a rejection.

test_that("hierarchical within-cluster permutation controls Type I under a clustered null", {
  skip_on_cran()
  set.seed(88L)
  n_clusters <- 10L
  per_cluster <- 40L
  n <- n_clusters * per_cluster
  cluster_id <- rep(seq_len(n_clusters), each = per_cluster)
  cluster_effect <- rep(stats::rnorm(n_clusters, sd = 1.5), each = per_cluster)
  z <- cbind(stats::rnorm(n), stats::rnorm(n))
  treatment <- stats::rbinom(n, 1L, stats::plogis(0.5 * z[, 1]))
  # Outcome carries the cluster effect and a confounder, but no treatment term.
  y <- cluster_effect + 0.8 * z[, 1] + stats::rnorm(n)

  res <- hierarchical_test(
    y, treatment, z,
    cluster_id = cluster_id, n_permutations = n_perm, seed = 21L
  )

  expect_gt(res$p_value, 0.05)
  expect_true(is.finite(res$statistic))
})

# Corner 9: full-pipeline reproducibility -------------------------------------
# Identical seeds must reproduce the statistic and p-value bit-for-bit across
# the permutation machinery; the verdict must be stable across seeds.

test_that("seeded results are reproducible and verdicts are seed-stable", {
  d <- make_confounded_dgp(400L, effect = 1.2, seed = 202L)

  first <- dr_date_test(d$y, d$treatment, d$z, n_permutations = 99L, seed = 999L)
  second <- dr_date_test(d$y, d$treatment, d$z, n_permutations = 99L, seed = 999L)
  expect_identical(first$statistic, second$statistic)
  expect_identical(first$p_value, second$p_value)

  set.seed(77L)
  x <- stats::rnorm(400L)
  y <- 0.7 * x + stats::rnorm(400L)
  hsic_a <- hsic_test(x, y, n_permutations = 99L, seed = 5L)
  hsic_b <- hsic_test(x, y, n_permutations = 99L, seed = 6L)
  # Different seeds, same scientific conclusion.
  expect_true(hsic_a$p_value < 0.05 && hsic_b$p_value < 0.05)
})

# Corner 10: permutation-null calibration -------------------------------------
# Under a sequence of true nulls the permutation p-values must not over-reject.
# Heavier than a single test; skipped on CRAN.

test_that("DR-DATE permutation p-values are calibrated under the null", {
  skip_on_cran()
  p_values <- vapply(
    seq_len(12L),
    function(i) {
      d <- make_confounded_dgp(300L, effect = 0, seed = 500L + i)
      dr_date_test(
        d$y, d$treatment, d$z,
        n_permutations = 199L, seed = i
      )$p_value
    },
    numeric(1)
  )
  # A valid (typically conservative) permutation test rejects well below the
  # Markov bound; a grossly inflated test would blow past it.
  expect_lt(mean(p_values < 0.05), 0.25)
})
