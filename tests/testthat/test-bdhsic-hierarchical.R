make_clustered_design <- function(seed = 1L, n_clust = 6L, n_per = 30L,
                                  causal_effect = 0.6, cluster_sd = 1.2) {
  set.seed(seed)
  n <- n_clust * n_per
  cluster_id <- rep(seq_len(n_clust), each = n_per)
  cluster_effect <- stats::rnorm(n_clust, sd = cluster_sd)[cluster_id]
  z <- matrix(stats::rnorm(n * 2L), n, 2L)
  x <- z[, 1L] + 0.4 * cluster_effect + stats::rnorm(n)
  y <- causal_effect * x + z[, 2L] + cluster_effect + stats::rnorm(n, sd = 0.5)
  list(x = x, y = y, z = z,
       cluster_id = factor(paste0("site", cluster_id)))
}


test_that("cluster_id activates within-cluster permutation by default", {
  d <- make_clustered_design(seed = 11L)
  res <- bd_hsic_test(
    d$x, d$y, d$z,
    cluster_id     = d$cluster_id,
    n_permutations = 99L,
    seed           = 11L
  )
  expect_s3_class(res, "kernel_test_result")
  expect_equal(res$method, "bd-HSIC")
  expect_equal(res$permutation_scheme, "within_cluster")
  expect_true(is.integer(res$cluster_id))
  expect_true(is.character(res$cluster_levels))
  expect_true(length(res$cluster_levels) >= 2L)
  expect_true(is.numeric(res$per_cluster_statistic))
  expect_equal(length(res$per_cluster_statistic),
               length(res$cluster_levels))
})

test_that("backwards compatible: cluster_id = NULL preserves original behaviour", {
  d <- make_clustered_design(seed = 12L)
  res <- bd_hsic_test(
    d$x, d$y, d$z,
    n_permutations = 99L,
    seed           = 12L
  )
  expect_equal(res$permutation_scheme, "propensity")
  expect_null(res$cluster_id)
  expect_null(res$per_cluster_statistic)
})

test_that("permutation = 'naive' works with and without cluster_id", {
  d <- make_clustered_design(seed = 13L)

  res_no_cl <- bd_hsic_test(
    d$x, d$y, d$z,
    permutation    = "naive",
    n_permutations = 99L,
    seed           = 13L
  )
  expect_equal(res_no_cl$permutation_scheme, "naive")

  res_with_cl <- bd_hsic_test(
    d$x, d$y, d$z,
    cluster_id     = d$cluster_id,
    permutation    = "naive",
    n_permutations = 99L,
    seed           = 13L
  )
  expect_equal(res_with_cl$permutation_scheme, "naive")
  # per_cluster_statistic still reported under naive when cluster_id supplied
  expect_true(is.numeric(res_with_cl$per_cluster_statistic))
})

test_that("permutation = 'within_cluster' requires cluster_id", {
  d <- make_clustered_design(seed = 14L)
  expect_error(
    bd_hsic_test(d$x, d$y, d$z,
                 permutation = "within_cluster",
                 n_permutations = 50L, seed = 14L),
    "requires"
  )
})

test_that("cluster_id length must match n", {
  d <- make_clustered_design(seed = 15L)
  bad_cl <- d$cluster_id[seq_len(length(d$cluster_id) - 1L)]
  expect_error(
    bd_hsic_test(d$x, d$y, d$z, cluster_id = bad_cl,
                 n_permutations = 50L, seed = 15L),
    "length"
  )
})

test_that("cluster_id must define at least 2 distinct clusters", {
  d <- make_clustered_design(seed = 16L)
  flat <- factor(rep("only_site", length(d$x)))
  expect_error(
    bd_hsic_test(d$x, d$y, d$z, cluster_id = flat,
                 n_permutations = 50L, seed = 16L),
    "at least 2"
  )
})

test_that("seed reproducibility holds under within-cluster permutation", {
  d <- make_clustered_design(seed = 17L)
  r1 <- bd_hsic_test(d$x, d$y, d$z, cluster_id = d$cluster_id,
                     n_permutations = 99L, seed = 17L)
  r2 <- bd_hsic_test(d$x, d$y, d$z, cluster_id = d$cluster_id,
                     n_permutations = 99L, seed = 17L)
  expect_equal(r1$statistic, r2$statistic)
  expect_equal(r1$p_value, r2$p_value)
  expect_equal(r1$per_cluster_statistic, r2$per_cluster_statistic)
})

test_that("under H0 (no causal effect) within-cluster permutation does not reject", {
  # Cluster effects exist but x has no causal effect on y
  set.seed(21L)
  n_clust <- 6L; n_per <- 25L; n <- n_clust * n_per
  cl_id <- rep(seq_len(n_clust), each = n_per)
  cl_eff <- stats::rnorm(n_clust, sd = 1.2)[cl_id]
  z <- matrix(stats::rnorm(n * 2L), n, 2L)
  x <- z[, 1L] + 0.4 * cl_eff + stats::rnorm(n)
  y <- z[, 2L] + cl_eff + stats::rnorm(n, sd = 0.5)  # x has no effect

  res <- bd_hsic_test(
    x, y, z,
    cluster_id     = factor(paste0("c", cl_id)),
    n_permutations = 199L,
    seed           = 21L
  )
  expect_gt(res$p_value, 0.05)
})
