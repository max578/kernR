test_that("hierarchical test runs with DR-DATE", {
  set.seed(42)
  n_clusters <- 10
  n_per <- 30
  n <- n_clusters * n_per
  cluster_id <- rep(1:n_clusters, each = n_per)

  x <- matrix(rnorm(n * 2), n, 2)
  t <- rbinom(n, 1, 0.5)
  y <- 0.5 * t + rnorm(n_clusters, sd = 0.5)[cluster_id] + x[, 1] + rnorm(n)

  result <- hierarchical_test(y, t, x, cluster_id,
    method = "dr-date",
    n_permutations = 30,
    seed = 1
  )

  expect_s3_class(result, "kernel_test_result")
  expect_true(grepl("Hierarchical", result$method))
  expect_true(!is.null(result$hierarchical))
  expect_equal(result$hierarchical$n_clusters, n_clusters)
})

test_that("hierarchical test respects weight methods", {
  set.seed(42)
  n_clusters <- 10
  n_per <- 30 # DR within-cluster sub-tests need >= 30 obs per cluster
  n <- n_clusters * n_per
  cluster_id <- rep(seq_len(n_clusters), each = n_per)
  x <- matrix(rnorm(n * 2), n, 2)
  t <- rbinom(n, 1, 0.5)
  y <- rnorm(n)

  r1 <- hierarchical_test(y, t, x, cluster_id,
    weight_method = "equal",
    n_permutations = 20,
    seed = 1
  )
  r2 <- hierarchical_test(y, t, x, cluster_id,
    weight_method = "within_only",
    n_permutations = 20,
    seed = 1
  )

  expect_s3_class(r1, "kernel_test_result")
  expect_s3_class(r2, "kernel_test_result")
  # The within-cluster component is actually computed -- guards against the
  # masked-failure bug where every within sub-test silently NA'd out.
  expect_false(all(is.na(r1$hierarchical$within_stats)))
})

test_that("hierarchical test rejects with too few clusters", {
  expect_error(
    hierarchical_test(rnorm(40), rbinom(40, 1, 0.5),
      matrix(rnorm(80), 40, 2),
      rep(1:2, each = 20),
      n_permutations = 10
    ),
    "3 clusters"
  )
})
