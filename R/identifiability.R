#' Latin-Hypercube Design Over Bounded Parameters
#'
#' Generates a Latin-hypercube sample of size `n` over the parameter
#' bounds in `bounds`. Each column is a random permutation of stratified
#' uniform draws (one per equal-width bin in `(0, 1]`), then scaled to the
#' supplied parameter range. The result is reproducible when `seed` is
#' supplied.
#'
#' This is a lightweight helper aimed at pre-PESTO screening: produce a
#' design matrix to feed an APSIM (or any) simulator, then pass the
#' resulting input/output pairs to [hsic_identifiability()] to flag
#' unidentifiable parameters before ensemble-smoother calibration.
#'
#' @param n Integer. Number of design points (rows). Must be `>= 2`.
#' @param bounds Two-column numeric matrix or data.frame of `[lower, upper]`
#'   bounds, with one row per parameter. If named, row names propagate to
#'   the column names of the returned design.
#' @param seed Integer or `NULL`. Random seed for reproducibility.
#'
#' @return A numeric matrix of dimension `n x nrow(bounds)`. Column names
#'   are inherited from `rownames(bounds)` when available, otherwise
#'   `theta1`, `theta2`, ....
#'
#' @references
#' McKay, M. D., Beckman, R. J., & Conover, W. J. (1979). A comparison of
#' three methods for selecting values of input variables in the analysis of
#' output from a computer code. *Technometrics*, 21(2), 239-245.
#'
#' @examples
#' bounds <- rbind(
#'   slope     = c(0.1, 2.0),
#'   intercept = c(-1, 1),
#'   noise_sd  = c(0.05, 0.5)
#' )
#' design <- lhs_design(50, bounds, seed = 1)
#' head(design)
#'
#' @seealso [hsic_identifiability()]
#' @export
lhs_design <- function(n, bounds, seed = NULL) {
  n <- as.integer(n)
  if (length(n) != 1L || is.na(n) || n < 2L) {
    stop("`n` must be an integer >= 2.", call. = FALSE)
  }

  bounds <- as.matrix(bounds)
  if (ncol(bounds) != 2L) {
    stop("`bounds` must have exactly two columns: [lower, upper].",
         call. = FALSE)
  }
  if (any(!is.finite(bounds))) {
    stop("`bounds` must contain only finite numeric values.", call. = FALSE)
  }
  if (any(bounds[, 2L] <= bounds[, 1L])) {
    stop("Each upper bound must be strictly greater than its lower bound.",
         call. = FALSE)
  }

  p <- nrow(bounds)
  param_names <- rownames(bounds)
  if (is.null(param_names)) param_names <- paste0("theta", seq_len(p))

  if (!is.null(seed)) set.seed(seed)

  design <- matrix(NA_real_, nrow = n, ncol = p,
                   dimnames = list(NULL, param_names))
  for (j in seq_len(p)) {
    u <- (sample.int(n) - stats::runif(n)) / n
    design[, j] <- bounds[j, 1L] + u * (bounds[j, 2L] - bounds[j, 1L])
  }
  design
}


#' HSIC-Based Identifiability Diagnostic
#'
#' Pre-PESTO (or pre-IES) screening: for each parameter `theta[, j]` and
#' each output `y[, k]`, computes an HSIC permutation test of independence
#' and flags parameters with no detectable association to any output as
#' unidentifiable. Useful for trimming the parameter space before
#' ensemble-smoother calibration of mechanistic ag-system models such as
#' APSIM.
#'
#' Kernel matrices are computed once per parameter and once per output, so
#' the total cost is `O((p + q) n^2)` (kernel construction) plus
#' `O(p q n_permutations n^2)` (permutation null), where `p` is the number
#' of parameters, `q` the number of outputs and `n` the design size.
#'
#' A parameter is **identifiable** at level `alpha` when its smallest
#' (optionally adjusted) p-value across outputs satisfies
#' `min_p <= alpha`. Across-grid p-value adjustment defaults to
#' Benjamini-Hochberg, which is the natural FDR control for screening
#' applications.
#'
#' @param theta Numeric matrix `n x p` of parameter design points (one row
#'   per simulator run, one column per parameter). Vectors are coerced via
#'   [as.matrix()].
#' @param y Numeric matrix `n x q` of simulator outputs. Vectors are
#'   coerced via [as.matrix()].
#' @param alpha Numeric in `(0, 1)`. Identifiability threshold on the
#'   adjusted minimum p-value. Default `0.05`.
#' @param p_adjust Character. Across-grid p-value adjustment method passed
#'   to [stats::p.adjust()]. Default `"BH"`. Use `"none"` to disable.
#' @param n_permutations Integer. Permutations per HSIC test. Default 500.
#' @param kernel_theta A [kernel_spec()] for parameter columns. Default
#'   RBF with per-column median heuristic.
#' @param kernel_y A [kernel_spec()] for output columns. Default RBF with
#'   per-column median heuristic.
#' @param seed Integer or `NULL`. Random seed for reproducibility.
#'
#' @return An object of class `"hsic_identifiability"` with components:
#'   \describe{
#'     \item{statistic}{`p x q` matrix of HSIC statistics.}
#'     \item{p_value}{`p x q` matrix of raw permutation p-values.}
#'     \item{p_value_adjusted}{`p x q` matrix of adjusted p-values (same as
#'       `p_value` when `p_adjust = "none"`).}
#'     \item{max_statistic}{Length-`p` vector: per-parameter maximum HSIC
#'       across outputs.}
#'     \item{min_p_value}{Length-`p` vector: per-parameter minimum adjusted
#'       p-value across outputs.}
#'     \item{identifiable}{Length-`p` logical: `min_p_value <= alpha`.}
#'     \item{rank}{Parameter indices ordered by descending `max_statistic`.}
#'     \item{alpha, p_adjust, n, n_permutations}{Inputs / metadata.}
#'     \item{param_names, output_names}{Character vectors.}
#'     \item{call}{The matched call.}
#'   }
#'
#' @references
#' Gretton, A., Fukumizu, K., Teo, C. H., Song, L., Scholkopf, B., &
#' Smola, A. J. (2008). A kernel statistical test of independence.
#' *NeurIPS*, 20.
#'
#' Da Veiga, S. (2015). Global sensitivity analysis with dependence
#' measures. *Journal of Statistical Computation and Simulation*, 85(7),
#' 1283-1305.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 60
#' # 3 active parameters + 1 inert
#' theta <- matrix(stats::runif(n * 4), nrow = n,
#'                 dimnames = list(NULL, paste0("p", 1:4)))
#' y1 <- theta[, 1] + 0.5 * theta[, 2]^2 + stats::rnorm(n, sd = 0.1)
#' y2 <- sin(2 * pi * theta[, 3]) + stats::rnorm(n, sd = 0.1)
#' y <- cbind(yield = y1, biomass = y2)
#' fit <- hsic_identifiability(theta, y, n_permutations = 199, seed = 1)
#' print(fit)
#' }
#'
#' @seealso [lhs_design()], [hsic_test()]
#' @export
hsic_identifiability <- function(theta, y,
                                 alpha = 0.05,
                                 p_adjust = c("BH", "holm", "hochberg",
                                              "bonferroni", "BY", "fdr",
                                              "none"),
                                 n_permutations = 500L,
                                 kernel_theta = kernel_spec(),
                                 kernel_y = kernel_spec(),
                                 seed = NULL) {
  cl <- match.call()
  p_adjust <- match.arg(p_adjust)

  theta <- as.matrix(theta)
  y <- as.matrix(y)
  n <- nrow(theta)

  if (nrow(y) != n) {
    stop("`theta` and `y` must have the same number of rows.", call. = FALSE)
  }
  if (n < 10L) {
    stop("At least 10 observations (design points) are required.",
         call. = FALSE)
  }
  if (!is.numeric(alpha) || length(alpha) != 1L ||
      alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a single number in (0, 1).", call. = FALSE)
  }
  n_permutations <- as.integer(n_permutations)
  if (length(n_permutations) != 1L || is.na(n_permutations) ||
      n_permutations < 1L) {
    stop("`n_permutations` must be a positive integer.", call. = FALSE)
  }

  p <- ncol(theta)
  q <- ncol(y)
  param_names <- colnames(theta)
  if (is.null(param_names)) param_names <- paste0("theta", seq_len(p))
  output_names <- colnames(y)
  if (is.null(output_names)) output_names <- paste0("y", seq_len(q))

  if (!is.null(seed)) set.seed(seed)

  # Precompute kernel matrices column-wise (re-used across the p x q grid).
  Kx_list <- vector("list", p)
  for (j in seq_len(p)) {
    xj <- theta[, j, drop = FALSE]
    kj <- resolve_bandwidth(kernel_theta, xj)
    Kx_list[[j]] <- kernel_matrix(xj, kernel = kj)
  }
  Ky_list <- vector("list", q)
  for (k in seq_len(q)) {
    yk <- y[, k, drop = FALSE]
    kk <- resolve_bandwidth(kernel_y, yk)
    Ky_list[[k]] <- kernel_matrix(yk, kernel = kk)
  }

  stat_mat <- matrix(NA_real_, p, q,
                     dimnames = list(param_names, output_names))
  pval_mat <- stat_mat

  for (j in seq_len(p)) {
    for (k in seq_len(q)) {
      stat <- hsic_stat_cpp(Kx_list[[j]], Ky_list[[k]])
      nd   <- permutation_hsic_cpp(Kx_list[[j]], Ky_list[[k]],
                                   n_permutations)
      stat_mat[j, k] <- stat
      pval_mat[j, k] <- (1 + sum(nd >= stat)) / (1 + n_permutations)
    }
  }

  pval_adj <- if (p_adjust == "none") {
    pval_mat
  } else {
    matrix(stats::p.adjust(as.vector(pval_mat), method = p_adjust),
           nrow = p, ncol = q,
           dimnames = dimnames(pval_mat))
  }

  max_stat <- apply(stat_mat, 1L, max)
  min_p    <- apply(pval_adj, 1L, min)
  identifiable <- min_p <= alpha
  rank <- order(max_stat, decreasing = TRUE)

  structure(
    list(
      statistic        = stat_mat,
      p_value          = pval_mat,
      p_value_adjusted = pval_adj,
      max_statistic    = max_stat,
      min_p_value      = min_p,
      identifiable     = identifiable,
      rank             = rank,
      alpha            = alpha,
      p_adjust         = p_adjust,
      n                = n,
      n_permutations   = n_permutations,
      param_names      = param_names,
      output_names     = output_names,
      call             = cl
    ),
    class = "hsic_identifiability"
  )
}


#' @export
print.hsic_identifiability <- function(x, digits = 4L, ...) {
  cat("\n  HSIC Identifiability Scan\n\n")
  cat("Parameters:  ", length(x$param_names), "\n")
  cat("Outputs:     ", length(x$output_names), "\n")
  cat("N:           ", x$n, "\n")
  cat("Permutations:", x$n_permutations, "\n")
  cat("Alpha:       ", x$alpha, "\n")
  cat("P-adjust:    ", x$p_adjust, "\n\n")

  ident   <- x$identifiable
  yes_nms <- x$param_names[ident]
  no_nms  <- x$param_names[!ident]

  cat("Identifiable (", length(yes_nms), "): ",
      if (length(yes_nms)) paste(yes_nms, collapse = ", ") else "<none>",
      "\n", sep = "")
  cat("Not identifiable (", length(no_nms), "): ",
      if (length(no_nms)) paste(no_nms, collapse = ", ") else "<none>",
      "\n\n", sep = "")

  ord <- x$rank
  tab <- data.frame(
    parameter    = x$param_names[ord],
    max_HSIC     = formatC(x$max_statistic[ord], digits = digits,
                           format = "g"),
    min_p        = formatC(x$min_p_value[ord], digits = digits,
                           format = "f"),
    identifiable = ifelse(x$identifiable[ord], "*", ""),
    stringsAsFactors = FALSE
  )
  cat("Per-parameter ranking (descending max HSIC):\n")
  print(tab, row.names = FALSE)
  cat("\n  (* = identifiable at alpha =", x$alpha, ")\n\n")
  invisible(x)
}


#' Plot an HSIC Identifiability Scan
#'
#' Bar plot of per-parameter maximum HSIC across outputs, ordered by
#' magnitude. Bars for identifiable parameters are coloured; non-
#' identifiable parameters are shown in grey.
#'
#' @param x An `hsic_identifiability` object.
#' @param col_yes,col_no Bar colours for identifiable / non-identifiable
#'   parameters.
#' @param ... Additional arguments passed to [graphics::barplot()].
#'
#' @return Invisibly returns `x`. Side effect: produces a base R plot.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 50
#' theta <- matrix(stats::runif(n * 3), nrow = n,
#'                 dimnames = list(NULL, c("active", "active2", "inert")))
#' y <- theta[, 1] + theta[, 2]^2 + stats::rnorm(n, sd = 0.1)
#' fit <- hsic_identifiability(theta, y, n_permutations = 199, seed = 1)
#' plot(fit)
#' }
#'
#' @export
plot.hsic_identifiability <- function(x,
                                      col_yes = "#0072B2",
                                      col_no  = "grey70",
                                      ...) {
  ord <- x$rank
  heights <- x$max_statistic[ord]
  cols <- ifelse(x$identifiable[ord], col_yes, col_no)
  graphics::barplot(
    heights,
    names.arg = x$param_names[ord],
    col       = cols,
    border    = NA,
    main      = "HSIC identifiability scan",
    ylab      = "max HSIC across outputs",
    las       = 2L,
    ...
  )
  graphics::legend(
    "topright",
    legend = c(
      paste0("identifiable (alpha=", x$alpha, ")"),
      "not identifiable"
    ),
    fill   = c(col_yes, col_no),
    border = NA,
    bty    = "n"
  )
  invisible(x)
}


#' @export
summary.hsic_identifiability <- function(object, ...) {
  cat("\n  HSIC Identifiability Scan - Summary\n\n")
  cat("Parameters:  ", length(object$param_names), "\n")
  cat("Outputs:     ", length(object$output_names), "\n")
  cat("N:           ", object$n, "\n")
  cat("Permutations:", object$n_permutations, "\n")
  cat("Alpha:       ", object$alpha, "\n")
  cat("P-adjust:    ", object$p_adjust, "\n\n")

  cat("HSIC statistics (parameter x output):\n")
  print(round(object$statistic, 4L))
  cat("\nAdjusted p-values (parameter x output):\n")
  print(round(object$p_value_adjusted, 4L))
  cat("\n")
  invisible(object)
}


#' @export
as.data.frame.hsic_identifiability <- function(x, row.names = NULL,
                                               optional = FALSE, ...) {
  p <- length(x$param_names)
  q <- length(x$output_names)
  data.frame(
    parameter        = rep(x$param_names, times = q),
    output           = rep(x$output_names, each  = p),
    statistic        = as.vector(x$statistic),
    p_value          = as.vector(x$p_value),
    p_value_adjusted = as.vector(x$p_value_adjusted),
    stringsAsFactors = FALSE
  )
}
