#' @export
print.kernel_test_result <- function(x, ...) {
  cat("\n")
  cat("  ", x$method, "Test\n\n")
  cat("Statistic:", formatC(x$statistic, digits = 6, format = "g"), "\n")
  cat("P-value:  ", formatC(x$p_value, digits = 4, format = "f"), "\n")
  cat("N:        ", x$n, "\n")
  cat("Perms:    ", x$n_permutations, "\n")
  if (!is.null(x$kernel_x)) {
    cat("Kernel X: ", x$kernel_x$type)
    if (x$kernel_x$type %in% c("rbf", "matern") && is.numeric(x$kernel_x$bandwidth)) {
      cat(" (bw = ", formatC(x$kernel_x$bandwidth, digits = 4, format = "g"), ")", sep = "")
    }
    cat("\n")
  }
  if (!is.null(x$kernel_y)) {
    cat("Kernel Y: ", x$kernel_y$type)
    if (x$kernel_y$type %in% c("rbf", "matern") && is.numeric(x$kernel_y$bandwidth)) {
      cat(" (bw = ", formatC(x$kernel_y$bandwidth, digits = 4, format = "g"), ")", sep = "")
    }
    cat("\n")
  }
  if (!is.na(x$ess)) {
    cat("ESS:      ", formatC(x$ess, digits = 1, format = "f"), "\n")
  }
  cat("\n")
  invisible(x)
}

#' @export
summary.kernel_test_result <- function(object, ...) {
  cat("\n")
  cat("  ", object$method, "Test - Summary\n\n")
  cat("Statistic:    ", formatC(object$statistic, digits = 6, format = "g"), "\n")
  cat("P-value:      ", formatC(object$p_value, digits = 4, format = "f"), "\n")
  cat("Sample size:  ", object$n, "\n")
  cat("Permutations: ", object$n_permutations, "\n")

  if (!is.na(object$ess)) {
    cat("Effective SS: ", formatC(object$ess, digits = 1, format = "f"), "\n")
  }

  # Null distribution summary
  nd <- object$null_distribution
  if (length(nd) > 0) {
    cat("\nNull distribution:\n")
    cat("  Min:    ", formatC(min(nd), digits = 4, format = "g"), "\n")
    cat("  Median: ", formatC(median(nd), digits = 4, format = "g"), "\n")
    cat("  Max:    ", formatC(max(nd), digits = 4, format = "g"), "\n")
    cat("  SD:     ", formatC(sd(nd), digits = 4, format = "g"), "\n")
  }

  cat("\nCall: ")
  print(object$call)
  cat("\n")
  invisible(object)
}

#' Plot a Kernel Test Result
#'
#' Plots the permutation null distribution with the observed statistic.
#'
#' @param x A `kernel_test_result` object.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisibly returns `x`. Side effect: produces a base R plot.
#'
#' @examples
#' set.seed(42)
#' x_data <- rnorm(100)
#' y_data <- x_data + rnorm(100, sd = 0.5)
#' res <- hsic_test(x_data, y_data)
#' plot(res)
#'
#' @export
plot.kernel_test_result <- function(x, ...) {
  nd <- x$null_distribution
  obs <- x$statistic

  hist(nd,
    breaks = 30,
    main = paste(x$method, "Permutation Null Distribution"),
    xlab = "Test Statistic",
    col = "grey80",
    border = "grey60",
    freq = FALSE
  )
  abline(v = obs, col = "#D55E00", lwd = 2, lty = 2)
  legend("topright",
    legend = c(
      paste("Observed =", formatC(obs, digits = 4, format = "g")),
      paste("p =", formatC(x$p_value, digits = 4, format = "f"))
    ),
    col = c("#D55E00", NA),
    lwd = c(2, NA),
    lty = c(2, NA),
    bty = "n"
  )

  invisible(x)
}

#' @export
print.taci_result <- function(x, ...) {
  cat("TACI mechanism-consistency test\n")
  if (!is.null(x$treatment_type)) {
    cat(sprintf("  treatment: %s%s\n", x$treatment_type,
                if (identical(x$treatment_type, "continuous")) {
                  sprintf(" (H0 baseline = %.4g)", x$baseline)
                } else {
                  ""
                }))
  }
  cat(sprintf("  statistic: %s bd-HSIC%s\n",
              if (isTRUE(x$adjusted)) "backdoor-ADJUSTED" else "unadjusted",
              if (isTRUE(x$adjusted)) {
                sprintf(" (density_ratio = %s)", x$density_ratio)
              } else {
                ""
              }))
  cat(sprintf("  observed bd-HSIC: %.4g\n", x$observed_statistic))
  cat(sprintf("  H0 tail p-value:  %.3f  (in tail: %s)\n", x$p_h0, x$in_tail))
  cat(sprintf("  H1 central [%.4g, %.4g]  obs at H1 pctile %.2f  consistent: %s%s\n",
              x$h1_interval[1], x$h1_interval[2], x$h1_percentile,
              x$h1_consistent, if (isTRUE(x$borderline)) "  [BORDERLINE]" else ""))
  cat(sprintf("  DECISION: %s%s\n", toupper(x$decision),
              if (isTRUE(x$borderline)) " (borderline -- label is fragile)" else ""))
  if (!is.null(x$grounding)) {
    cat(sprintf("  GROUNDING: %s%s\n", x$grounding,
                if (identical(x$grounding, "[unverified]"))
                  " (mechanism provenance not declared -- verdict not grounded)"
                else ""))
  }
  if (!isTRUE(x$posterior_adequacy$ok)) {
    cat(sprintf("  [!] posterior-adequacy WARNING: %s\n",
                x$posterior_adequacy$reason))
  }
  invisible(x)
}
