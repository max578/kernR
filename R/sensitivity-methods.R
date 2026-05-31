# sensitivity-methods.R -- S3 methods for hsic_sensitivity objects.
#
# print, plot, and as.data.frame methods for the objects returned by
# hsic_sensitivity(). Split out of sensitivity.R to keep the estimator and its
# display logic in separate files; the constructor and its documentation live
# in R/sensitivity.R.

#' @export
print.hsic_sensitivity <- function(x, digits = 3L, ...) {
  cat("\n  HSIC-Sensitivity Indices\n\n")
  cat("Parameters:  ", length(x$param_names), "\n")
  cat("Outputs:     ", length(x$output_names), "\n")
  cat("N:           ", x$n, "\n")
  if (!is.null(x$p_value)) {
    cat("Permutations:", x$n_permutations, "\n")
    cat("P-adjust:    ", x$p_adjust, "\n")
  } else {
    cat("Permutations: skipped (p_value = FALSE)\n")
  }
  cat("Total-order: ", if (isTRUE(x$total_order)) "yes" else "no",
      "\n\n", sep = "")

  ord <- x$rank
  if (!is.null(x$p_value_adjusted)) {
    min_p <- apply(x$p_value_adjusted[ord, , drop = FALSE], 1L, min)
    p_str <- formatC(min_p, digits = 4L, format = "f")
  } else {
    p_str <- rep("--", length(ord))
  }

  first_max <- formatC(x$total_index[ord], digits = digits, format = "f")

  if (isTRUE(x$total_order)) {
    total_max <- apply(x$index_total_order[ord, , drop = FALSE], 1L, max)
    interaction <- total_max - x$total_index[ord]
    tab <- data.frame(
      parameter   = x$param_names[ord],
      S_first_max = first_max,
      T_total_max = formatC(total_max, digits = digits, format = "f"),
      interaction = formatC(interaction, digits = digits, format = "f"),
      min_p_first = p_str,
      stringsAsFactors = FALSE
    )
    if (isTRUE(x$total_order_ci)) {
      lo <- apply(x$ci_total_order_lower[ord, , drop = FALSE], 1L, min)
      hi <- apply(x$ci_total_order_upper[ord, , drop = FALSE], 1L, max)
      tab$T_CI <- paste0(
        "[", formatC(lo, digits = digits, format = "f"), ", ",
        formatC(hi, digits = digits, format = "f"), "]"
      )
    }
    if (identical(x$total_order_test, "cond_perm")) {
      pT_adj <- x$p_value_total_order_adjusted
      tab$min_p_total <- formatC(
        apply(pT_adj[ord, , drop = FALSE], 1L, min),
        digits = 4L, format = "f"
      )
    }
    cat("Per-parameter ranking (descending S, max across outputs):\n")
    cat("  S = first-order index   T = total-order index   interaction = T - S\n")
    if (isTRUE(x$total_order_ci)) {
      cat("  T_CI = pair-bootstrap ",
          formatC(100 * x$ci_level, digits = 0, format = "f"),
          "% percentile CI on T (B = ", x$n_bootstrap,
          " resamples). NOT a significance test for T = 0; see ",
          "?hsic_sensitivity Details.\n", sep = "")
    }
    if (identical(x$total_order_test, "cond_perm")) {
      cat("  min_p_total = conditional-permutation p-value for ",
          "H_0: X_j _||_ Y | X_{~j} (", x$n_permutations,
          " permutations, ", x$p_adjust, "-adjusted across grid).\n",
          sep = "")
    }
    cat("\n")
  } else {
    tab <- data.frame(
      parameter   = x$param_names[ord],
      S_first_max = first_max,
      min_p_first = p_str,
      stringsAsFactors = FALSE
    )
    cat("Per-parameter ranking (descending S, max across outputs):\n")
    cat("  S = first-order index\n\n")
  }
  print(tab, row.names = FALSE)
  cat("\n")
  invisible(x)
}


#' Plot HSIC-Sensitivity Indices
#'
#' Bar plot of per-parameter HSIC-Sensitivity Index, ordered by
#' first-order index magnitude. When `total_order` was set on the fit,
#' `which = "total"` shows the total-order indices and `which = "both"`
#' shows side-by-side bars for `S` and `T`.
#'
#' @param x An `hsic_sensitivity` object.
#' @param which Character. `"first"` (default), `"total"`, or `"both"`.
#'   `"total"` / `"both"` require `x$total_order` to be `TRUE`.
#' @param alpha Numeric in `(0, 1)`. Significance level for colour
#'   coding (first-order only; ignored if the object carries no
#'   p-values or under `which = "total"`).
#' @param col_sig,col_nonsig,col_total Bar colours.
#' @param ... Additional arguments passed to [graphics::barplot()].
#'
#' @return Invisibly returns `x`. Side effect: produces a base R plot.
#' @export
plot.hsic_sensitivity <- function(x,
                                  which      = c("first", "total", "both"),
                                  alpha      = 0.05,
                                  col_sig    = "#0072B2",
                                  col_nonsig = "grey70",
                                  col_total  = "#D55E00",
                                  ...) {
  which <- match.arg(which)
  if (which %in% c("total", "both") && !isTRUE(x$total_order)) {
    stop("`which = \"",
         which,
         "\"` requires the fit was made with `total_order = TRUE`.",
         call. = FALSE)
  }
  ord <- x$rank
  first_max <- x$total_index[ord]
  param_lab <- x$param_names[ord]

  if (which == "first") {
    cols <- if (!is.null(x$p_value_adjusted)) {
      min_p <- apply(x$p_value_adjusted[ord, , drop = FALSE], 1L, min)
      ifelse(min_p <= alpha, col_sig, col_nonsig)
    } else {
      rep(col_sig, length(ord))
    }
    graphics::barplot(
      first_max, names.arg = param_lab,
      col = cols, border = NA,
      main = "HSIC-Sensitivity Index (first-order)",
      ylab = "S^HSIC", ylim = c(0, max(first_max) * 1.1),
      las = 2L, ...
    )
    if (!is.null(x$p_value_adjusted)) {
      graphics::legend(
        "topright",
        legend = c(paste0("significant (alpha=", alpha, ")"),
                   "not significant"),
        fill = c(col_sig, col_nonsig), border = NA, bty = "n"
      )
    }
  } else if (which == "total") {
    total_max <- apply(x$index_total_order[ord, , drop = FALSE], 1L, max)
    graphics::barplot(
      total_max, names.arg = param_lab,
      col = col_total, border = NA,
      main = "HSIC-Sensitivity Index (total-order)",
      ylab = "T^HSIC", ylim = c(0, max(total_max) * 1.1),
      las = 2L, ...
    )
  } else {  # both
    total_max <- apply(x$index_total_order[ord, , drop = FALSE], 1L, max)
    mat <- rbind(first = first_max, total = total_max)
    colnames(mat) <- param_lab
    ymax <- max(mat) * 1.15
    graphics::barplot(
      mat, beside = TRUE,
      col = c(col_sig, col_total), border = NA,
      main = "HSIC-Sensitivity Indices (first-order vs total-order)",
      ylab = "Index value", ylim = c(0, ymax),
      las = 2L, ...
    )
    graphics::legend(
      "topright",
      legend = c("S (first-order)", "T (total-order)"),
      fill = c(col_sig, col_total), border = NA, bty = "n"
    )
  }
  invisible(x)
}


#' @export
as.data.frame.hsic_sensitivity <- function(x, row.names = NULL,
                                           optional = FALSE, ...) {
  p <- length(x$param_names)
  q <- length(x$output_names)
  out <- data.frame(
    parameter        = rep(x$param_names, times = q),
    output           = rep(x$output_names, each  = p),
    index            = as.vector(x$index),
    statistic        = as.vector(x$statistic),
    stringsAsFactors = FALSE
  )
  if (!is.null(x$p_value)) {
    out$p_value          <- as.vector(x$p_value)
    out$p_value_adjusted <- as.vector(x$p_value_adjusted)
  }
  if (isTRUE(x$total_order)) {
    out$index_total_order     <- as.vector(x$index_total_order)
    out$statistic_total_order <- as.vector(x$statistic_total_order)
    if (isTRUE(x$total_order_ci)) {
      out$ci_total_order_lower <- as.vector(x$ci_total_order_lower)
      out$ci_total_order_upper <- as.vector(x$ci_total_order_upper)
    }
    if (identical(x$total_order_test, "cond_perm")) {
      out$p_value_total_order          <- as.vector(x$p_value_total_order)
      out$p_value_total_order_adjusted <-
        as.vector(x$p_value_total_order_adjusted)
    }
  }
  out
}
