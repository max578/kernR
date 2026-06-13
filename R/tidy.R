# tidy.R -- broom-style tidiers for kernR test-result objects.
#
# kernR's verdict objects (kernel_test_result and its subclasses, taci_result)
# are S3 lists with hand-named fields. Downstream scripts that wanted a flat
# one-row-per-term summary had to reach into those fields by name, which is
# fragile -- the field is `p_value`, not the broom-canonical `p.value`, so two
# consuming scripts guessed wrong. These methods register against the generic
# `tidy()` (re-exported by broom) so `broom::tidy()` and `generics::tidy()`
# both return a stable, documented data.frame for every kernR verdict.

#' @importFrom generics tidy
#' @export
generics::tidy

#' Tidy a kernR Test Result Into a One-Row-Per-Term Data Frame
#'
#' Turns a kernR verdict object into a flat, `broom`-style summary
#' `data.frame`: one row per tested term, with stable, documented columns. The
#' method is registered against the `tidy()` generic (re-exported by `broom`),
#' so `broom::tidy(result)` and `generics::tidy(result)` both dispatch here.
#'
#' This is the supported accessor for downstream code. Reaching into the result
#' list by field name is fragile: kernR's native field is `p_value` (with an
#' underscore), whereas the `broom` convention -- followed here -- is `p.value`
#' (with a dot). Use `tidy()` and read `p.value` rather than guessing at the
#' raw field name.
#'
#' @param x A `kernel_test_result` (the class returned by [hsic_test()],
#'   [mmd_test()], [bd_hsic_test()], and friends) or a subclass such as
#'   `mmd_ppc`.
#' @param ... Currently unused; present for generic compatibility.
#'
#' @returns A one-row `data.frame` with columns:
#'   \describe{
#'     \item{term}{Character. The quantity tested -- the test `method`
#'       (for example `"bd-HSIC"`, `"MMD"`).}
#'     \item{statistic}{Numeric. The observed test statistic.}
#'     \item{p.value}{Numeric. The permutation p-value. Note the
#'       `broom`-canonical dot, distinct from the result's native `p_value`
#'       field.}
#'     \item{n}{Integer. The sample size the statistic was computed on.}
#'     \item{n_permutations}{Integer. Number of permutations in the null.}
#'     \item{ess}{Numeric. Effective sample size of the importance weights, or
#'       `NA` for unweighted tests.}
#'   }
#'   Subclasses contribute extra columns where they carry extra fields: an
#'   `mmd_ppc` result adds `surprise_bits` and `reject`.
#'
#' @seealso [hsic_test()], [mmd_test()], [bd_hsic_test()], [mmd_ppc()]
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(200L), 100L, 2L)
#' y <- x[, 1L] + rnorm(100L)
#' res <- hsic_test(x, y, n_permutations = 199L, seed = 1L)
#' tidy(res)
#' @export
tidy.kernel_test_result <- function(x, ...) {
  out <- data.frame(
    term = as.character(x$method),
    statistic = as.numeric(x$statistic),
    p.value = as.numeric(x$p_value),
    n = as.integer(x$n),
    n_permutations = as.integer(x$n_permutations),
    ess = if (is.null(x$ess)) NA_real_ else as.numeric(x$ess),
    stringsAsFactors = FALSE
  )

  # Subclass extras --------------------------------------------------------
  # An mmd_ppc result carries the posterior-predictive surprise and the
  # reject/accept verdict; surface them so a PPC tidies to a self-contained row.
  if (inherits(x, "mmd_ppc")) {
    out$surprise_bits <- as.numeric(x$surprise_bits)
    out$reject <- as.logical(x$reject)
  }

  out
}

#' Tidy a TACI Mechanism-Consistency Result
#'
#' Turns a [taci_test()] verdict into a flat, `broom`-style one-row summary.
#' Registered against the `tidy()` generic (re-exported by `broom`), so
#' `broom::tidy()` and `generics::tidy()` both dispatch here.
#'
#' The `p.value` column carries the broom-canonical dot and reports the H0 tail
#' p-value (`p_h0`). The `grounding` column threads through the Independent
#' Oracle Principle label, so a downstream summary cannot silently drop the
#' fact that a verdict built on an un-declared mechanism is `"[unverified]"`.
#'
#' @param x A `taci_result` object from [taci_test()].
#' @param ... Currently unused; present for generic compatibility.
#'
#' @returns A one-row `data.frame` with columns:
#'   \describe{
#'     \item{term}{Character. Always `"taci"`.}
#'     \item{statistic}{Numeric. The observed bd-HSIC statistic.}
#'     \item{p.value}{Numeric. The H0 tail p-value (`p_h0`).}
#'     \item{decision}{Character. The three-way decision enum
#'       (`"no_effect"`, `"mechanism_consistent_effect"`,
#'       `"mechanism_inconsistent_effect"`).}
#'     \item{grounding}{Character. `"grounded"` when the mechanism's provenance
#'       was declared, else `"[unverified]"` (Independent Oracle Principle).}
#'     \item{h1_percentile}{Numeric. Where the observed statistic sits in the
#'       model-implied H1 band, in `[0, 1]`.}
#'     \item{borderline}{Logical. Whether the consistency label is fragile.}
#'     \item{n}{Integer. The sample size.}
#'   }
#'
#' @seealso [taci_test()]
#' @examples
#' set.seed(1)
#' n <- 60L
#' nrate <- runif(n, 0, 200)
#' yield <- 1.1 + 4.2 * (1 - exp(-0.018 * nrate)) + rnorm(n, 0, 0.25)
#' post <- cbind(ymax = rnorm(150L, 4.2, 0.30),
#'               rate = rnorm(150L, 0.018, 0.004),
#'               y0   = rnorm(150L, 1.1, 0.15))
#' mitscherlich <- function(theta, X, t) {
#'   theta[3L] + theta[1L] * (1 - exp(-theta[2L] * t))
#' }
#' res <- taci_test(post, mitscherlich, X = matrix(1, n, 1L),
#'                  treatment = nrate, outcome = yield,
#'                  n_perm = 80L, seed = 1L)
#' tidy(res)
#' @export
tidy.taci_result <- function(x, ...) {
  data.frame(
    term = "taci",
    statistic = as.numeric(x$observed_statistic),
    p.value = as.numeric(x$p_h0),
    decision = as.character(x$decision),
    grounding = as.character(x$grounding),
    h1_percentile = as.numeric(x$h1_percentile),
    borderline = as.logical(x$borderline),
    n = as.integer(x$n),
    stringsAsFactors = FALSE
  )
}
