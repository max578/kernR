# abstention.R -- the typed "reliability gate failed" marker (K4 residual,
# closing K8).
#
# kernR's three reliability gates -- the ESS floor (`bd_hsic_test()`,
# `dr_date_test()`, `dr_dett_test()`, `aggregate_downscale()`), the proxymix
# density-ratio fit-quality gate (`bd_hsic_test()`), and the TACI
# posterior-adequacy guard (`taci_test()`) -- already computed a reliability
# flag (`ess_warning`, `density_ratio_warning`, `posterior_adequacy$ok`) and
# an ordinary `warning()`, but attached nothing a cross-member gate could
# recognise without knowing kernR's field names. The leader-node
# `is_orchestra_decline()` predicate (`ORCHESTRA_dev/integration/
# refusal_contract.R`) names "kernR ESS-floor" as an abstention it
# recognises via the `_(refusal|abstention)$` class-name convention (the same
# route `gpfield_abstention` uses); until now nothing in kernR carried that
# class, so the predicate returned FALSE for every kernR reliability
# failure.
#
# Stamping does not change WHAT is refused: the underlying flag and
# `warning()` are untouched, and the object still carries its full result
# (a degraded verdict, not a hard refusal, remains inspectable by hand).
# This only prepends a class name so a cross-member gate can ask "did this
# member decline to stand behind its result?" once, generically.

#' Stamp the typed `kernR_abstention` marker onto a result
#'
#' Prepends `"kernR_abstention"` to `x`'s class, idempotently, so
#' `inherits(x, "kernR_abstention")` -- and the federation's
#' `is_orchestra_decline()`, which matches on the `_(refusal|abstention)$`
#' suffix -- recognise a kernR result whose reliability gate failed. Does not
#' alter any field on `x`; the raw reliability flags this function is called
#' from (`ess_warning`, `density_ratio_warning`,
#' `posterior_adequacy$ok`) remain exactly as computed.
#'
#' @param x A `kernel_test_result` or `taci_result` whose reliability gate
#'   has already failed.
#'
#' @returns `x` with `"kernR_abstention"` prepended to its class.
#' @noRd
#' @keywords internal
.mark_kernR_abstention <- function(x) {
  if (!any(class(x) == "kernR_abstention")) {
    class(x) <- c("kernR_abstention", class(x))
  }
  x
}

#' Did a kernR result decline to stand behind its own verdict?
#'
#' `TRUE` when a reliability gate failed (the ESS floor, the proxymix
#' density-ratio fit-quality gate, or the TACI posterior-adequacy guard) and
#' the result was stamped with the typed `kernR_abstention` marker. This is
#' kernR's local mirror of the federation's `is_orchestra_decline()`
#' predicate (`ORCHESTRA_dev/integration/refusal_contract.R`), scoped to the
#' member's own result classes so calling code inside kernR does not need to
#' depend on the composition-layer contract file.
#'
#' @param x Any object.
#'
#' @returns A length-1 logical.
#'
#' @examples
#' set.seed(42L)
#' n <- 60L
#' z <- matrix(rnorm(n * 2L), n, 2L)
#' x <- z[, 1L] + rnorm(n, sd = 0.5)
#' y <- 0.5 * x + rnorm(n, sd = 0.5)
#' res <- suppressWarnings(bd_hsic_test(x, y, z, min_ess_fraction = 0.999))
#' is_kernR_abstention(res)  # TRUE -- the ESS floor tripped
#'
#' @export
is_kernR_abstention <- function(x) {
  inherits(x, "kernR_abstention")
}
