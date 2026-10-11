#' Deprecated: use mechanism_consistency_test()
#'
#' `taci_test()` is the former name of [mechanism_consistency_test()]. It
#' forwards all arguments and gives a deprecation warning.
#'
#' @param ... Arguments passed to [mechanism_consistency_test()].
#' @return A `mechanism_consistency_result`.
#' @keywords internal
#' @export
taci_test <- function(...) {
  .Deprecated("mechanism_consistency_test", package = "kernR")
  mechanism_consistency_test(...)
}
