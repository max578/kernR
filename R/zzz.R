.onLoad <- function(libname, pkgname) {
  # Register S3 methods for PESTO's package-qualified S7 class name.
  # S7 sets `class(x)` = c("PESTO::pesto_ensemble_manifest",
  # "S7_object"); the bare-name S3 method file
  # `mmd_ppc.pesto_ensemble_manifest` is not reachable through
  # standard UseMethod() dispatch because R cannot parse `::` in an
  # S3-method identifier. registerS3method() wires it up at load time.
  ns <- asNamespace(pkgname)
  registerS3method(
    "mmd_ppc", "PESTO::pesto_ensemble_manifest",
    get("mmd_ppc.pesto_ensemble_manifest", envir = ns),
    envir = ns
  )
  invisible()
}
