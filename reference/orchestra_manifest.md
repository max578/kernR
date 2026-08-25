# The orchestra ensemble-manifest contract (kernR-side implementation)

A versioned, hashed, provenance-complete S7 result object. It is not
merely property-compatible with the federation's reference
`orchestra_manifest` – it is that class: declared with `package = NULL`,
so its S7 identity is the bare name the reference contract dispatches on
rather than a `kernR`-namespaced look-alike a consumer would refuse. A
`taci_result` or a `kernel_test_result` is adapted into one through
[`as_orchestra_manifest()`](https://max578.github.io/kernR/reference/as_orchestra_manifest.md):
the verdict rides in the typed `summary`
(`inferential_target = "treatment_effects"`), the test provenance rides
in `metadata`, and the payload is integrity-hashed so a tampered verdict
is detected downstream.

## Value

An S7 object of class `orchestra_manifest`.

## See also

[`as_orchestra_manifest()`](https://max578.github.io/kernR/reference/as_orchestra_manifest.md),
[`verify_manifest()`](https://max578.github.io/kernR/reference/verify_manifest.md)

Other orchestra manifest:
[`as_orchestra_manifest()`](https://max578.github.io/kernR/reference/as_orchestra_manifest.md),
[`verify_manifest()`](https://max578.github.io/kernR/reference/verify_manifest.md)
