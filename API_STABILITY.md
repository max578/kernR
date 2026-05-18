# API stability policy

`kernR` follows a published two-phase policy. The phase boundary is
the 1.0.0 release.

## Pre-1.0 (current)

Versions `0.x.y` follow **additive-by-intent** evolution. Every minor
release (`0.x → 0.(x+1)`) adds new exports without breaking existing
signatures; patch releases fix bugs and refresh documentation.

`0.1.0` is the first tagged public release. Pre-release development
proceeded under `0.0.0.9001 → 0.0.0.9015`; `NEWS.md` preserves the
per-feature changelog of that cycle as historical context.

The policy is "additive-by-intent" rather than "frozen": pre-1.0
reserves the right to break an existing signature when a design flaw
surfaces, but every such change must be:

1. Listed under a `## Breaking changes` heading in `NEWS.md`, first.
2. Justified in the release notes.
3. Where feasible, accompanied by a temporary back-compat shim.

The 0.0.0.9013 retraction of the un-calibrated Pick-Freeze
`total_order_p_value` mode is the precedent: the argument was made
defunct with an error pointing at both the safe CI option
(`total_order_ci`) and the replacement significance test
(`total_order_test = "cond_perm"` in 0.0.0.9014).

If you depend on `kernR` pre-1.0, pin the version in `renv.lock` or
`DESCRIPTION` (`Imports: kernR (>= 0.1.0)`).

## 1.0 and after

From `1.0.0`, `kernR` adopts **strict frozen-API additive evolution** —
the same policy as `glmnet`, `mgcv`, and other long-lived statistical
inference packages:

- Signatures of exported functions never change in a
  backwards-incompatible way across major versions.
- New capability arrives via new entry points (new exports), never
  via changes to existing ones.
- If an existing function genuinely needs to be retired, it is
  marked with `lifecycle::deprecate_warn()`, kept for ≥ 2 minor
  versions, then promoted to `lifecycle::deprecate_stop()` for
  ≥ 1 more minor version, then removed in the next major release.
  Successors are signposted in the deprecation message.

This policy is chosen because `kernR` provides hypothesis-testing
machinery that downstream consumers (the PESTO ⇄ kernR ⇄ proxymix
chain; agricultural pipeline code; methods-paper supplementary
material) wire into long-lived calibration and verdict workflows.
Silent breakage across versions would defeat that contract.

## Cross-package contracts

`kernR` imports `PESTO (>= 0.3.0)` and dispatches on
`PESTO::pesto_ensemble_manifest` via S3 methods registered at
`.onLoad()`. The contract is owned jointly with PESTO; any change to
the S7 class shape on PESTO's side requires a coordinated kernR
release (and a corresponding `Imports:` lower bound bump).

`proxymix (>= 0.3.0)` in `Suggests` provides the
`density_ratio = "proxymix"` backend. The interface is the
`proxymix::fit_proxymix()` / `proxymix::dgmm()` pair; backwards-
incompatible changes there will be detected at run-time via the
`requireNamespace()` guard.

## Versioning and tags

`kernR` uses [Semantic Versioning 2.0.0](https://semver.org/). Every
release is git-tagged `vX.Y.Z` on `main`; the tag is annotated and
carries the `NEWS.md` entry as its message.

## Reporting an unintended break

If you discover that a `kernR` release has silently broken your
pipeline, open an issue at
<https://github.com/AAGI-AUS/kernR/issues> with the version you
upgraded from and to plus a small reproducible example. Unintended
breaks at any stage are treated as bugs.
