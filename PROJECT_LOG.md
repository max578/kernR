# PROJECT_LOG.md — kernR

## Compressed Project Context

**Objective**: Build an R package (`kernR`) implementing kernel-based causal distributional tests from Hu/Sejdinovic/Evans (JMLR 2024) and Fawkes/Hu/Evans/Sejdinovic (TMLR 2024), with novel hierarchical extensions.
**Core methods**: bd-HSIC (backdoor-adjusted HSIC for causal association), DR-DATE/DR-DETT (doubly robust distributional treatment effect tests).
**Technical approach**: Rcpp/RcppArmadillo kernel engine, R-level density ratio and propensity estimation, permutation inference. S3 classes, pipe-friendly API, formula interface.
**Key assumptions**: Backdoor criterion holds (no unmeasured confounding); characteristic kernels (RBF, Matern); valid permutation via clustering/binning.
**Novel contribution**: Hierarchical/nested data support via within/between decomposition of test statistics.
**Status**: Phases 0-5 complete. 69 tests passing. 4 vignettes built. R CMD check: 0 errors.
**Open issues**: System FLIBS broken (workaround documented); placeholder ORCIDs; no Git repo yet; Phases 6-7 (performance, polish, release) pending.

---

## Log Entries

### 2026-04-09 — v0.0.0.9000 — Initial Development

**Summary**: Created complete package skeleton through Phase 5 (hierarchical extensions).

**Key Decisions & Rationale**:
- Chose S3 over S4/R5 for simplicity and pipe-friendliness
- Dropped `osqp` from Suggests (not available on system; add when KMM is implemented)
- Removed OpenMP from Makevars (link errors on macOS; defer to Phase 6)
- Density ratio: logistic NCE via glm/ranger instead of neural NCE (lighter dependency)
- Hierarchical test uses within/between decomposition with configurable weighting (equal, ICC, within-only)
- FLIBS workaround: local builds strip $(FLIBS) from Makevars; restore for CRAN compliance

**Technical Updates**:
- 5 Rcpp source files: kernels.cpp, hsic.cpp, mmd.cpp, permutation.cpp, density_ratio.cpp
- 16 R source files covering full API
- 69 tests across 7 test files, all passing
- 4 vignettes (quickstart, bdhsic, drtest, hierarchical) — all building
- R CMD check: 0 errors, 1 warning (FLIBS portability), 2 notes (dev version, ORCID)

**Next Steps**:
- Phase 6: Nystrom approximation, RFF, torch backend, adaptive permutation, theory + performance vignettes, pkgdown
- Phase 7: Real data examples (Lalonde, IHDP), comparison benchmarks, JOSS paper, CRAN submission
- Fix ORCIDs (get Max's real ORCID, verify Dino's)
- Set up Git repository and CI/CD
