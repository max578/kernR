# Changelog

## kernR 0.3.1

### Testing

- New end-to-end analytical-correctness test (`test-end-to-end.R`). One
  confounded data-generating process is run through the whole pipeline
  and every analytical corner is checked for a meaningful, correct
  verdict: marginal
  [`hsic_test()`](https://max578.github.io/kernR/reference/hsic_test.md)
  power and independence;
  [`bd_hsic_test()`](https://max578.github.io/kernR/reference/bd_hsic_test.md)
  removing a purely confounder-induced association while detecting a
  genuine causal one; propensity recovery;
  [`dr_date_test()`](https://max578.github.io/kernR/reference/dr_date_test.md)
  power, Type I control, and double robustness (AIPW and IPW-only
  fallback agree);
  [`dr_dett_test()`](https://max578.github.io/kernR/reference/dr_dett_test.md);
  two-sample
  [`mmd_test()`](https://max578.github.io/kernR/reference/mmd_test.md);
  Nystrom and RFF agreement with exact HSIC; hierarchical within-cluster
  Type I control; full-pipeline seed reproducibility; and
  permutation-null calibration.
- New coverage for the public density-ratio API and weight diagnostics
  (`test-density-ratio-api.R`):
  [`fit_density_ratio()`](https://max578.github.io/kernR/reference/fit_density_ratio.md)
  /
  [`predict_density_ratio()`](https://max578.github.io/kernR/reference/predict_density_ratio.md)
  round-trip and reproducibility,
  [`plot_weights()`](https://max578.github.io/kernR/reference/plot_weights.md),
  and the new
  [`print.cme_fit()`](https://max578.github.io/kernR/reference/print.cme_fit.md)
  method.

### Minor improvements and fixes

- New [`print()`](https://rdrr.io/r/base/print.html) method for
  [`fit_cme()`](https://max578.github.io/kernR/reference/fit_cme.md)
  objects
  ([`print.cme_fit()`](https://max578.github.io/kernR/reference/print.cme_fit.md)).
- `@family` tags added across the exported API so the documentation
  cross-links related functions;
  [`dr_date_scenario()`](https://max578.github.io/kernR/reference/dr_date_scenario.md)’s
  example now runs (`\donttest` rather than `\dontrun`).
- `R/sensitivity.R` split: the
  [`hsic_sensitivity()`](https://max578.github.io/kernR/reference/hsic_sensitivity.md)
  S3 methods now live in `R/sensitivity-methods.R`.
- Dropped the `Remotes: github::max578/PESTO` line now that PESTO (\>=
  0.4.1) is served from the max578 r-universe; the CI
  `extra-repositories` entry resolves it.

## kernR 0.3.0

### Correctness

- **DR-DATE and DR-DETT are now genuinely doubly robust.**
  [`dr_date_test()`](https://max578.github.io/kernR/reference/dr_date_test.md)
  previously fitted a conditional mean embedding outcome model and then
  discarded it, returning an inverse-probability-weighted statistic
  regardless of `outcome_model` – so `outcome_model = "krr"` and
  `"zero"` gave identical results despite the documented
  double-robustness claim. The statistic now forms the augmented (AIPW)
  counterfactual mean embeddings, consistent if *either* the propensity
  or the outcome model is correctly specified (Fawkes, Hu, Evans &
  Sejdinovic, 2024).
- **DR-DETT control reweighting corrected.** The effect-on-the-treated
  control counterfactual was reweighting controls by the inverted
  treatment odds `(1 - e) / e`; it now uses the correct treatment odds
  `e / (1 - e)` with an augmented outcome-model correction on the
  control arm.
- **`seed=` now makes permutation tests reproducible.** The C++
  permutation routines drew from Armadillo’s internal RNG, which ignores
  [`set.seed()`](https://rdrr.io/r/base/Random.html); they now draw
  through R’s RNG, so a fixed `seed` reproduces the null distribution
  and p-value of
  [`hsic_test()`](https://max578.github.io/kernR/reference/hsic_test.md),
  [`mmd_test()`](https://max578.github.io/kernR/reference/mmd_test.md),
  [`bd_hsic_test()`](https://max578.github.io/kernR/reference/bd_hsic_test.md),
  [`dr_date_test()`](https://max578.github.io/kernR/reference/dr_date_test.md),
  [`dr_dett_test()`](https://max578.github.io/kernR/reference/dr_dett_test.md)
  and the permutation-based sensitivity paths exactly.
- **[`hierarchical_test()`](https://max578.github.io/kernR/reference/hierarchical_test.md)
  no longer silently swallows within-cluster failures.** Failed
  within-cluster sub-tests are counted and surfaced via a warning (or an
  error if every eligible cluster fails). Clusters too small for the
  chosen sub-test (DR sub-tests need at least 30 observations) are now
  skipped with a clear message and *excluded* from the within-cluster
  average, rather than silently contributing zero.
- Integer-overflow hardening in the C++ HSIC normaliser
  (`(double)(n * n)` to `(double)n * n`).

### Features

- [`dr_date_test()`](https://max578.github.io/kernR/reference/dr_date_test.md)
  and
  [`dr_dett_test()`](https://max578.github.io/kernR/reference/dr_dett_test.md)
  gain `cross_fit` (default `TRUE`) and `n_folds` arguments: both
  nuisances are cross-fitted and evaluated out-of-fold, as the doubly
  robust theory requires under flexible nuisance estimators
  (Chernozhukov et al., 2018).
- [`dr_date_test()`](https://max578.github.io/kernR/reference/dr_date_test.md)
  and
  [`dr_dett_test()`](https://max578.github.io/kernR/reference/dr_dett_test.md)
  gain `min_ess_fraction` and now report the per-arm effective sample
  size (`ess`, `ess_warning`), warning when the inverse-probability
  weights collapse.
- [`estimate_propensity()`](https://max578.github.io/kernR/reference/estimate_propensity.md)
  gains a `seed` argument for reproducible cross-fitting folds.

### Documentation and packaging

- Dropped unused `Suggests` (`future`, `future.apply`, `ggplot2`,
  `viridisLite`); raised the `PESTO` floor to the tested `>= 0.4.1`.
- Removed an unpublishable talk citation from
  [`aggregate_downscale()`](https://max578.github.io/kernR/reference/aggregate_downscale.md).
- Corrected a mislabelled “two-sided” comment on the (correct, one-sided
  upper-tail) HSIC permutation p-value.

## kernR 0.2.0

First publish of the local development line to AAGI-AUS. Lands the
0.0.0.9001 → 0.0.0.9015 cycle on top of the existing AAGI v0.1.x
lineage; subsequent 0.0.0.x history below is retained as the per-feature
changelog. The version-number jump (origin v0.1.1 → v0.2.0) signals the
substantial new public surface arriving in this release.

### Public surface (31 exports)

- **Kernel primitives**:
  [`kernel_spec()`](https://max578.github.io/kernR/reference/kernel_spec.md),
  [`kernel_matrix()`](https://max578.github.io/kernR/reference/kernel_matrix.md),
  [`select_bandwidth()`](https://max578.github.io/kernR/reference/select_bandwidth.md).
- **Independence and two-sample tests**:
  [`hsic_test()`](https://max578.github.io/kernR/reference/hsic_test.md),
  [`mmd_test()`](https://max578.github.io/kernR/reference/mmd_test.md),
  [`hsic_test_nystrom()`](https://max578.github.io/kernR/reference/hsic_test_nystrom.md).
- **Causal association**:
  [`bd_hsic_test()`](https://max578.github.io/kernR/reference/bd_hsic_test.md)
  (backdoor-adjusted HSIC; Hu, Sejdinovic & Evans, JMLR 2024).
- **Distributional treatment effects**:
  [`dr_date_test()`](https://max578.github.io/kernR/reference/dr_date_test.md),
  [`dr_dett_test()`](https://max578.github.io/kernR/reference/dr_dett_test.md),
  [`dr_date_scenario()`](https://max578.github.io/kernR/reference/dr_date_scenario.md),
  [`kernel_causal_test()`](https://max578.github.io/kernR/reference/kernel_causal_test.md)
  (Fawkes, Hu, Evans & Sejdinovic, TMLR 2024).
- **Hierarchical / clustered designs**:
  [`hierarchical_test()`](https://max578.github.io/kernR/reference/hierarchical_test.md)
  with within-cluster permutation.
- **Sensitivity and identifiability**:
  [`hsic_identifiability()`](https://max578.github.io/kernR/reference/hsic_identifiability.md),
  [`hsic_sensitivity()`](https://max578.github.io/kernR/reference/hsic_sensitivity.md)
  (Da Veiga 2015; conditional-permutation null for total-order
  significance).
- **Density-ratio and propensity**:
  [`fit_density_ratio()`](https://max578.github.io/kernR/reference/fit_density_ratio.md),
  [`predict_density_ratio()`](https://max578.github.io/kernR/reference/predict_density_ratio.md),
  [`estimate_density_ratio()`](https://max578.github.io/kernR/reference/estimate_density_ratio.md)
  (logistic / ranger / xgboost / proxymix backends);
  [`estimate_propensity()`](https://max578.github.io/kernR/reference/estimate_propensity.md),
  [`assess_overlap()`](https://max578.github.io/kernR/reference/assess_overlap.md),
  [`plot_weights()`](https://max578.github.io/kernR/reference/plot_weights.md),
  [`effective_sample_size()`](https://max578.github.io/kernR/reference/effective_sample_size.md).
- **Low-rank acceleration**:
  [`nystrom_factor()`](https://max578.github.io/kernR/reference/nystrom_factor.md)
  (Williams & Seeger 2001),
  [`rff_features()`](https://max578.github.io/kernR/reference/rff_features.md)
  (Rahimi & Recht 2007).
- **Kernel downscaling and distribution regression**:
  [`kernel_downscale()`](https://max578.github.io/kernR/reference/kernel_downscale.md)
  (Park, Muandet, Fukumizu & Sejdinovic 2013),
  [`fit_cme()`](https://max578.github.io/kernR/reference/fit_cme.md),
  [`dist_regression()`](https://max578.github.io/kernR/reference/dist_regression.md)
  (Szabó, Sriperumbudur, Póczos & Gretton 2016),
  [`aggregate_downscale()`](https://max578.github.io/kernR/reference/aggregate_downscale.md),
  [`posterior_sample_aggregate()`](https://max578.github.io/kernR/reference/posterior_sample_aggregate.md).
- **Posterior-predictive check + PESTO contract**:
  [`mmd_ppc()`](https://max578.github.io/kernR/reference/mmd_ppc.md)
  (consumes
  [`PESTO::pesto_ensemble_manifest`](https://rdrr.io/pkg/PESTO/man/pesto_ensemble_manifest.html)
  via S3 dispatch).
- **Design**:
  [`lhs_design()`](https://max578.github.io/kernR/reference/lhs_design.md)
  Latin-hypercube helper.

### Cross-package contracts

- Imports `PESTO (>= 0.3.0)` for the `pesto_ensemble_manifest` S7 class.
- Optional `proxymix (>= 0.3.0)` in `Suggests` as a density-ratio
  backend
  ([`requireNamespace()`](https://rdrr.io/r/base/ns-load.html)-guarded).

### Documentation

- 12 vignettes covering: quick start, bd-HSIC tutorial, DR-DATE/ DETT
  tutorial, hierarchical data, HSIC identifiability, HSIC sensitivity,
  MMD posterior-predictive check, DR-DATE scenario via PESTO manifest,
  hierarchical bd-HSIC, Nyström acceleration, proxymix binding, kernel
  downscaling.

### R CMD check posture

`R CMD check --as-cran` on the AAGI / CI Linux environment is expected
to report 0 errors / 0 warnings / a small number of environmental NOTEs
(new-submission boilerplate, HTML Tidy version on macOS). Two
Apple-clang-21 toolchain WARNINGs surface only on the maintainer’s local
machine (R’s own `R_ext/Boolean.h` and a personal `~/.R/Makevars`); both
are absent on CRAN’s build farm and on standard GitHub Actions runners.

### Historical development log

All historical entries below are retained for full traceability
(0.0.0.9001 sole-authorship consolidation → 0.0.0.9015 third downscaling
method).

### kernR 0.0.0.9015

#### Third downscaling method: `aggregate_downscale()`

Closes the orchestra-completion gap surfaced after the 2026-05-16
deferral closeout: a downscaling method for the *no-paired-training-
data* regime, where only the aggregate observation and a known
aggregator are available. Companion to the proxymix Tier-2 stub
`from_aggregate_likelihood()`, hosted on the kernR (consumption) side so
proxymix’s CRAN pre-submission stays untouched.

- New export `aggregate_downscale(y, aggregator, latent_prior, ...)`.
  Dispatches on the aggregator’s class:
  - **Linear matrix `A`** — closed-form per-component Kalman update plus
    mixture-weight reweighting by per-component evidence
    `N(y | A mu_k, A Sigma_k A^T + sigma_y^2 I)`.
  - **Function `T(x)`** — importance sampling within each prior
    component (`n_samples_per_component` draws; ESS-floor
    `min_ess_fraction` reliability gate with explicit
    [`warning()`](https://rdrr.io/r/base/warning.html) on collapse).
- New export
  [`posterior_sample_aggregate()`](https://max578.github.io/kernR/reference/posterior_sample_aggregate.md)
  for drawing from the posterior mixture (downstream uncertainty
  propagation).
- Accepts the latent prior either as a list
  `(means, covariances, weights)` or as any object exposing those slots
  — including
  [`proxymix::fit_proxymix()`](https://rdrr.io/pkg/proxymix/man/fit_proxymix.html)
  results (gated via
  [`methods::slot()`](https://rdrr.io/r/methods/slot.html), no hard
  dependency added).
- Vignette `kernR-downscaling.Rmd` reframed from “Two flavours” to
  “Three flavours”; full comparison table.
- 8 new test blocks at `tests/testthat/test-aggregate-downscale.R`:
  single-component closed-form recovery (Kalman update exact);
  two-component mixture reweighting toward the likely cluster;
  non-linear IS path runs with valid moments and p.s.d. covariance;
  reproducibility under seed; ESS-floor warning fires under collapse;
  [`posterior_sample_aggregate()`](https://max578.github.io/kernR/reference/posterior_sample_aggregate.md)
  recovers the posterior mean; input validation; proxymix `gmm_fit`
  prior accepted via slot extraction (skipped when proxymix absent).

This is the **third** kernR downscaling method, structurally distinct
from
[`kernel_downscale()`](https://max578.github.io/kernR/reference/kernel_downscale.md)
(CME, paired training data) and
[`dist_regression()`](https://max578.github.io/kernR/reference/dist_regression.md)
(bag-of-points, distribution-level regression). Each occupies a
different cell of the (training-data-shape × aggregator-knowledge)
matrix.

### kernR 0.0.0.9014

#### Deferral closeout (post 2026-05-16 critical review)

Closes the three kernR-side tickets named as deferred in
`DRAINSTORMING/reports/orchestra_critical_review_response_2026-05-16.md`:
density-ratio fit/predict refactor, backend diagnostics, and properly
null-calibrated total-order significance test.

##### Density-ratio fit/predict refactor (closes P0 [\#2](https://github.com/max578/kernR/issues/2))

- New exports: `fit_density_ratio(x, z, method, ...)` and
  `predict_density_ratio(object, new_x, new_z, type = c("log_ratio", "weight", "ratio"))`.
  All four backends (`logistic`, `ranger`, `xgboost`, `proxymix`) return
  a `density_ratio_fit` object that can be applied to held-out rows.
- [`bd_hsic_test()`](https://max578.github.io/kernR/reference/bd_hsic_test.md)
  now **honours the documented train/test split**: fits density-ratio on
  the training half, predicts on the held-out test half. The runtime
  warning about the sample-split leak from 0.0.0.9013 is therefore
  retired. Closes P0 [\#2](https://github.com/max578/kernR/issues/2) of
  the critical review.
- Internal ratio computation is now in **log-space** end-to-end
  (`log_ratio = log(p_joint) - log(p_marg)` for proxymix;
  `log_ratio = log(p) - log(1-p) + log(n_noise)` for classifier
  backends). Numerically stable on extreme tails.
- `bd_hsic_test()$density_ratio_fit` carries the fitted model so callers
  can inspect backend diagnostics.
- [`estimate_density_ratio()`](https://max578.github.io/kernR/reference/estimate_density_ratio.md)
  retained as a thin backwards-compatible wrapper that fits and predicts
  on the same data; new code should prefer the explicit fit/predict
  pair.

##### Backend diagnostics for proxymix (closes P1 [\#3](https://github.com/max578/kernR/issues/3))

- `fit_density_ratio(method = "proxymix")` now surfaces per-GMM
  convergence diagnostics on `fit$diagnostics`: `joint_converged`,
  `marg_converged`, `joint_loglik`, `marg_loglik`, `joint_bic`,
  `marg_bic`, `joint_aic`, `marg_aic`, `joint_iterations`,
  `marg_iterations`, `n_components`.
- `print.density_ratio_fit()` summarises components + BIC + per-GMM
  convergence on screen.

##### Conditional-permutation total-order significance (closes P0 [\#1](https://github.com/max578/kernR/issues/1) remaining work)

- New arg `total_order_test = c("none", "cond_perm")`. When
  `"cond_perm"`: cluster the design points by `X_{~j}` similarity
  (k-means, `n_clusters_cp` bins; `"auto"` chooses
  `min(floor(n / 5), 20)`); within each cluster permute `Y`; recompute
  `T_j` on each of `n_permutations` permuted designs; report
  `p = (1 + #{T_perm ≥ T_obs}) / (1 + n_permutations)`.
- This is a **properly null-calibrated** test against
  `H_0: X_j ⫫ Y | X_{~j}`, distinct from (and replacing) the retracted
  0.0.0.9012 pair-bootstrap method. The `total_order_test` flag on the
  result distinguishes the new mode from the retracted one.
- Result list re-populates `p_value_total_order` and
  `p_value_total_order_adjusted` under the new mode. The grid-wide
  `p_adjust` applies as for first-order p-values.
- [`print()`](https://rdrr.io/r/base/print.html) surfaces a
  `min_p_total` column when `total_order_test = "cond_perm"`, with an
  honest annotation naming the null tested.
- Regression test asserts that under pure-noise `Y`, the new method does
  **not** systematically reject — at least one parameter’s raw p-value
  remains ≥ 0.05 (the 0.0.0.9012 failure mode is pinned).
- Defunct-arg error for `total_order_p_value` updated to point at both
  `total_order_ci` (uncertainty) and `total_order_test = "cond_perm"`
  (significance).

##### Test counts

Test suites at L99-coverage: density-ratio fit/predict + diagnostics +
refactored bd-HSIC + cond_perm calibration on additive and pure-noise
designs.

### kernR 0.0.0.9013

#### Pick-Freeze p-value retraction (post critical-review, 2026-05-16)

The 0.0.0.9012 `total_order_p_value = TRUE` mode was found by critical
review to be **not null-calibrated**. The pair-bootstrap samples the
empirical joint distribution, not a null-of-no-effect, so under
pure-noise `Y` every parameter was assigned a tiny p-value
(`p ≈ 1/(1 + B)`) by mechanical bootstrap geometry rather than real
signal. The orchestra smoke output `active min-p = inert min-p = 0.010`
was the visible failure.

Changes:

- `total_order_p_value` arg is **defunct**. Passing any non-`NULL` value
  now errors with a pointer to `total_order_ci`. Code that used the
  field must migrate.
- New arg `total_order_ci = TRUE` activates a pair-bootstrap percentile
  CI on the index `T_j` itself — uncertainty quantification, not a
  hypothesis test.
- Result list loses `p_value_total_order` and
  `p_value_total_order_adjusted`. The CI fields `ci_total_order_lower` /
  `ci_total_order_upper` remain (still valid).
- [`print()`](https://rdrr.io/r/base/print.html) drops the misleading
  `min_p_total` column; instead shows a `T_CI` range when
  `total_order_ci = TRUE`, with an explicit “NOT a significance test”
  caveat.
- [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) loses
  the p-value columns; retains CI columns.
- Test suite picks up a critical-review regression test that asserts no
  significance fields are returned under pure-noise `Y`, plus a
  defunct-arg-error test.
- Docs
  ([`?hsic_sensitivity`](https://max578.github.io/kernR/reference/hsic_sensitivity.md)
  Details + `vignettes/kernR-sensitivity.Rmd`) rewritten to reflect the
  retraction and the still-open future-work item: a properly
  null-calibrated total-order significance test (candidate path:
  conditional-independence rather than the marginal complement
  formulation).

A separate `feedback_total_order_calibration.md` memory entry has been
written to ensure the lesson — “bootstrap-around-empirical is not a
null” — persists across sessions.

### kernR 0.0.0.9012

#### Pick-Freeze bootstrap p-values for total-order HSIC sensitivity

- [`hsic_sensitivity()`](https://max578.github.io/kernR/reference/hsic_sensitivity.md)
  gains three new arguments: `total_order_p_value` (logical, default
  `FALSE` — backwards-compatible), `n_bootstrap` (integer, default
  `200L`), and `ci_level` (numeric, default `0.95`). When activated,
  computes Pick-Freeze pair-bootstrap p-values for the null
  `H_0: T_j = 0` and percentile CIs for each total-order index.
- Result list picks up four new fields when active:
  `p_value_total_order`, `p_value_total_order_adjusted`,
  `ci_total_order_lower`, `ci_total_order_upper`. The grid-wide
  `p_adjust` method (default BH) applies to the total-order p-value grid
  as well as the first-order grid.
- [`print()`](https://rdrr.io/r/base/print.html) method surfaces a
  `min_p_total` column when total-order p-values are present;
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) exposes
  the four new fields as long-format columns alongside the existing
  first-order ones.
- Closes the v0.0.0.9006 deferred item (“Total-order p-values explicitly
  deferred and documented”) from the B5 total-order extension.
  Permutation null is not the right inference for `T_j = 0` — under
  marginal independence `T_j` is `1`, not `0` — hence the bootstrap
  formulation (Da Veiga 2015 §4).
- New tests at `tests/testthat/test-sensitivity.R` (8 new blocks):
  population, backwards-compat, validation, reproducibility,
  additive-case CI claim, kwarg validation.

### kernR 0.0.0.9011

#### Cross-package contract formalisation

- `proxymix (>= 0.3.0)` formally declared in `Suggests:` to match the
  optional density-ratio backend already wired in 0.0.0.9010. Closes the
  DESCRIPTION-vs-NEWS gap noted during the 2026-05-16 cross-workspace
  orchestra audit.
- New vignette `kernR-proxymix-binding.Rmd` — single-screen demo of all
  four
  [`estimate_density_ratio()`](https://max578.github.io/kernR/reference/estimate_density_ratio.md)
  backends (`logistic`, `ranger`, `xgboost`, `proxymix`) on one toy
  confounded problem, with ESS, weight-range, and bd-HSIC p-value
  tabulated side-by-side. Chunks are guarded by
  [`requireNamespace("proxymix")`](https://github.com/max578/proxymix);
  the vignette renders unconditionally and degrades gracefully when
  proxymix is unavailable.
- No behavioural change to any exported function. Soft dependency only;
  CRAN-bound builds without proxymix remain green.

### kernR 0.0.0.9010

#### proxymix density-ratio backend (UQ ag-stack roadmap §C1)

- [`estimate_density_ratio()`](https://max578.github.io/kernR/reference/estimate_density_ratio.md)
  gains a fourth backend, `method = "proxymix"`. Fits Gaussian-mixture
  proxies (Hoek & Elliott, 2024) to the joint and product-of-marginals
  sample clouds via `proxymix::fit_proxymix(regime = "sample")` and
  computes pointwise density ratios from
  [`proxymix::dgmm()`](https://rdrr.io/pkg/proxymix/man/dgmm.html)
  evaluations. Useful when the underlying densities are multimodal or
  when NCE-classifier calibration is unreliable.
- [`bd_hsic_test()`](https://max578.github.io/kernR/reference/bd_hsic_test.md)
  exposes the new option via `density_ratio = "proxymix"` — the
  cross-package wedge between kernR (verdict) and proxymix
  (density-ratio bridge) in the UQ ag stack.
- New argument `proxymix_components =` (default `2L`) sets the
  Gaussian-mixture component count per density.
- `proxymix (>= 0.3.0)` added to `Suggests` (soft dependency,
  GRDC-firewalled, MIT).
  [`requireNamespace()`](https://rdrr.io/r/base/ns-load.html) guard in
  the dispatch emits a clear install hint when proxymix is unavailable.
- New tests at `tests/testthat/test-density-ratio-proxymix.R`
  (`skip_if_not_installed("proxymix")` so CRAN’s farm and downstream
  users without proxymix are unaffected).
- DR-DATE /
  [`dr_date_scenario()`](https://max578.github.io/kernR/reference/dr_date_scenario.md)
  still uses logistic / ranger / xgboost propensity backends.
  Proxymix-via-propensity is a different statistical task (modelling
  `P(T = 1 | X)` rather than a density ratio between two sample clouds)
  and is tracked as future work.

#### Contract symmetry: `mmd_ppc()` for the PESTO 0.3.0 manifest

- New
  [`mmd_ppc.pesto_ensemble_manifest()`](https://max578.github.io/kernR/reference/mmd_ppc.md)
  S3 method —
  [`mmd_ppc()`](https://max578.github.io/kernR/reference/mmd_ppc.md) now
  consumes a
  [`PESTO::pesto_ensemble_manifest`](https://rdrr.io/pkg/PESTO/man/pesto_ensemble_manifest.html)
  directly, completing the v0.3.0 cross-package contract symmetry
  alongside
  [`dr_date_scenario()`](https://max578.github.io/kernR/reference/dr_date_scenario.md).
  Posterior-predictive sample comes from `m@outputs`; the user must
  supply held-out `observed` (the manifest’s `obs_target` slot is a
  single nobs-dim point — the data the posterior was fit to — and is
  unsuitable as a two-sample comparator).
- New `outputs =` argument (parallel to
  [`dr_date_scenario()`](https://max578.github.io/kernR/reference/dr_date_scenario.md)’s
  convention) lets the user focus the check on specific observation
  columns.
- `R/zzz.R` `.onLoad()` now calls
  [`registerS3method()`](https://rdrr.io/r/base/ns-internal.html) to
  wire the method up under PESTO’s package-qualified S7 class string
  (`"PESTO::pesto_ensemble_manifest"`) — standard
  [`UseMethod()`](https://rdrr.io/r/base/UseMethod.html) dispatch can’t
  reach it via the bare-name function file name because R can’t parse
  `::` in an S3-method identifier.
- `kernR-ppc` vignette extended with a cross-package handoff section
  demonstrating both true out-of-sample and retrodictive use against the
  new contract.

### kernR 0.0.0.9009

#### Fixes to §B2 (`dr_date_scenario()`) shipped earlier in 0.0.0.9007

- Class detection: S7 sets the S3 class attribute to the
  **package-qualified** `"PESTO::pesto_ensemble_manifest"` (plus
  `"S7_object"`), not the bare class name. `.validate_manifest_pair()`
  now accepts both qualified and bare forms via a new
  `.is_pesto_manifest()` helper.
- `PESTO` is now referenced via
  `@importFrom PESTO pesto_ensemble_manifest as_manifest` in
  `R/kernR-package.R`, clearing the “Namespace in Imports field not
  imported from: PESTO” NOTE without moving the dep to Suggests.
- Test helper rewritten: the canonical ag-scenario use case is **one
  PESTO posterior, forward-simulated under two scenarios** (same params,
  different outputs) — not two independent PESTO calibrations on
  different data (which would produce divergent posteriors and violate
  DR-DATE’s positivity assumption, returning `p ≈ 1`). Test and vignette
  now both reflect the canonical pattern, with a documented escape hatch
  for the different-posteriors case.
- `R CMD check --as-cran` now returns 0 errors / 0 notes / 2 baseline
  WARNINGs (unchanged local-env structural — gate held).

### kernR 0.0.0.9008

#### New features — kernel-based downscaling

Two complementary downscaling primitives, both rooted in the
Park-Muandet-Fukumizu-Sejdinovic / Szabó-Sriperumbudur-Póczos-Gretton
family of RKHS regression methods.

- `kernel_downscale(coarse, fine, new_coarse, ...)`: vector-in /
  vector-out kernel-based downscaling via conditional mean embedding
  (Park, Muandet, Fukumizu, Sejdinovic 2013). Trains a CME on paired
  `(coarse, fine)` data and predicts fine-resolution outputs at new
  coarse inputs. Multi-output `fine` supported. Lambda by LOO-CV by
  default. S3: `print`, `as.data.frame`. Returns the `n_new x n_train`
  weight matrix on demand.
- `dist_regression(bags, y, ...)`: bag-in / vector-out distribution
  regression (Szabó, Sriperumbudur, Póczos, Gretton 2016). Each input is
  a *bag of points* mapped to its empirical mean embedding;
  `outer = "linear"` (inner-product of embeddings) or `"rbf"` (Gaussian
  over embedding-space distance). Variable bag sizes supported.
  Multivariate `y` supported.
  [`predict()`](https://rdrr.io/r/stats/predict.html) method for
  out-of-bag prediction. S3: `print`, `as.data.frame`.
- [`fit_cme()`](https://max578.github.io/kernR/reference/fit_cme.md) and
  [`predict.cme_fit()`](https://max578.github.io/kernR/reference/predict.cme_fit.md)
  are now exported (previously `@keywords internal`); they remain the
  lower-level building block, with
  [`kernel_downscale()`](https://max578.github.io/kernR/reference/kernel_downscale.md)
  as the user-facing wrapper.
- New vignette `kernR-downscaling`: covers both methods with worked
  ag-systems examples (coarse climate → paddock yield; bag-of-soil-cores
  → paddock yield).

### kernR 0.0.0.9007

#### New features (UQ ag-stack roadmap §B2)

- `dr_date_scenario(baseline, intervention, ...)` — DR-DATE
  distributional treatment-effect test for the **two-scenario APSIM use
  case**, where `baseline` and `intervention` are
  [`PESTO::pesto_ensemble_manifest`](https://rdrr.io/pkg/PESTO/man/pesto_ensemble_manifest.html)
  objects (the v0.3.0 cross-package S7 contract). Pools parameters as
  covariates, outputs as the outcome, scenario label as binary
  treatment; dispatches to the existing
  [`dr_date_test()`](https://max578.github.io/kernR/reference/dr_date_test.md)
  machinery. Returns a `dr_date_scenario` (subclass of
  `kernel_test_result`) carrying both run-ids plus a directly actionable
  verdict line.
- `print.dr_date_scenario()` — verdict-focused printer.
- New vignette: `kernR-drdate-scenario` — synthetic linear-Gaussian
  scenario contrast end-to-end (PESTO IES → manifest → DR-DATE).
- `PESTO (>= 0.3.0)` becomes a hard `Imports:` so the S7 contract
  resolves at install time.

#### Notes

- This is the kernR-side counterpart of PESTO §A5 (which shipped the
  `pesto_ensemble_manifest` S7 class). The legacy lightweight
  [`pesto_ensemble()`](https://max578.github.io/kernR/reference/pesto_ensemble.md)
  S3 class for
  [`mmd_ppc()`](https://max578.github.io/kernR/reference/mmd_ppc.md) is
  unchanged; a future release will add a
  [`mmd_ppc.pesto_ensemble_manifest()`](https://max578.github.io/kernR/reference/mmd_ppc.md)
  method for consistency with the new contract.
- Proxymix density-ratio backend (roadmap §C1) will become a fourth
  `propensity_model` option once landed; tracked as future work.

### kernR 0.0.0.9006

#### New features

- [`hsic_sensitivity()`](https://max578.github.io/kernR/reference/hsic_sensitivity.md)
  gains `total_order = FALSE` argument (default; backwards-compatible).
  When `TRUE`, computes total-order indices via Da Veiga’s complement
  formulation
  `T_j = 1 - HSIC(X_{~j}, Y) / sqrt(HSIC(X_{~j}, X_{~j}) HSIC(Y, Y))`,
  where `X_{~j}` is the design matrix with column `j` removed.
  `T_j - S_j` (returned in the result and shown in `print`) quantifies
  the contribution of `X_j` through *interactions* with other
  parameters. Result gains `index_total_order` (`p x q` matrix),
  `statistic_total_order`, and `total_order` (logical flag) fields.
- [`plot.hsic_sensitivity()`](https://max578.github.io/kernR/reference/plot.hsic_sensitivity.md)
  gains `which = c("first", "total", "both")` argument; `"both"`
  produces side-by-side first-vs-total bars.
- `as.data.frame.hsic_sensitivity()` includes total-order columns when
  present.
- B5 vignette `kernR-sensitivity` extended with a pure-interaction
  example (Y = X1 \* X2) demonstrating S ~ 0 but T strong; and a
  near-additive contrast where T ~ S.

#### Notes

- Total-order permutation p-values are intentionally not computed. The
  natural null for `T_j = 0` is conditional independence of `X_j` and
  `Y` given `X_{~j}` – genuinely harder than the marginal-HSIC
  permutation. Indices are interpreted directly; future work may add
  Pick-Freeze-style p-values.
- Nystrom acceleration for total-order is also deferred. The naive
  materialisation `F F^T` is slower than exact computation at typical
  kernR scales; the proper unblock is a factor-only HSIC primitive
  (Nystrom-on-Nystrom), which is now the natural next acceleration item.

### kernR 0.0.0.9005

#### New features

- [`nystrom_factor()`](https://max578.github.io/kernR/reference/nystrom_factor.md):
  Nystrom low-rank kernel factorisation (Williams & Seeger, 2001).
  Returns an `n x m` factor `F` with `F %*% t(F) \approx K` for any
  \[kernel_spec()\]. `O(n m^2)` construction, `O(n m)` storage. Honours
  DESCRIPTION’s Nystrom claim.
- [`rff_features()`](https://max578.github.io/kernR/reference/rff_features.md):
  Random Fourier Features (Rahimi & Recht, 2007) for RBF kernels.
  Returns an `n x D` feature matrix with `Phi %*% t(Phi) \approx K`.
  Data-independent random projection. Honours DESCRIPTION’s RFF claim.
- [`hsic_test_nystrom()`](https://max578.github.io/kernR/reference/hsic_test_nystrom.md):
  Drop-in accelerated HSIC independence test via low-rank factorisation.
  `method = "nystrom"` (default) or `"rff"`. `O(n m_x m_y)` per
  permutation (vs `O(n^2)` for the exact test); verdict-equivalent to
  [`hsic_test()`](https://max578.github.io/kernR/reference/hsic_test.md)
  at moderate `m`. Uses the biased HSIC estimator (the form that factors
  cleanly through low-rank approximations).
- New vignette `kernR-nystrom`: correctness check + scaling benchmark +
  when-to-use guide.

#### Release-gate investigation (no package change)

- Confirmed that both observed `R CMD check --as-cran` WARNINGs on the
  development machine are local-environment artifacts:
  - `-Wno-unused-command-line-argument` and `-mcpu=native` come from the
    user’s global `~/.R/Makevars`; bypassing it via
    `R_MAKEVARS_USER=/dev/null` removes the “compilation flags” WARNING.
  - The install WARNING is from R’s own `R_ext/Boolean.h:62` using a
    `#pragma clang diagnostic ignored "-Wfixed-enum-extension"` that
    bleeding-edge Apple clang 21.0.0 (MacOSX26.4.1.sdk) does not
    recognise. Not a kernR issue; will not appear on CRAN’s build farm.
    Documented in PROJECT_LOG.

### kernR 0.0.0.9004

#### New features

- [`bd_hsic_test()`](https://max578.github.io/kernR/reference/bd_hsic_test.md)
  gains `cluster_id =` and `permutation =` arguments, implementing
  roadmap item **B4** (hierarchical bd-HSIC). When `cluster_id` is
  supplied, the permutation null is built by within-cluster reshuffling
  of `y` (the safer default for clustered ag designs); the result
  carries `permutation_scheme`, `cluster_id`, `cluster_levels`, and a
  `per_cluster_statistic` stratified breakdown. Backwards-compatible:
  `cluster_id = NULL` preserves the original Hu/Sejdinovic/Evans
  propensity-cluster behaviour.
- [`hsic_sensitivity()`](https://max578.github.io/kernR/reference/hsic_sensitivity.md):
  First-order HSIC-Sensitivity Index per Da Veiga 2015 – normalised HSIC
  bounded in `[0, 1]`, Sobol-comparable in scale, captures
  distributional effects Sobol misses (variance, skewness, tails). S3
  methods: `print`, `plot`, `as.data.frame`. Implements roadmap item
  **B5** (first-order; total-order intentionally deferred until Nystrom
  acceleration lands). Reuses the per-column kernel-matrix cache pattern
  from
  [`hsic_identifiability()`](https://max578.github.io/kernR/reference/hsic_identifiability.md).
- New vignettes:
  - `kernR-hierarchical-bdhsic`: multi-site stub design comparing naive
    vs within-cluster permutation; per-site stratified contributions.
  - `kernR-sensitivity`: HSIC-SI on a stub APSIM-like simulator with
    mean / variance / tail / inert effects; demonstrates the
    distributional-effect catch.

#### Release-gate cleanup

- Added `@importFrom stats weighted.mean` (clears NOTE on missing global
  function definition).
- Escaped `{` / `}` in math expressions in three Rcpp roxygen blocks
  (`mmd2_unbiased_cpp`, `rulsif_solve_cpp`, `weighted_hsic_stat_cpp`)
  using `\deqn{...}{...}` form (clears NOTE on Rd `Lost braces`).
- `R CMD check --as-cran` is now **0 errors / 2 WARNINGs / 0 NOTES**;
  the remaining WARNINGs are the toolchain-level non-portable Makevars
  flags (`-Wno-unused-command-line-argument`, `-mcpu=native`), restored
  at the FLIBS portability release-gate ritual per workspace CLAUDE.md.

### kernR 0.0.0.9003

#### New features

- [`mmd_ppc()`](https://max578.github.io/kernR/reference/mmd_ppc.md):
  Posterior-predictive check via MMD two-sample test. Returns the
  standard `kernel_test_result` plus a Shannon-information *surprise*
  diagnostic (`-log2(p_value)`) and an explicit reject/accept verdict at
  level `alpha`. Implements roadmap item **B3**.
- [`pesto_ensemble()`](https://max578.github.io/kernR/reference/pesto_ensemble.md):
  Lightweight constructor + S3 class (`pesto_ensemble`) bundling a
  posterior-predictive sample matrix, optional held-out observations,
  and free-form metadata. Defines the kernR-side schema of the
  cross-package PESTO -\> kernR contract until PESTO’s native emitter
  lands.
  [`mmd_ppc()`](https://max578.github.io/kernR/reference/mmd_ppc.md)
  dispatches on it.
- New vignette `kernR-ppc`: walkthrough on a stubbed ensemble across
  calibrated / narrow-variance / mean-shifted scenarios.

### kernR 0.0.0.9002

#### New features

- [`hsic_identifiability()`](https://max578.github.io/kernR/reference/hsic_identifiability.md):
  HSIC-based pre-IES screening that flags unidentifiable APSIM (or any
  simulator) parameters before ensemble-smoother calibration. Returns a
  `p x q` grid of HSIC statistics + permutation p-values across
  parameters and outputs, with Benjamini-Hochberg adjustment by default.
  S3 methods: `print`, `summary`, `plot`, `as.data.frame`. Implements
  roadmap item **B1**.
- [`lhs_design()`](https://max578.github.io/kernR/reference/lhs_design.md):
  Lightweight Latin-hypercube design helper over bounded parameter
  ranges; base-R only.
- New vignette `kernR-identifiability`: walkthrough of the pre-IES
  screening workflow on a stubbed APSIM archetype.

### kernR 0.0.0.9001

#### Authorship & administrative

- Sole authorship by Max Moldovan (`aut`, `cre`, `cph`). D. Sejdinovic
  removed from `Authors@R`; paper citations to Hu/Sejdinovic/Evans (JMLR
  2024. and Fawkes/Hu/Evans/Sejdinovic (TMLR 2024) retained throughout
        vignettes, man pages, and source comments as scientific
        attribution.
- Maintainer email corrected to `max.moldovan@adelaide.edu.au` (was
  placeholder).
- Vignette `author:` YAML headers updated.

### kernR 0.0.0.9000

- Initial development version.
- **Kernel engine**: RBF, Matern, linear, polynomial kernels with
  RcppArmadillo backend.
- **Bandwidth selection**: Median heuristic (Rcpp), Scott’s rule.
- **Base tests**:
  [`hsic_test()`](https://max578.github.io/kernR/reference/hsic_test.md)
  (independence),
  [`mmd_test()`](https://max578.github.io/kernR/reference/mmd_test.md)
  (two-sample).
- **Causal tests**:
  - [`bd_hsic_test()`](https://max578.github.io/kernR/reference/bd_hsic_test.md):
    Backdoor-adjusted HSIC for causal association testing.
  - [`dr_date_test()`](https://max578.github.io/kernR/reference/dr_date_test.md):
    Doubly robust distributional average treatment effect.
  - [`dr_dett_test()`](https://max578.github.io/kernR/reference/dr_dett_test.md):
    Doubly robust distributional effect on the treated.
- **Hierarchical**:
  [`hierarchical_test()`](https://max578.github.io/kernR/reference/hierarchical_test.md)
  for nested/clustered data with within/between decomposition.
- **Unified interface**:
  [`kernel_causal_test()`](https://max578.github.io/kernR/reference/kernel_causal_test.md)
  with formula syntax `y ~ treatment | confounders`.
- **Density ratio estimation**: Logistic NCE via
  `glm`/`ranger`/`xgboost`; RuLSIF (kernel-based).
- **Propensity scores**: Cross-fitted estimation with trimming and
  diagnostics.
- **Diagnostics**:
  [`assess_overlap()`](https://max578.github.io/kernR/reference/assess_overlap.md),
  [`plot_weights()`](https://max578.github.io/kernR/reference/plot_weights.md),
  [`effective_sample_size()`](https://max578.github.io/kernR/reference/effective_sample_size.md).
- **Vignettes**: Quick start, bd-HSIC tutorial, DR-DATE/DETT tutorial,
  hierarchical data.
