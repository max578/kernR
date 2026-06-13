# kernR 0.8.1

## Bug fixes

* The `.validate_manifest_pair()` obs-schema test (`test-drdate-scenario.R`)
  now guards on whether the installed PESTO actually exports the development
  symbol `pesto_obs_schema`, rather than on a PESTO version string. A `.9000`
  development version cannot distinguish a PESTO build that exports the symbol
  from one that does not, so continuous integration -- which installs PESTO
  from its published `Remotes` source -- erroneously ran the test against a
  build lacking the symbol and reported an error. The test is now skipped where
  the capability is absent, in line with the federation rule that a member's
  tests must not hard-fail on a sibling's unreleased symbol. No package code
  changed.

# kernR 0.8.0

## New features

* New exported `relative_entropy()`, `relative_entropy_ensemble()`, and
  `cir_objective()`: the kernR-owned statistic of Assimilative Causal Inference
  (ACI; Andreou, Chen & Bollt, *Nat. Commun.* 2026). `relative_entropy()`
  computes the closed-form Kullback-Leibler divergence between two Gaussians --
  the smoother-versus-filter posterior contrast that identifies a hidden cause;
  `relative_entropy_ensemble()` is the moment-matched wrapper for posterior
  ensembles; `cir_objective()` reduces a divergence-versus-lag profile to a
  single threshold-free causal influence range. The filter/smoother engine
  lives upstream (kalmix / PESTO); kernR owns the distributional measure and the
  range integral. This release unblocks kalmix v0.2.0, which consumes the ACI
  statistic.

* New `tidy()` methods for kernR verdict objects (`tidy.kernel_test_result()`,
  `tidy.taci_result()`), registered against the `broom`-style `tidy()` generic
  (re-exported from `generics`). They return a stable one-row-per-term
  `data.frame` so downstream code stops reaching into result fields by name.
  The summary uses the `broom`-canonical `p.value` column (with a dot),
  distinct from the result's native `p_value` field (with an underscore) -- the
  exact mismatch that previously tripped two consuming scripts. An `mmd_ppc`
  result additionally tidies its `surprise_bits` and `reject`; a `taci_result`
  threads its `decision` and `grounding` (Independent Oracle Principle).

* `taci_test()` gained `mechanism_provenance` / `posterior_provenance`
  arguments and now labels its verdict's grounding (Independent Oracle
  Principle): TACI builds its entire reference band from the caller-supplied
  `mechanism`, which it cannot itself verify against reality, so the result
  carries a `grounding` token (`"grounded"` only when `mechanism_provenance` is
  declared, else `"[unverified]"`) and a human-facing `verdict` string suffixed
  `[unverified]` when un-grounded. The three-way `decision` enum is unchanged
  (consumers unaffected); the grounding is pure metadata.

* `dr_date_scenario()` now grounds a scenario comparison against the
  Independent Oracle Principle: `.validate_manifest_pair()` refuses a pair built
  against incompatible APSIM major versions (a silent coefficient-drift hazard
  behind an identical schema) and, when both manifests carry a
  `pesto_ensemble_manifest` `obs_schema` (PESTO schema 1.1.0+), refuses when a
  shared output column disagrees on unit or quantity. Opt-in-by-presence: a
  pre-1.1.0 manifest with no `obs_schema` skips the unit check gracefully.

* New exported `taci_test()` implements theory-anchored causal inference: a
  mechanism-consistency test that asks whether an observed treatment effect
  agrees with the effect a calibrated mechanistic model predicts. The reference
  distribution is built from the model's own posterior-predictive draws rather
  than a permutation null, and the statistic is the same weighted bd-HSIC the
  `bd_hsic_test()` engine uses, so the verdict is three-way -- consistent with
  the model-implied effect, an effect the model does not predict, or no effect.
  A posterior-adequacy guard flags a degenerate H1 band when the posterior pins
  the model-implied effect too precisely. Backdoor adjustment reuses the
  density-ratio machinery; binary and continuous treatments are both supported.

* New exported helpers `weighted_hsic_stat()` and `resolve_bandwidth()` expose
  the kernel-matrix and weighted-HSIC primitives that underpin `bd_hsic_test()`.
  `weighted_hsic_stat()` is a validated R wrapper over the internal compiled
  engine; `resolve_bandwidth()` fills a `kernel_spec`'s median-heuristic
  bandwidth from data. Together they let downstream methods (for example
  theory-anchored causal inference) assemble a bd-HSIC-compatible statistic
  against a custom reference distribution without reaching into kernR internals.

* `bd_hsic_test()` now gates its verdict on the `proxymix` density-ratio
  backend's fit quality (C4). `fit_density_ratio(method = "proxymix")` surfaces
  a single `fit_quality` pass-through flag (`ok` / `status` / `reason`) reduced
  from per-GMM convergence; when a mixture proxy fails to converge,
  `bd_hsic_test()` warns and sets `result$density_ratio_warning = TRUE` (`FALSE`
  for clean fits and all other backends). Mirrors the existing ESS-floor gate:
  an untrustworthy verdict is flagged, never reported silently.

* New `joint_coverage_test()`: a joint (multivariate) calibration diagnostic.
  Where `coverage_test()` assesses each output dimension separately (marginal
  calibration) and is blind to the dependence structure, this builds a
  multivariate rank histogram -- each held-out observation is reduced, with the
  ensemble, to one multivariate rank through a pre-rank function, and the
  histogram is tested for uniformity. The default band-depth pre-rank
  (Thorarinsdottir et al. 2016) is sensitive to correlation / dependence
  miscalibration that marginal and average-rank methods miss; an ensemble with
  correct margins but a wrong correlation -- which `coverage_test()` passes -- is
  caught here. The average-rank pre-rank (Gneiting et al. 2008) gives the
  familiar U-shape = under-dispersed reading. This is the natural validation for
  PESTO 0.6.0's covariance inflation / localisation: it asks whether the
  ensemble's *covariance*, not just its margins, is calibrated.

## Minor improvements and fixes

* `bd_hsic_test()` no longer imposes a hard floor of 20 observations. Small-N
  field trials (a handful of paddocks by seasons) are now accepted down to the
  genuine mathematical minimum of `6` -- the train/test split must leave at
  least two test observations for the weighted HSIC and two propensity
  clusters. The statistical risk a small sample carries is surfaced by the
  existing ESS-floor reliability gate (a warning), not refused outright.

* Documented the verdict-object p-value accessor explicitly. The field on a
  `kernel_test_result` (and on an `mmd_ppc`) is `p_value` with an underscore,
  **not** `p.value`; the `?mmd_test` and `?mmd_ppc` Value sections now say so
  and point at `tidy()` for the dotted `p.value` column.

* Added a "declare mechanism provenance" worked example to `?taci_test`,
  showing how naming where the mechanism's calibration came from
  (`mechanism_provenance`) moves a verdict's `grounding` from `"[unverified]"`
  to `"grounded"` (Independent Oracle Principle).

## Dependencies

* Raise the `PESTO` dependency floor to `>= 0.6.0`. PESTO 0.6.0 adds covariance
  inflation and localisation against ensemble under-dispersion and records
  per-iteration spread-ESS / inflation / localisation diagnostics in the
  ensemble manifest; the C2 manifest-consumption surface
  (`dr_date_scenario()`, `mmd_ppc()`, `coverage_test()`) is verified against it.

* New `Imports: generics` -- a small (`methods`-only) dependency supplying the
  `tidy()` generic the new tidiers register against, so `broom::tidy()` and
  `generics::tidy()` both dispatch to them.

# kernR 0.7.0

## New features

* New `concordance_test_nystrom()`: a low-rank counterpart to
  `concordance_test()` for large ensembles. The pooled sample is factorised
  once -- by the Nystrom method (default) or random Fourier features -- and the
  summed pairwise unbiased MMD and its joint-permutation null are computed from
  the `n x m` factor in `O(n m)` per permutation rather than `O(n^2)`, with the
  rank `m` controlling the speed / accuracy trade-off. The verdict object and
  the pairwise discrepancy matrix that localises the offending source are
  identical to the exact test.

* New `ksd_test_nystrom()`: a low-rank counterpart to `ksd_test()` for large
  samples. The `n x n` Stein-kernel matrix is never materialised; it is
  replaced by a Nystrom factor with which the kernel Stein discrepancy and its
  wild-bootstrap null are computed in `O(n m)`. Because the Langevin Stein
  kernel is itself positive semi-definite, the factorised kernel is a valid
  Stein kernel in its own right, so the accelerated test is the exact procedure
  applied to a rank-`m` kernel: the degenerate-U-statistic calibration is
  preserved (empirical Type-I error stays at the nominal level), and the
  approximation trades statistical power -- not test validity -- for speed.

* `ksd_test()` and `concordance_test()` gain an `n_exact_max` argument
  (default `5000`). Above this sample-size ceiling the call is delegated to the
  low-rank counterpart, so a large sample stays tractable without an `O(n^2)`
  blow-up. The switch is never silent: it is announced by a message, recorded
  in the returned object (`approximation` and `m`, so the result is
  reproducible from the object), and escapable with `n_exact_max = Inf` to
  force the exact test at any size.

# kernR 0.6.0

## New features

* New `coverage_test()`: a graded coverage / calibration diagnostic for a
  predictive ensemble against held-out observations. Where `mmd_ppc()` and
  `ksd_test()` give a binary reject/accept verdict, `coverage_test()` reports
  *how* and *which way* an ensemble is mis-calibrated -- empirical coverage at
  nominal interval levels, a signed dispersion ratio (`Var(PIT)/(1/12)`: above
  one is under-dispersion / over-confidence, below one is over-dispersion), a
  bias indicator, and a rank-histogram uniformity test -- classifying the
  ensemble as calibrated, under-dispersed, over-dispersed, or biased. It is the
  interpretable companion to the kernel verdicts, and quantifies the
  ensemble-under-dispersion that the cross-member PESTO validation surfaced
  (e.g. "your 90% intervals cover 56% -- under-dispersed by a factor of 1.9").
  Dispatches on a numeric matrix, a `pesto_ensemble`, or a
  `pesto_ensemble_manifest` (threading fidelity provenance), mirroring
  `mmd_ppc()`. Marginal (per-output) calibration; pure R.

# kernR 0.5.0

## Fidelity-provenance awareness when consuming PESTO manifests

* `dr_date_scenario()` now inspects the `fidelity` slot of the two input
  `PESTO::pesto_ensemble_manifest` objects. A provenance mismatch -- one
  scenario single-fidelity and the other multi-fidelity, or two
  multi-fidelity runs with different stack shapes / final levels -- raises
  a `warning` by default (the two ensembles may sit at different physical
  fidelities, confounding the distributional contrast), escalated to an
  error with the new `strict_fidelity = TRUE` argument. The fidelity
  provenance (`list(baseline, intervention)`) is threaded into the result
  and shown by `print()`.
* `mmd_ppc()` on a manifest records the manifest's `fidelity` provenance in
  `result$pesto_metadata$fidelity`, so a PPC verdict traces back to the
  fidelity the producing ensemble was calibrated at.
* Both checks are forward-compatible: manifests from PESTO versions that do
  not populate the slot read as `NULL` and pass. Full fidelity provenance
  is populated by PESTO's multi-fidelity `pesto_ies_callback()` runs.

# kernR 0.4.0

## New features

* New `ksd_test()`: a one-sample kernel Stein discrepancy goodness-of-fit
  test. Where `mmd_test()` compares two samples, `ksd_test()` compares a
  sample against a *distribution* supplied through its score
  (the gradient of its log density), so the target may be unnormalised and
  no reference sample is needed. The calibration framing is direct -- given
  posterior or ensemble draws and the score of the distribution they claim
  to represent, the test asks whether the draws actually follow it. It is
  sensitive to mean, variance, and tail mis-specification, and calibrates
  the degenerate U-statistic null with a wild bootstrap (Chwialkowski et
  al., 2016).
* The default base kernel is the inverse multi-quadric (IMQ), which detects
  non-convergence in regimes where the Gaussian kernel is blind as dimension
  grows (Gorham & Mackey, 2017); the Gaussian (RBF) base kernel remains
  available via `kernel = "rbf"`.
* New `gaussian_score()`: a score-function factory for a multivariate-normal
  target, for use as the `score` argument of `ksd_test()`.
* New `concordance_test()`: a kernel k-sample concordance test asking whether
  several samples -- posterior draws from different inference engines, or
  scenario ensembles from different simulators -- come from a common
  distribution. The statistic is the summed pairwise Maximum Mean Discrepancy
  with a single joint-permutation null (so the family-wise error is
  controlled), and the result carries the full pairwise discrepancy matrix, so
  a rejection localises *which* source diverges and on which margin.
* New `numeric_score()`: a finite-difference score adapter that turns any
  (possibly unnormalised) log density into the score `ksd_test()` needs, so a
  target need only be expressible as a log-density function, not
  hand-differentiated.
* New vignette *Calibration and concordance: kernR's validation layer* tying
  `ksd_test()`, `concordance_test()`, and `numeric_score()` together.

# kernR 0.3.1

## Testing

* New end-to-end analytical-correctness test (`test-end-to-end.R`). One
  confounded data-generating process is run through the whole pipeline and
  every analytical corner is checked for a meaningful, correct verdict:
  marginal `hsic_test()` power and independence; `bd_hsic_test()` removing a
  purely confounder-induced association while detecting a genuine causal one;
  propensity recovery; `dr_date_test()` power, Type I control, and double
  robustness (AIPW and IPW-only fallback agree); `dr_dett_test()`;
  two-sample `mmd_test()`; Nystrom and RFF agreement with exact HSIC;
  hierarchical within-cluster Type I control; full-pipeline seed
  reproducibility; and permutation-null calibration.
* New coverage for the public density-ratio API and weight diagnostics
  (`test-density-ratio-api.R`): `fit_density_ratio()` /
  `predict_density_ratio()` round-trip and reproducibility, `plot_weights()`,
  and the new `print.cme_fit()` method.

## Minor improvements and fixes

* New `print()` method for `fit_cme()` objects (`print.cme_fit()`).
* `@family` tags added across the exported API so the documentation
  cross-links related functions; `dr_date_scenario()`'s example now runs
  (`\donttest` rather than `\dontrun`).
* `R/sensitivity.R` split: the `hsic_sensitivity()` S3 methods now live in
  `R/sensitivity-methods.R`.
* Dropped the `Remotes: github::max578/PESTO` line now that PESTO (>= 0.4.1)
  is served from the max578 r-universe; the CI `extra-repositories` entry
  resolves it.

# kernR 0.3.0

## Correctness

* **DR-DATE and DR-DETT are now genuinely doubly robust.** `dr_date_test()`
  previously fitted a conditional mean embedding outcome model and then
  discarded it, returning an inverse-probability-weighted statistic
  regardless of `outcome_model` -- so `outcome_model = "krr"` and `"zero"`
  gave identical results despite the documented double-robustness claim.
  The statistic now forms the augmented (AIPW) counterfactual mean
  embeddings, consistent if *either* the propensity or the outcome model
  is correctly specified (Fawkes, Hu, Evans & Sejdinovic, 2024).
* **DR-DETT control reweighting corrected.** The effect-on-the-treated
  control counterfactual was reweighting controls by the inverted
  treatment odds `(1 - e) / e`; it now uses the correct treatment odds
  `e / (1 - e)` with an augmented outcome-model correction on the control
  arm.
* **`seed=` now makes permutation tests reproducible.** The C++ permutation
  routines drew from Armadillo's internal RNG, which ignores `set.seed()`;
  they now draw through R's RNG, so a fixed `seed` reproduces the null
  distribution and p-value of `hsic_test()`, `mmd_test()`,
  `bd_hsic_test()`, `dr_date_test()`, `dr_dett_test()` and the
  permutation-based sensitivity paths exactly.
* **`hierarchical_test()` no longer silently swallows within-cluster
  failures.** Failed within-cluster sub-tests are counted and surfaced via
  a warning (or an error if every eligible cluster fails). Clusters too
  small for the chosen sub-test (DR sub-tests need at least 30
  observations) are now skipped with a clear message and *excluded* from
  the within-cluster average, rather than silently contributing zero.
* Integer-overflow hardening in the C++ HSIC normaliser
  (`(double)(n * n)` to `(double)n * n`).

## Features

* `dr_date_test()` and `dr_dett_test()` gain `cross_fit` (default `TRUE`)
  and `n_folds` arguments: both nuisances are cross-fitted and evaluated
  out-of-fold, as the doubly robust theory requires under flexible
  nuisance estimators (Chernozhukov et al., 2018).
* `dr_date_test()` and `dr_dett_test()` gain `min_ess_fraction` and now
  report the per-arm effective sample size (`ess`, `ess_warning`),
  warning when the inverse-probability weights collapse.
* `estimate_propensity()` gains a `seed` argument for reproducible
  cross-fitting folds.

## Documentation and packaging

* Dropped unused `Suggests` (`future`, `future.apply`, `ggplot2`,
  `viridisLite`); raised the `PESTO` floor to the tested `>= 0.4.1`.
* Removed an unpublishable talk citation from `aggregate_downscale()`.
* Corrected a mislabelled "two-sided" comment on the (correct, one-sided
  upper-tail) HSIC permutation p-value.

# kernR 0.2.0

First publish of the local development line to AAGI-AUS. Lands the
0.0.0.9001 → 0.0.0.9015 cycle on top of the existing AAGI v0.1.x
lineage; subsequent 0.0.0.x history below is retained as the
per-feature changelog. The version-number jump (origin v0.1.1 →
v0.2.0) signals the substantial new public surface arriving in this
release.

## Public surface (31 exports)

- **Kernel primitives**: `kernel_spec()`, `kernel_matrix()`,
  `select_bandwidth()`.
- **Independence and two-sample tests**: `hsic_test()`,
  `mmd_test()`, `hsic_test_nystrom()`.
- **Causal association**: `bd_hsic_test()` (backdoor-adjusted HSIC;
  Hu, Sejdinovic & Evans, JMLR 2024).
- **Distributional treatment effects**: `dr_date_test()`,
  `dr_dett_test()`, `dr_date_scenario()`, `kernel_causal_test()`
  (Fawkes, Hu, Evans & Sejdinovic, TMLR 2024).
- **Hierarchical / clustered designs**: `hierarchical_test()` with
  within-cluster permutation.
- **Sensitivity and identifiability**: `hsic_identifiability()`,
  `hsic_sensitivity()` (Da Veiga 2015; conditional-permutation null
  for total-order significance).
- **Density-ratio and propensity**: `fit_density_ratio()`,
  `predict_density_ratio()`, `estimate_density_ratio()` (logistic /
  ranger / xgboost / proxymix backends); `estimate_propensity()`,
  `assess_overlap()`, `plot_weights()`, `effective_sample_size()`.
- **Low-rank acceleration**: `nystrom_factor()` (Williams & Seeger
  2001), `rff_features()` (Rahimi & Recht 2007).
- **Kernel downscaling and distribution regression**:
  `kernel_downscale()` (Park, Muandet, Fukumizu & Sejdinovic 2013),
  `fit_cme()`, `dist_regression()` (Szabó, Sriperumbudur, Póczos &
  Gretton 2016), `aggregate_downscale()`,
  `posterior_sample_aggregate()`.
- **Posterior-predictive check + PESTO contract**: `mmd_ppc()`
  (consumes `PESTO::pesto_ensemble_manifest` via S3 dispatch).
- **Design**: `lhs_design()` Latin-hypercube helper.

## Cross-package contracts

- Imports `PESTO (>= 0.3.0)` for the `pesto_ensemble_manifest` S7
  class.
- Optional `proxymix (>= 0.3.0)` in `Suggests` as a density-ratio
  backend (`requireNamespace()`-guarded).

## Documentation

- 12 vignettes covering: quick start, bd-HSIC tutorial, DR-DATE/
  DETT tutorial, hierarchical data, HSIC identifiability, HSIC
  sensitivity, MMD posterior-predictive check, DR-DATE scenario
  via PESTO manifest, hierarchical bd-HSIC, Nyström acceleration,
  proxymix binding, kernel downscaling.

## R CMD check posture

`R CMD check --as-cran` on the AAGI / CI Linux environment is
expected to report 0 errors / 0 warnings / a small number of
environmental NOTEs (new-submission boilerplate, HTML Tidy version
on macOS). Two Apple-clang-21 toolchain WARNINGs surface only on
the maintainer's local machine (R's own `R_ext/Boolean.h` and a
personal `~/.R/Makevars`); both are absent on CRAN's build farm
and on standard GitHub Actions runners.

## Historical development log

All historical entries below are retained for full traceability
(0.0.0.9001 sole-authorship consolidation → 0.0.0.9015 third
downscaling method).

## kernR 0.0.0.9015

### Third downscaling method: `aggregate_downscale()`

Closes the orchestra-completion gap surfaced after the 2026-05-16
deferral closeout: a downscaling method for the *no-paired-training-
data* regime, where only the aggregate observation and a known
aggregator are available. Companion to the proxymix Tier-2 stub
`from_aggregate_likelihood()`, hosted on the kernR (consumption)
side so proxymix's CRAN pre-submission stays untouched.

* New export `aggregate_downscale(y, aggregator, latent_prior, ...)`.
  Dispatches on the aggregator's class:
  * **Linear matrix `A`** — closed-form per-component Kalman update
    plus mixture-weight reweighting by per-component evidence
    `N(y | A mu_k, A Sigma_k A^T + sigma_y^2 I)`.
  * **Function `T(x)`** — importance sampling within each prior
    component (`n_samples_per_component` draws; ESS-floor
    `min_ess_fraction` reliability gate with explicit `warning()`
    on collapse).
* New export `posterior_sample_aggregate()` for drawing from the
  posterior mixture (downstream uncertainty propagation).
* Accepts the latent prior either as a list
  `(means, covariances, weights)` or as any object exposing those
  slots — including `proxymix::fit_proxymix()` results (gated via
  `methods::slot()`, no hard dependency added).
* Vignette `kernR-downscaling.Rmd` reframed from "Two flavours" to
  "Three flavours"; full comparison table.
* 8 new test blocks at `tests/testthat/test-aggregate-downscale.R`:
  single-component closed-form recovery (Kalman update exact);
  two-component mixture reweighting toward the likely cluster;
  non-linear IS path runs with valid moments and p.s.d. covariance;
  reproducibility under seed; ESS-floor warning fires under
  collapse; `posterior_sample_aggregate()` recovers the posterior
  mean; input validation; proxymix `gmm_fit` prior accepted via
  slot extraction (skipped when proxymix absent).

This is the **third** kernR downscaling method, structurally
distinct from `kernel_downscale()` (CME, paired training data) and
`dist_regression()` (bag-of-points, distribution-level regression).
Each occupies a different cell of the (training-data-shape ×
aggregator-knowledge) matrix.

## kernR 0.0.0.9014

### Deferral closeout (post 2026-05-16 critical review)

Closes the three kernR-side tickets named as deferred in
`DRAINSTORMING/reports/orchestra_critical_review_response_2026-05-16.md`:
density-ratio fit/predict refactor, backend diagnostics, and
properly null-calibrated total-order significance test.

#### Density-ratio fit/predict refactor (closes P0 #2)

* New exports: `fit_density_ratio(x, z, method, ...)` and
  `predict_density_ratio(object, new_x, new_z, type = c("log_ratio",
  "weight", "ratio"))`. All four backends (`logistic`, `ranger`,
  `xgboost`, `proxymix`) return a `density_ratio_fit` object that can
  be applied to held-out rows.
* `bd_hsic_test()` now **honours the documented train/test split**:
  fits density-ratio on the training half, predicts on the held-out
  test half. The runtime warning about the sample-split leak from
  0.0.0.9013 is therefore retired. Closes P0 #2 of the critical
  review.
* Internal ratio computation is now in **log-space** end-to-end
  (`log_ratio = log(p_joint) - log(p_marg)` for proxymix;
  `log_ratio = log(p) - log(1-p) + log(n_noise)` for classifier
  backends). Numerically stable on extreme tails.
* `bd_hsic_test()$density_ratio_fit` carries the fitted model so
  callers can inspect backend diagnostics.
* `estimate_density_ratio()` retained as a thin
  backwards-compatible wrapper that fits and predicts on the same
  data; new code should prefer the explicit fit/predict pair.

#### Backend diagnostics for proxymix (closes P1 #3)

* `fit_density_ratio(method = "proxymix")` now surfaces per-GMM
  convergence diagnostics on `fit$diagnostics`:
  `joint_converged`, `marg_converged`, `joint_loglik`, `marg_loglik`,
  `joint_bic`, `marg_bic`, `joint_aic`, `marg_aic`,
  `joint_iterations`, `marg_iterations`, `n_components`.
* `print.density_ratio_fit()` summarises components + BIC + per-GMM
  convergence on screen.

#### Conditional-permutation total-order significance (closes P0 #1 remaining work)

* New arg `total_order_test = c("none", "cond_perm")`. When
  `"cond_perm"`: cluster the design points by `X_{~j}` similarity
  (k-means, `n_clusters_cp` bins; `"auto"` chooses
  `min(floor(n / 5), 20)`); within each cluster permute `Y`;
  recompute `T_j` on each of `n_permutations` permuted designs;
  report `p = (1 + #{T_perm ≥ T_obs}) / (1 + n_permutations)`.
* This is a **properly null-calibrated** test against
  `H_0: X_j ⫫ Y | X_{~j}`, distinct from (and replacing) the
  retracted 0.0.0.9012 pair-bootstrap method. The `total_order_test`
  flag on the result distinguishes the new mode from the retracted
  one.
* Result list re-populates `p_value_total_order` and
  `p_value_total_order_adjusted` under the new mode. The grid-wide
  `p_adjust` applies as for first-order p-values.
* `print()` surfaces a `min_p_total` column when `total_order_test
  = "cond_perm"`, with an honest annotation naming the null tested.
* Regression test asserts that under pure-noise `Y`, the new method
  does **not** systematically reject — at least one parameter's raw
  p-value remains ≥ 0.05 (the 0.0.0.9012 failure mode is pinned).
* Defunct-arg error for `total_order_p_value` updated to point at
  both `total_order_ci` (uncertainty) and
  `total_order_test = "cond_perm"` (significance).

#### Test counts

Test suites at L99-coverage: density-ratio fit/predict + diagnostics
+ refactored bd-HSIC + cond_perm calibration on additive and
pure-noise designs.

## kernR 0.0.0.9013

### Pick-Freeze p-value retraction (post critical-review, 2026-05-16)

The 0.0.0.9012 `total_order_p_value = TRUE` mode was found by
critical review to be **not null-calibrated**. The pair-bootstrap
samples the empirical joint distribution, not a null-of-no-effect,
so under pure-noise `Y` every parameter was assigned a tiny p-value
(`p ≈ 1/(1 + B)`) by mechanical bootstrap geometry rather than real
signal. The orchestra smoke output `active min-p = inert min-p = 0.010`
was the visible failure.

Changes:

* `total_order_p_value` arg is **defunct**. Passing any non-`NULL`
  value now errors with a pointer to `total_order_ci`. Code that
  used the field must migrate.
* New arg `total_order_ci = TRUE` activates a pair-bootstrap
  percentile CI on the index `T_j` itself — uncertainty
  quantification, not a hypothesis test.
* Result list loses `p_value_total_order` and
  `p_value_total_order_adjusted`. The CI fields
  `ci_total_order_lower` / `ci_total_order_upper` remain (still
  valid).
* `print()` drops the misleading `min_p_total` column; instead
  shows a `T_CI` range when `total_order_ci = TRUE`, with an
  explicit "NOT a significance test" caveat.
* `as.data.frame()` loses the p-value columns; retains CI columns.
* Test suite picks up a critical-review regression test that asserts
  no significance fields are returned under pure-noise `Y`, plus a
  defunct-arg-error test.
* Docs (`?hsic_sensitivity` Details + `vignettes/kernR-sensitivity.Rmd`)
  rewritten to reflect the retraction and the still-open future-work
  item: a properly null-calibrated total-order significance test
  (candidate path: conditional-independence rather than the marginal
  complement formulation).

A separate `feedback_total_order_calibration.md` memory entry has
been written to ensure the lesson — "bootstrap-around-empirical is
not a null" — persists across sessions.

## kernR 0.0.0.9012

### Pick-Freeze bootstrap p-values for total-order HSIC sensitivity

* `hsic_sensitivity()` gains three new arguments: `total_order_p_value`
  (logical, default `FALSE` — backwards-compatible), `n_bootstrap`
  (integer, default `200L`), and `ci_level` (numeric, default `0.95`).
  When activated, computes Pick-Freeze pair-bootstrap p-values for the
  null `H_0: T_j = 0` and percentile CIs for each total-order index.
* Result list picks up four new fields when active:
  `p_value_total_order`, `p_value_total_order_adjusted`,
  `ci_total_order_lower`, `ci_total_order_upper`. The grid-wide
  `p_adjust` method (default BH) applies to the total-order p-value
  grid as well as the first-order grid.
* `print()` method surfaces a `min_p_total` column when total-order
  p-values are present; `as.data.frame()` exposes the four new fields
  as long-format columns alongside the existing first-order ones.
* Closes the v0.0.0.9006 deferred item ("Total-order p-values
  explicitly deferred and documented") from the B5 total-order
  extension. Permutation null is not the right inference for
  `T_j = 0` — under marginal independence `T_j` is `1`, not `0` —
  hence the bootstrap formulation (Da Veiga 2015 §4).
* New tests at `tests/testthat/test-sensitivity.R` (8 new blocks):
  population, backwards-compat, validation, reproducibility,
  additive-case CI claim, kwarg validation.

## kernR 0.0.0.9011

### Cross-package contract formalisation

* `proxymix (>= 0.3.0)` formally declared in `Suggests:` to match the
  optional density-ratio backend already wired in 0.0.0.9010. Closes the
  DESCRIPTION-vs-NEWS gap noted during the 2026-05-16 cross-workspace
  orchestra audit.
* New vignette `kernR-proxymix-binding.Rmd` — single-screen demo of all
  four `estimate_density_ratio()` backends (`logistic`, `ranger`,
  `xgboost`, `proxymix`) on one toy confounded problem, with ESS,
  weight-range, and bd-HSIC p-value tabulated side-by-side. Chunks are
  guarded by `requireNamespace("proxymix")`; the vignette renders
  unconditionally and degrades gracefully when proxymix is unavailable.
* No behavioural change to any exported function. Soft dependency only;
  CRAN-bound builds without proxymix remain green.

## kernR 0.0.0.9010

### proxymix density-ratio backend (UQ ag-stack roadmap §C1)

* `estimate_density_ratio()` gains a fourth backend, `method =
  "proxymix"`. Fits Gaussian-mixture proxies (Hoek & Elliott, 2024)
  to the joint and product-of-marginals sample clouds via
  `proxymix::fit_proxymix(regime = "sample")` and computes pointwise
  density ratios from `proxymix::dgmm()` evaluations. Useful when the
  underlying densities are multimodal or when NCE-classifier
  calibration is unreliable.
* `bd_hsic_test()` exposes the new option via `density_ratio =
  "proxymix"` — the cross-package wedge between kernR (verdict) and
  proxymix (density-ratio bridge) in the UQ ag stack.
* New argument `proxymix_components =` (default `2L`) sets the
  Gaussian-mixture component count per density.
* `proxymix (>= 0.3.0)` added to `Suggests` (soft dependency,
  GRDC-firewalled, MIT). `requireNamespace()` guard in the dispatch
  emits a clear install hint when proxymix is unavailable.
* New tests at `tests/testthat/test-density-ratio-proxymix.R`
  (`skip_if_not_installed("proxymix")` so CRAN's farm and downstream
  users without proxymix are unaffected).
* DR-DATE / `dr_date_scenario()` still uses logistic / ranger /
  xgboost propensity backends. Proxymix-via-propensity is a different
  statistical task (modelling `P(T = 1 | X)` rather than a density
  ratio between two sample clouds) and is tracked as future work.

### Contract symmetry: `mmd_ppc()` for the PESTO 0.3.0 manifest

* New `mmd_ppc.pesto_ensemble_manifest()` S3 method — `mmd_ppc()`
  now consumes a `PESTO::pesto_ensemble_manifest` directly, completing
  the v0.3.0 cross-package contract symmetry alongside
  `dr_date_scenario()`. Posterior-predictive sample comes from
  `m@outputs`; the user must supply held-out `observed` (the manifest's
  `obs_target` slot is a single nobs-dim point — the data the posterior
  was fit to — and is unsuitable as a two-sample comparator).
* New `outputs =` argument (parallel to `dr_date_scenario()`'s
  convention) lets the user focus the check on specific observation
  columns.
* `R/zzz.R` `.onLoad()` now calls `registerS3method()` to wire the
  method up under PESTO's package-qualified S7 class string
  (`"PESTO::pesto_ensemble_manifest"`) — standard `UseMethod()`
  dispatch can't reach it via the bare-name function file name
  because R can't parse `::` in an S3-method identifier.
* `kernR-ppc` vignette extended with a cross-package handoff section
  demonstrating both true out-of-sample and retrodictive use against
  the new contract.

## kernR 0.0.0.9009

### Fixes to §B2 (`dr_date_scenario()`) shipped earlier in 0.0.0.9007

* Class detection: S7 sets the S3 class attribute to the
  **package-qualified** `"PESTO::pesto_ensemble_manifest"` (plus
  `"S7_object"`), not the bare class name. `.validate_manifest_pair()`
  now accepts both qualified and bare forms via a new
  `.is_pesto_manifest()` helper.
* `PESTO` is now referenced via `@importFrom PESTO
  pesto_ensemble_manifest as_manifest` in `R/kernR-package.R`,
  clearing the "Namespace in Imports field not imported from: PESTO"
  NOTE without moving the dep to Suggests.
* Test helper rewritten: the canonical ag-scenario use case is **one
  PESTO posterior, forward-simulated under two scenarios** (same
  params, different outputs) — not two independent PESTO calibrations
  on different data (which would produce divergent posteriors and
  violate DR-DATE's positivity assumption, returning `p ≈ 1`).
  Test and vignette now both reflect the canonical pattern, with a
  documented escape hatch for the different-posteriors case.
* `R CMD check --as-cran` now returns 0 errors / 0 notes / 2
  baseline WARNINGs (unchanged local-env structural — gate held).

## kernR 0.0.0.9008

### New features — kernel-based downscaling

Two complementary downscaling primitives, both rooted in the
Park-Muandet-Fukumizu-Sejdinovic / Szabó-Sriperumbudur-Póczos-Gretton
family of RKHS regression methods.

* `kernel_downscale(coarse, fine, new_coarse, ...)`: vector-in /
  vector-out kernel-based downscaling via conditional mean embedding
  (Park, Muandet, Fukumizu, Sejdinovic 2013). Trains a CME on paired
  `(coarse, fine)` data and predicts fine-resolution outputs at new
  coarse inputs. Multi-output `fine` supported. Lambda by LOO-CV by
  default. S3: `print`, `as.data.frame`. Returns the
  `n_new x n_train` weight matrix on demand.
* `dist_regression(bags, y, ...)`: bag-in / vector-out distribution
  regression (Szabó, Sriperumbudur, Póczos, Gretton 2016). Each input
  is a *bag of points* mapped to its empirical mean embedding;
  `outer = "linear"` (inner-product of embeddings) or `"rbf"`
  (Gaussian over embedding-space distance). Variable bag sizes
  supported. Multivariate `y` supported. `predict()` method for
  out-of-bag prediction. S3: `print`, `as.data.frame`.
* `fit_cme()` and `predict.cme_fit()` are now exported (previously
  `@keywords internal`); they remain the lower-level building block,
  with `kernel_downscale()` as the user-facing wrapper.
* New vignette `kernR-downscaling`: covers both methods with
  worked ag-systems examples (coarse climate → paddock yield;
  bag-of-soil-cores → paddock yield).

## kernR 0.0.0.9007

### New features (UQ ag-stack roadmap §B2)

* `dr_date_scenario(baseline, intervention, ...)` — DR-DATE distributional
  treatment-effect test for the **two-scenario APSIM use case**, where
  `baseline` and `intervention` are `PESTO::pesto_ensemble_manifest`
  objects (the v0.3.0 cross-package S7 contract). Pools parameters as
  covariates, outputs as the outcome, scenario label as binary
  treatment; dispatches to the existing `dr_date_test()` machinery.
  Returns a `dr_date_scenario` (subclass of `kernel_test_result`)
  carrying both run-ids plus a directly actionable verdict line.
* `print.dr_date_scenario()` — verdict-focused printer.
* New vignette: `kernR-drdate-scenario` — synthetic linear-Gaussian
  scenario contrast end-to-end (PESTO IES → manifest → DR-DATE).
* `PESTO (>= 0.3.0)` becomes a hard `Imports:` so the S7 contract
  resolves at install time.

### Notes

* This is the kernR-side counterpart of PESTO §A5 (which shipped the
  `pesto_ensemble_manifest` S7 class). The legacy lightweight
  `pesto_ensemble()` S3 class for `mmd_ppc()` is unchanged; a future
  release will add a `mmd_ppc.pesto_ensemble_manifest()` method for
  consistency with the new contract.
* Proxymix density-ratio backend (roadmap §C1) will become a fourth
  `propensity_model` option once landed; tracked as future work.

## kernR 0.0.0.9006

### New features

* `hsic_sensitivity()` gains `total_order = FALSE` argument
  (default; backwards-compatible). When `TRUE`, computes total-order
  indices via Da Veiga's complement formulation
  `T_j = 1 - HSIC(X_{~j}, Y) / sqrt(HSIC(X_{~j}, X_{~j}) HSIC(Y, Y))`,
  where `X_{~j}` is the design matrix with column `j` removed.
  `T_j - S_j` (returned in the result and shown in `print`) quantifies
  the contribution of `X_j` through *interactions* with other
  parameters. Result gains `index_total_order` (`p x q` matrix),
  `statistic_total_order`, and `total_order` (logical flag) fields.
* `plot.hsic_sensitivity()` gains `which = c("first", "total", "both")`
  argument; `"both"` produces side-by-side first-vs-total bars.
* `as.data.frame.hsic_sensitivity()` includes total-order columns when
  present.
* B5 vignette `kernR-sensitivity` extended with a pure-interaction
  example (Y = X1 * X2) demonstrating S ~ 0 but T strong; and a
  near-additive contrast where T ~ S.

### Notes

* Total-order permutation p-values are intentionally not computed.
  The natural null for `T_j = 0` is conditional independence of `X_j`
  and `Y` given `X_{~j}` -- genuinely harder than the marginal-HSIC
  permutation. Indices are interpreted directly; future work may add
  Pick-Freeze-style p-values.
* Nystrom acceleration for total-order is also deferred. The naive
  materialisation `F F^T` is slower than exact computation at typical
  kernR scales; the proper unblock is a factor-only HSIC primitive
  (Nystrom-on-Nystrom), which is now the natural next acceleration
  item.

## kernR 0.0.0.9005

### New features

* `nystrom_factor()`: Nystrom low-rank kernel factorisation
  (Williams & Seeger, 2001). Returns an `n x m` factor `F` with
  `F %*% t(F) \approx K` for any [kernel_spec()]. `O(n m^2)`
  construction, `O(n m)` storage. Honours DESCRIPTION's Nystrom claim.
* `rff_features()`: Random Fourier Features (Rahimi & Recht, 2007) for
  RBF kernels. Returns an `n x D` feature matrix with
  `Phi %*% t(Phi) \approx K`. Data-independent random projection.
  Honours DESCRIPTION's RFF claim.
* `hsic_test_nystrom()`: Drop-in accelerated HSIC independence test
  via low-rank factorisation. `method = "nystrom"` (default) or
  `"rff"`. `O(n m_x m_y)` per permutation (vs `O(n^2)` for the exact
  test); verdict-equivalent to `hsic_test()` at moderate `m`. Uses the
  biased HSIC estimator (the form that factors cleanly through
  low-rank approximations).
* New vignette `kernR-nystrom`: correctness check + scaling
  benchmark + when-to-use guide.

### Release-gate investigation (no package change)

* Confirmed that both observed `R CMD check --as-cran` WARNINGs on the
  development machine are local-environment artifacts:
  - `-Wno-unused-command-line-argument` and `-mcpu=native` come from
    the user's global `~/.R/Makevars`; bypassing it via
    `R_MAKEVARS_USER=/dev/null` removes the "compilation flags"
    WARNING.
  - The install WARNING is from R's own `R_ext/Boolean.h:62` using a
    `#pragma clang diagnostic ignored "-Wfixed-enum-extension"` that
    bleeding-edge Apple clang 21.0.0 (MacOSX26.4.1.sdk) does not
    recognise. Not a kernR issue; will not appear on CRAN's build
    farm. Documented in PROJECT_LOG.

## kernR 0.0.0.9004

### New features

* `bd_hsic_test()` gains `cluster_id =` and `permutation =` arguments,
  implementing roadmap item **B4** (hierarchical bd-HSIC). When
  `cluster_id` is supplied, the permutation null is built by
  within-cluster reshuffling of `y` (the safer default for clustered ag
  designs); the result carries `permutation_scheme`, `cluster_id`,
  `cluster_levels`, and a `per_cluster_statistic` stratified breakdown.
  Backwards-compatible: `cluster_id = NULL` preserves the original
  Hu/Sejdinovic/Evans propensity-cluster behaviour.
* `hsic_sensitivity()`: First-order HSIC-Sensitivity Index per
  Da Veiga 2015 -- normalised HSIC bounded in `[0, 1]`, Sobol-comparable
  in scale, captures distributional effects Sobol misses (variance,
  skewness, tails). S3 methods: `print`, `plot`, `as.data.frame`.
  Implements roadmap item **B5** (first-order; total-order intentionally
  deferred until Nystrom acceleration lands). Reuses the per-column
  kernel-matrix cache pattern from `hsic_identifiability()`.
* New vignettes:
  - `kernR-hierarchical-bdhsic`: multi-site stub design comparing
    naive vs within-cluster permutation; per-site stratified
    contributions.
  - `kernR-sensitivity`: HSIC-SI on a stub APSIM-like simulator with
    mean / variance / tail / inert effects; demonstrates the
    distributional-effect catch.

### Release-gate cleanup

* Added `@importFrom stats weighted.mean` (clears NOTE on missing global
  function definition).
* Escaped `{` / `}` in math expressions in three Rcpp roxygen blocks
  (`mmd2_unbiased_cpp`, `rulsif_solve_cpp`, `weighted_hsic_stat_cpp`)
  using `\deqn{...}{...}` form (clears NOTE on Rd `Lost braces`).
* `R CMD check --as-cran` is now **0 errors / 2 WARNINGs / 0 NOTES**;
  the remaining WARNINGs are the toolchain-level non-portable Makevars
  flags (`-Wno-unused-command-line-argument`, `-mcpu=native`), restored
  at the FLIBS portability release-gate ritual per workspace
  CLAUDE.md.

## kernR 0.0.0.9003

### New features

* `mmd_ppc()`: Posterior-predictive check via MMD two-sample test.
  Returns the standard `kernel_test_result` plus a Shannon-information
  *surprise* diagnostic (`-log2(p_value)`) and an explicit reject/accept
  verdict at level `alpha`. Implements roadmap item **B3**.
* `pesto_ensemble()`: Lightweight constructor + S3 class
  (`pesto_ensemble`) bundling a posterior-predictive sample matrix,
  optional held-out observations, and free-form metadata. Defines the
  kernR-side schema of the cross-package PESTO -> kernR contract until
  PESTO's native emitter lands. `mmd_ppc()` dispatches on it.
* New vignette `kernR-ppc`: walkthrough on a stubbed ensemble across
  calibrated / narrow-variance / mean-shifted scenarios.

## kernR 0.0.0.9002

### New features

* `hsic_identifiability()`: HSIC-based pre-IES screening that flags
  unidentifiable APSIM (or any simulator) parameters before
  ensemble-smoother calibration. Returns a `p x q` grid of HSIC
  statistics + permutation p-values across parameters and outputs, with
  Benjamini-Hochberg adjustment by default. S3 methods: `print`,
  `summary`, `plot`, `as.data.frame`. Implements roadmap item **B1**.
* `lhs_design()`: Lightweight Latin-hypercube design helper over bounded
  parameter ranges; base-R only.
* New vignette `kernR-identifiability`: walkthrough of the pre-IES
  screening workflow on a stubbed APSIM archetype.

## kernR 0.0.0.9001

### Authorship & administrative

* Sole authorship by Max Moldovan (`aut`, `cre`, `cph`). D. Sejdinovic
  removed from `Authors@R`; paper citations to Hu/Sejdinovic/Evans (JMLR
  2024) and Fawkes/Hu/Evans/Sejdinovic (TMLR 2024) retained throughout
  vignettes, man pages, and source comments as scientific attribution.
* Maintainer email corrected to `max.moldovan@adelaide.edu.au` (was
  placeholder).
* Vignette `author:` YAML headers updated.

## kernR 0.0.0.9000

* Initial development version.
* **Kernel engine**: RBF, Matern, linear, polynomial kernels with RcppArmadillo backend.
* **Bandwidth selection**: Median heuristic (Rcpp), Scott's rule.
* **Base tests**: `hsic_test()` (independence), `mmd_test()` (two-sample).
* **Causal tests**:
  - `bd_hsic_test()`: Backdoor-adjusted HSIC for causal association testing.
  - `dr_date_test()`: Doubly robust distributional average treatment effect.
  - `dr_dett_test()`: Doubly robust distributional effect on the treated.
* **Hierarchical**: `hierarchical_test()` for nested/clustered data with within/between decomposition.
* **Unified interface**: `kernel_causal_test()` with formula syntax `y ~ treatment | confounders`.
* **Density ratio estimation**: Logistic NCE via `glm`/`ranger`/`xgboost`; RuLSIF (kernel-based).
* **Propensity scores**: Cross-fitted estimation with trimming and diagnostics.
* **Diagnostics**: `assess_overlap()`, `plot_weights()`, `effective_sample_size()`.
* **Vignettes**: Quick start, bd-HSIC tutorial, DR-DATE/DETT tutorial, hierarchical data.
