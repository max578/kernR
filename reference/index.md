# Package index

## Kernel primitives

Kernel specifications, Gram matrices, and bandwidth selection.

- [`kernel_spec()`](https://max578.github.io/kernR/reference/kernel_spec.md)
  : Create a Kernel Specification
- [`kernel_matrix()`](https://max578.github.io/kernR/reference/kernel_matrix.md)
  : Compute a Kernel Matrix
- [`select_bandwidth()`](https://max578.github.io/kernR/reference/select_bandwidth.md)
  : Select Kernel Bandwidth

## Independence and two-sample tests

HSIC independence test and MMD two-sample test, with low-rank
acceleration.

- [`hsic_test()`](https://max578.github.io/kernR/reference/hsic_test.md)
  : HSIC Independence Test
- [`mmd_test()`](https://max578.github.io/kernR/reference/mmd_test.md) :
  MMD Two-Sample Test
- [`hsic_test_nystrom()`](https://max578.github.io/kernR/reference/hsic_test_nystrom.md)
  : HSIC Independence Test via Low-Rank Factorisation
- [`plot(`*`<kernel_test_result>`*`)`](https://max578.github.io/kernR/reference/plot.kernel_test_result.md)
  : Plot a Kernel Test Result

## Causal association (bd-HSIC)

Backdoor-adjusted HSIC for testing causal association after confounder
adjustment.

- [`bd_hsic_test()`](https://max578.github.io/kernR/reference/bd_hsic_test.md)
  : Backdoor-HSIC Test for Causal Association

## Distributional treatment effects

Doubly robust kernel tests for distributional treatment effects,
including a PESTO-scenario convenience.

- [`dr_date_test()`](https://max578.github.io/kernR/reference/dr_date_test.md)
  : Doubly Robust Distributional Average Treatment Effect Test (DR-DATE)
- [`dr_dett_test()`](https://max578.github.io/kernR/reference/dr_dett_test.md)
  : Doubly Robust Distributional Effect on the Treated Test (DR-DETT)
- [`dr_date_scenario()`](https://max578.github.io/kernR/reference/dr_date_scenario.md)
  : DR-DATE for Two PESTO Ensemble Scenarios
- [`kernel_causal_test()`](https://max578.github.io/kernR/reference/kernel_causal_test.md)
  : Unified Kernel Causal Test

## Hierarchical and clustered designs

Within-cluster permutation for nested data (farms, hospitals, schools).

- [`hierarchical_test()`](https://max578.github.io/kernR/reference/hierarchical_test.md)
  : Hierarchical Kernel Causal Test

## Sensitivity and identifiability

First- and total-order kernel sensitivity; pre-IES identifiability
screening.

- [`hsic_identifiability()`](https://max578.github.io/kernR/reference/hsic_identifiability.md)
  : HSIC-Based Identifiability Diagnostic
- [`hsic_sensitivity()`](https://max578.github.io/kernR/reference/hsic_sensitivity.md)
  : HSIC-Based Distributional Sensitivity Index
- [`plot(`*`<hsic_identifiability>`*`)`](https://max578.github.io/kernR/reference/plot.hsic_identifiability.md)
  : Plot an HSIC Identifiability Scan
- [`plot(`*`<hsic_sensitivity>`*`)`](https://max578.github.io/kernR/reference/plot.hsic_sensitivity.md)
  : Plot HSIC-Sensitivity Indices

## Density-ratio and propensity

Plug-in density-ratio estimation with four backends; propensity-score
diagnostics.

- [`fit_density_ratio()`](https://max578.github.io/kernR/reference/fit_density_ratio.md)
  : Fit a Density-Ratio Model
- [`predict_density_ratio()`](https://max578.github.io/kernR/reference/predict_density_ratio.md)
  : Predict from a Fitted Density-Ratio Model
- [`estimate_density_ratio()`](https://max578.github.io/kernR/reference/estimate_density_ratio.md)
  : Estimate Density Ratios (backwards-compatible wrapper)
- [`estimate_propensity()`](https://max578.github.io/kernR/reference/estimate_propensity.md)
  : Estimate Propensity Scores
- [`assess_overlap()`](https://max578.github.io/kernR/reference/assess_overlap.md)
  : Assess Propensity Score Overlap
- [`plot_weights()`](https://max578.github.io/kernR/reference/plot_weights.md)
  : Plot Weight Diagnostics
- [`effective_sample_size()`](https://max578.github.io/kernR/reference/effective_sample_size.md)
  : Compute Effective Sample Size

## Low-rank acceleration

Nyström and Random Fourier Feature factorisations for large-n kernels.

- [`nystrom_factor()`](https://max578.github.io/kernR/reference/nystrom_factor.md)
  : Nystrom Low-Rank Kernel Factorisation
- [`rff_features()`](https://max578.github.io/kernR/reference/rff_features.md)
  : Random Fourier Features for the RBF Kernel

## Kernel downscaling and distribution regression

Conditional mean embedding, distribution regression, and aggregate
downscaling.

- [`kernel_downscale()`](https://max578.github.io/kernR/reference/kernel_downscale.md)
  : Kernel-Based Statistical Downscaling
- [`fit_cme()`](https://max578.github.io/kernR/reference/fit_cme.md) :
  Estimate Conditional Mean Embedding via Kernel Ridge Regression
- [`dist_regression()`](https://max578.github.io/kernR/reference/dist_regression.md)
  : Kernel Distribution Regression
- [`aggregate_downscale()`](https://max578.github.io/kernR/reference/aggregate_downscale.md)
  : Aggregate-Likelihood Downscaling
- [`posterior_sample_aggregate()`](https://max578.github.io/kernR/reference/posterior_sample_aggregate.md)
  : Sample from the posterior of an aggregate-downscale fit
- [`predict(`*`<cme_fit>`*`)`](https://max578.github.io/kernR/reference/predict.cme_fit.md)
  : Predict Conditional Mean Embedding Weights at New Points
- [`print(`*`<cme_fit>`*`)`](https://max578.github.io/kernR/reference/print.cme_fit.md)
  : Print a Conditional Mean Embedding Fit
- [`predict(`*`<dist_regression>`*`)`](https://max578.github.io/kernR/reference/predict.dist_regression.md)
  : Predict from a Fitted Distribution Regression Model

## Posterior-predictive check

MMD-based posterior-predictive check; consumes the PESTO cross-package
contract.

- [`mmd_ppc()`](https://max578.github.io/kernR/reference/mmd_ppc.md) :
  MMD Posterior-Predictive Check
- [`pesto_ensemble()`](https://max578.github.io/kernR/reference/pesto_ensemble.md)
  : PESTO Ensemble Manifest (Constructor)

## Design

Latin-hypercube design helper.

- [`lhs_design()`](https://max578.github.io/kernR/reference/lhs_design.md)
  : Latin-Hypercube Design Over Bounded Parameters
