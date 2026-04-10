# kernR (development version)

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
