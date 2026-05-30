# Pre-IES Identifiability Screening with HSIC

## Why screen before calibrating?

Iterative Ensemble-Smoother (IES) calibration of a mechanistic
agricultural model such as APSIM is expensive: each ensemble member is a
full simulator run. Spending realisations on parameters that have no
detectable effect on the outputs of interest is wasted budget.

A pre-IES screen answers a sharper question than variance-based Sobol
analysis: *does this parameter influence the distribution of any output
at all?* HSIC – the Hilbert-Schmidt Independence Criterion – detects
both linear and non-linear, mean-shift and tail-shift effects, including
those Sobol misses when variance is preserved.

The workflow is three lines:

1.  Build a Latin-hypercube design across plausible parameter ranges.
2.  Run the simulator on the design.
3.  Pass the design + outputs to
    [`hsic_identifiability()`](https://max578.github.io/kernR/reference/hsic_identifiability.md).

## A stubbed APSIM archetype

Real APSIM workflows are heavy, so this vignette uses a small simulator
stub with four parameters – two genuinely identifiable, one weakly
nonlinear, and one inert – and two outputs (yield and biomass).

``` r

library(kernR)

stub_apsim <- function(theta) {
  # theta: n x 4 matrix with columns slope, curvature, weak, inert
  yield   <- 1.5 * theta[, "slope"] +
             0.3 * stats::rnorm(nrow(theta), sd = 0.1)
  biomass <- theta[, "curvature"]^2 +
             0.4 * theta[, "weak"] +
             stats::rnorm(nrow(theta), sd = 0.1)
  cbind(yield = yield, biomass = biomass)
}
```

## Latin-hypercube design

``` r

bounds <- rbind(
  slope     = c(0.0, 2.0),   # active, linear in yield
  curvature = c(-1.0, 1.0),  # active, quadratic in biomass
  weak     = c(0.0, 1.0),    # weak linear effect on biomass
  inert    = c(0.0, 1.0)     # no effect
)
design <- lhs_design(n = 80L, bounds = bounds, seed = 11L)
head(design)
#>          slope   curvature      weak      inert
#> [1,] 0.8390559 -0.73444672 0.1898935 0.08378885
#> [2,] 1.3794481 -0.31872283 0.5788989 0.15216185
#> [3,] 0.6082601  0.01189946 0.2574184 0.22065620
#> [4,] 0.3924308 -0.86662113 0.6744266 0.35879365
#> [5,] 0.9017098  0.84070026 0.8906003 0.79415223
#> [6,] 1.4887851 -0.11886940 0.9224578 0.63657637
```

## Simulate and screen

``` r

outputs <- stub_apsim(design)
fit <- hsic_identifiability(
  theta          = design,
  y              = outputs,
  alpha          = 0.05,
  p_adjust       = "BH",
  n_permutations = 299L,
  seed           = 11L
)
fit
#> 
#>   HSIC Identifiability Scan
#> 
#> Parameters:   4 
#> Outputs:      2 
#> N:            80 
#> Permutations: 299 
#> Alpha:        0.05 
#> P-adjust:     BH 
#> 
#> Identifiable (3): slope, curvature, weak
#> Not identifiable (1): inert
#> 
#> Per-parameter ranking (descending max HSIC):
#>  parameter max_HSIC  min_p identifiable
#>      slope  0.09783 0.0089            *
#>  curvature  0.02458 0.0089            *
#>       weak  0.01057 0.0089            *
#>      inert 0.001757 0.9493             
#> 
#>   (* = identifiable at alpha = 0.05 )
```

The print method ranks parameters by their maximum HSIC across outputs.
`slope` and `curvature` should be flagged identifiable with small
p-values; `weak` may or may not survive depending on noise; `inert`
should fall below the threshold.

## Visual diagnostic

``` r

plot(fit)
```

![](kernR-identifiability_files/figure-html/plot-1.png)

## Handing off to PESTO

The identifiable subset is precisely the set you would forward to
[`pesto_ies()`](https://rdrr.io/pkg/PESTO/man/pesto_ies.html) (PESTO’s
IES wrapper) as the parameter prior. In pseudo-R:

    identifiable_params <- names(fit$identifiable)[fit$identifiable]
    pesto::pesto_ies(..., params = identifiable_params, ...)

For the cross-package workflow specification (kernR -\> PESTO -\> APSIM
ensemble loop), see the architecture document
`uq_ag_stack_design_v0.md`.

## Tabular form

For downstream reporting (or to feed into a manuscript table), call
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html):

``` r

df <- as.data.frame(fit)
head(df)
#>   parameter  output    statistic     p_value p_value_adjusted
#> 1     slope   yield 0.0978281488 0.003333333      0.008888889
#> 2 curvature   yield 0.0033398926 0.206666667      0.413333333
#> 3      weak   yield 0.0004566826 0.980000000      0.980000000
#> 4     inert   yield 0.0017572227 0.593333333      0.949333333
#> 5     slope biomass 0.0012291333 0.783333333      0.980000000
#> 6 curvature biomass 0.0245816130 0.003333333      0.008888889
```

## Notes on practice

- **Design size.** `n` should comfortably exceed the number of
  parameters; we recommend `n >= 10p` as a starting point for
  ag-systems-scale designs.
- **Permutations.** 199-499 permutations is usually sufficient; the
  minimum achievable p-value at `n_permutations` is
  `1 / (n_permutations + 1)`, so `n_permutations >= 199` is the floor
  for testing at `alpha = 0.05`.
- **Multiple testing.** With many parameters and several outputs,
  unadjusted p-values will produce spurious “identifiable” flags. The
  default Benjamini-Hochberg adjustment controls the false discovery
  rate across the `p x q` grid.
- **Cost.** Kernel matrices are reused across the grid, so cost scales
  as `O((p + q) n^2)` for kernel construction plus `O(p q B n^2)` for
  the permutation null, where `B = n_permutations`.

## References

- Hu, Z., Sejdinovic, D., & Evans, R. J. (2024). *A kernel-based
  statistical test for causal inference with backdoor adjustment.*
  Journal of Machine Learning Research, 25.
- Gretton, A., Fukumizu, K., Teo, C. H., Song, L., Schölkopf, B., &
  Smola, A. J. (2008). A kernel statistical test of independence.
  *NeurIPS*, 20.
- Da Veiga, S. (2015). Global sensitivity analysis with dependence
  measures. *Journal of Statistical Computation and Simulation*, 85(7),
  1283-1305.
