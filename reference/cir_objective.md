# Objective Causal Influence Range (threshold-free)

Reduces a divergence-versus-lag profile to a single, threshold-free
causal influence range (CIR) – the effective horizon over which a cause
keeps informing the estimate of its effect (Andreou et al. 2026, eqs.
8-9). Given the relative entropy `divergence` between the **lagged**
smoother (future of the effect included only up to a lag `L`) and the
**complete** smoother (all future included), the divergence falls from
\\M\\ at lag 0 (the filter, no future) towards 0 as `L` grows. The
subjective CIR at tolerance \\\varepsilon\\ is \\\tau\_\varepsilon =
\inf\\L : D(L) \le \varepsilon\\\\; integrating it out, \$\$\tau =
\frac{1}{M}\int_0^M \tau\_\varepsilon \\ d\varepsilon =
\frac{1}{M}\int_0^{L\_{\max}} D(L)\\ dL,\$\$ the second equality by
parts. The result is a decorrelation-time analogue: a lead time in the
units of `lag`, free of any cut-off threshold.

## Usage

``` r
cir_objective(lag, divergence)
```

## Arguments

- lag:

  Numeric vector of non-negative, increasing lags `L` (units of time).
  Should start at (or include) 0.

- divergence:

  Numeric vector, same length as `lag`: the relative entropy \\D(L)\\ of
  the lag-`L` smoother from the complete smoother. The value at the
  smallest lag is taken as \\M\\ (the normalising maximum).

## Value

A single non-negative numeric: the objective CIR (a lead time in the
units of `lag`). Zero when there is no recoverable future information
(\\M = 0\\).

## Details

The producer of the `divergence` profile (the lagged/complete smoother
passes) is the upstream filtering engine; this function owns only the
range integral.

## References

Andreou, M., Chen, N. & Bollt, E. (2026). Assimilative causal inference.
*Nature Communications* 17, 1854.

## Examples

``` r
lag <- seq(0, 10, by = 0.5)
D <- 2 * exp(-0.6 * lag)            # divergence decays with lag
cir_objective(lag, D)
#> [1] 1.674986
```
