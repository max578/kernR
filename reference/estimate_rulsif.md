# Estimate Density Ratios via RuLSIF

Relative unconstrained Least-Squares Importance Fitting. Kernel-based
closed-form density ratio estimation.

## Usage

``` r
estimate_rulsif(x_num, x_den, kernel = kernel_spec(), lambda = 0.1, alpha = 0)
```

## Arguments

- x_num:

  Numeric matrix. Numerator samples.

- x_den:

  Numeric matrix. Denominator samples.

- kernel:

  Kernel specification. Default is RBF with median heuristic.

- lambda:

  Regularisation parameter. Default is 0.1.

- alpha:

  Relative parameter (0 = LSIF, 0.5 = symmetric). Default is 0.

## Value

Named list with `weights` and `ess`.
