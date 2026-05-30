# Parse Causal Formula

Parses `y ~ treatment | confounders` into component matrices.

## Usage

``` r
parse_causal_formula(formula, data)
```

## Arguments

- formula:

  Formula.

- data:

  data.frame.

## Value

List with `y`, `treatment`, `covariates`.
