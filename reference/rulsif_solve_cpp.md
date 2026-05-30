# RuLSIF kernel density ratio estimation (core solver)

Solves the RuLSIF optimisation \$\$\theta = (H + \lambda I)^{-1} h\$\$
with non-negativity constraint (clamp negatives to 0).

## Usage

``` r
rulsif_solve_cpp(H, h, lambda)
```

## Arguments

- H:

  Gram matrix (n_basis x n_basis).

- h:

  Mean kernel vector (n_basis).

- lambda:

  Regularisation parameter.

## Value

Coefficient vector theta.
