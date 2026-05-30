# Validate and Coerce Input Data

Converts input to a numeric matrix with validation.

## Usage

``` r
validate_input(x, name = "x", min_n = 1L, min_d = 1L)
```

## Arguments

- x:

  Input data (vector, matrix, data.frame, data.table).

- name:

  Name of the argument (for error messages).

- min_n:

  Minimum number of observations.

- min_d:

  Minimum number of columns.

## Value

A numeric matrix.
