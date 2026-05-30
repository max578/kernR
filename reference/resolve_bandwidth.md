# Resolve Kernel Bandwidth

If `kernel$bandwidth` is `"median"`, compute the median heuristic from
the data. Otherwise return the fixed bandwidth.

## Usage

``` r
resolve_bandwidth(kernel, x)
```

## Arguments

- kernel:

  A `kernel_spec` object.

- x:

  Numeric matrix (n x d).

## Value

A `kernel_spec` with resolved numeric bandwidth.
