# Compute the biased MMD^2 statistic

Compute the biased MMD^2 statistic

## Usage

``` r
mmd2_biased_cpp(Kxx, Kyy, Kxy)
```

## Arguments

- Kxx:

  n x n kernel matrix for sample X.

- Kyy:

  m x m kernel matrix for sample Y.

- Kxy:

  n x m kernel matrix between X and Y.

## Value

Scalar biased MMD^2 value.
