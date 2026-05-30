# Compute the unbiased MMD^2 statistic

Unbiased estimator of MMD^2: \$\$\frac{1}{n(n-1)} \sum\_{i \ne j}
K\_{xx}(i,j) + \frac{1}{m(m-1)} \sum\_{i \ne j} K\_{yy}(i,j) -
\frac{2}{nm} \sum\_{i,j} K\_{xy}(i,j)\$\$

## Usage

``` r
mmd2_unbiased_cpp(Kxx, Kyy, Kxy)
```

## Arguments

- Kxx:

  n x n kernel matrix for sample X.

- Kyy:

  m x m kernel matrix for sample Y.

- Kxy:

  n x m kernel matrix between X and Y.

## Value

Scalar unbiased MMD^2 value.
