# Polynomial kernel matrix

K(x, y) = (x %\*% t(y) + offset)^degree

## Usage

``` r
polynomial_kernel_matrix_cpp(x, y, degree, offset)
```

## Arguments

- x:

  Numeric matrix (n x d).

- y:

  Numeric matrix (m x d).

- degree:

  Integer degree.

- offset:

  Scalar offset.

## Value

n x m kernel matrix.
