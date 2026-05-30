# Plot HSIC-Sensitivity Indices

Bar plot of per-parameter HSIC-Sensitivity Index, ordered by first-order
index magnitude. When `total_order` was set on the fit,
`which = "total"` shows the total-order indices and `which = "both"`
shows side-by-side bars for `S` and `T`.

## Usage

``` r
# S3 method for class 'hsic_sensitivity'
plot(
  x,
  which = c("first", "total", "both"),
  alpha = 0.05,
  col_sig = "#0072B2",
  col_nonsig = "grey70",
  col_total = "#D55E00",
  ...
)
```

## Arguments

- x:

  An `hsic_sensitivity` object.

- which:

  Character. `"first"` (default), `"total"`, or `"both"`. `"total"` /
  `"both"` require `x$total_order` to be `TRUE`.

- alpha:

  Numeric in `(0, 1)`. Significance level for colour coding (first-order
  only; ignored if the object carries no p-values or under
  `which = "total"`).

- col_sig, col_nonsig, col_total:

  Bar colours.

- ...:

  Additional arguments passed to
  [`graphics::barplot()`](https://rdrr.io/r/graphics/barplot.html).

## Value

Invisibly returns `x`. Side effect: produces a base R plot.
