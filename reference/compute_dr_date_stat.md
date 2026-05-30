# Compute DR-DATE Test Statistic (augmented IPW)

Compute DR-DATE Test Statistic (augmented IPW)

## Usage

``` r
compute_dr_date_stat(Ky, treatment, e_hat, C1, C0)
```

## Arguments

- Ky:

  n x n outcome kernel matrix.

- treatment:

  Binary treatment vector.

- e_hat:

  Propensity scores.

- C1, C0:

  n x n conditional-mean-embedding coefficient matrices for the treated
  and control arms. Row `i` holds the coefficients of \\\hat m_a(x_i)\\
  over the n outcome embeddings; zero matrices give the
  inverse-probability-weighted statistic.

## Value

Scalar test statistic.
