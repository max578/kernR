# Relative Entropy Between Two Gaussians (the ACI causal metric)

Computes the Kullback-Leibler divergence \\\mathcal{P}(p, q) = \int p
\log(p / q)\\ between two multivariate Gaussian distributions \\p =
\mathcal{N}(\mu_p, \Sigma_p)\\ and \\q = \mathcal{N}(\mu_q, \Sigma_q)\\.
This is the operational statistic of Assimilative Causal Inference
(ACI): with `p` the **smoother** posterior of a hidden cause (using the
future of the observed effect) and `q` the **filter** posterior (past
only), a non-zero value identifies the hidden variable as a cause of the
observed effect at that instant.

## Usage

``` r
relative_entropy(mu_p, sigma_p, mu_q, sigma_q)
```

## Arguments

- mu_p:

  Numeric vector. Mean of `p` (the smoother posterior, for ACI).

- sigma_p:

  Covariance of `p`: matrix, variance vector, or scalar.

- mu_q:

  Numeric vector. Mean of `q` (the filter posterior, for ACI).

- sigma_q:

  Covariance of `q`: matrix, variance vector, or scalar. Must be
  positive definite (it is inverted).

## Value

A single non-negative numeric: the relative entropy \\\mathcal{P}(p,
q)\\.

## Details

The closed form for `k`-dimensional Gaussians is
\$\$\tfrac{1}{2}\left\[\operatorname{tr}(\Sigma_q^{-1}\Sigma_p) +
(\mu_q-\mu_p)^\top \Sigma_q^{-1} (\mu_q-\mu_p) - k + \log\frac{\det
\Sigma_q}{\det \Sigma_p}\right\].\$\$ The measure is non-negative, zero
if and only if the two Gaussians coincide, and asymmetric in `p` and `q`
(it is a divergence, not a distance). It captures differences in both
mean and covariance and is invariant under any smooth invertible
reparameterisation of the state.

Covariance arguments accept a `k x k` matrix, a length-`k` variance
vector (diagonal covariance), or a scalar (isotropic). The reference
covariance \\\Sigma_q\\ must be positive definite.

## References

Andreou, M., Chen, N. & Bollt, E. (2026). Assimilative causal inference.
*Nature Communications* 17, 1854.

## Examples

``` r
relative_entropy(0, 1, 0, 1)            # identical -> 0
#> [1] 0
relative_entropy(1, 1, 0, 1)            # mean shift of 1 sd -> 0.5
#> [1] 0.5
relative_entropy(c(0, 0), diag(2), c(1, 0), diag(2))
#> [1] 0.5
```
