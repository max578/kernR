# Calibration and concordance: kernR's validation layer

`kernR` is a verdict layer: it does not fit models, it decides whether a
fitted model’s output should be trusted. Two questions sit at the centre
of that role, and this vignette covers the two tests that answer them.

- **Calibration.** Does a sample actually follow the distribution it
  claims to represent? This is a one-sample, score-based question,
  answered by
  [`ksd_test()`](https://max578.github.io/kernR/reference/ksd_test.md) –
  the kernel Stein discrepancy goodness-of-fit test.
- **Concordance.** Do several sources – posterior draws from different
  inference engines, or scenario ensembles from different simulators –
  agree with each other? This is a k-sample, sample-based question,
  answered by
  [`concordance_test()`](https://max578.github.io/kernR/reference/concordance_test.md).

The two are complementary. Calibration compares a sample against a
*density*; concordance compares samples against *each other*. Together
they let `kernR` act as a falsification gate over an upstream inference
pipeline.

## Calibration: `ksd_test()`

The kernel Stein discrepancy needs only the *score* of the target – the
gradient of its log density – so the target may be unnormalised and no
reference sample is required. For a multivariate-normal target the score
is available in closed form via
[`gaussian_score()`](https://max578.github.io/kernR/reference/gaussian_score.md).

``` r

set.seed(1)

# A sample that genuinely follows the standard normal
x_ok <- matrix(rnorm(400), ncol = 2)
ksd_test(x_ok, n_boot = 199, seed = 1)
#> 
#>    KSD Test
#> 
#> Statistic: -0.00225082 
#> P-value:   0.6300 
#> N:         200 
#> Perms:     199 
#> Kernel X:  imq
#> 
#> Goodness-of-fit verdict
#>   Stein kernel: imq (beta = -0.5)
#>   Bandwidth:    1.608 (median heuristic)
#>   Bootstrap:    wild, B = 199
#>   Surprise:     0.667 bits
#>   Verdict:      consistent with target
```

The verdict is *consistent with target*: the statistic sits near zero
and the p-value is large. A mis-specified sample – here shifted in mean
– is caught.

``` r

x_bad <- x_ok + 0.7
ksd_test(x_bad, n_boot = 199, seed = 1)
#> 
#>    KSD Test
#> 
#> Statistic: 0.477917 
#> P-value:   0.0050 
#> N:         200 
#> Perms:     199 
#> Kernel X:  imq
#> 
#> Goodness-of-fit verdict
#>   Stein kernel: imq (beta = -0.5)
#>   Bandwidth:    1.608 (median heuristic)
#>   Bootstrap:    wild, B = 199
#>   Surprise:     7.644 bits
#>   Verdict:      REJECT (sample inconsistent with target)
```

When the target has no convenient closed-form score, supply it as a log
density and let
[`numeric_score()`](https://max578.github.io/kernR/reference/numeric_score.md)
take the gradient by finite differences. Any additive normalising
constant cancels, so an unnormalised log density is enough.

``` r

log_density <- function(z) -0.5 * rowSums(z^2)   # standard normal, unnormalised
ksd_test(x_ok, score = numeric_score(log_density), n_boot = 199, seed = 1)
#> 
#>    KSD Test
#> 
#> Statistic: -0.00225082 
#> P-value:   0.6300 
#> N:         200 
#> Perms:     199 
#> Kernel X:  imq
#> 
#> Goodness-of-fit verdict
#>   Stein kernel: imq (beta = -0.5)
#>   Bandwidth:    1.608 (median heuristic)
#>   Bootstrap:    wild, B = 199
#>   Surprise:     0.667 bits
#>   Verdict:      consistent with target
```

This is the adapter that lets
[`ksd_test()`](https://max578.github.io/kernR/reference/ksd_test.md)
consume any externally supplied log-posterior evaluator: wrap the
evaluator in
[`numeric_score()`](https://max578.github.io/kernR/reference/numeric_score.md)
and the test asks whether posterior draws are calibrated against that
posterior, with no dependency on the producer of the evaluator.

The default base kernel is the inverse multi-quadric, which stays
sensitive to mis-specification in higher dimensions where the Gaussian
kernel can lose power; pass `kernel = "rbf"` for the Gaussian
alternative.

## Concordance: `concordance_test()`

Where calibration checks one sample against a density, concordance
checks several samples against one another. The input is a list – one
element per source – and a named list labels the sources in the output.

``` r

set.seed(2)
draws <- list(
  engine_a = matrix(rnorm(400), ncol = 2),
  engine_b = matrix(rnorm(400), ncol = 2),
  engine_c = matrix(rnorm(400), ncol = 2)
)
concordance_test(draws, n_permutations = 199, seed = 1)
#> 
#>    Concordance Test
#> 
#> Statistic: 0.00131634 
#> P-value:   0.3650 
#> N:         600 
#> Perms:     199 
#> Kernel X:  rbf (bw = 1.682)
#> 
#> Concordance verdict
#>   Sources:    3 (engine_a, engine_b, engine_c)
#>   Verdict:    concordant
#>   Pairwise MMD-squared:
#>          engine_a engine_b  engine_c 
#> engine_a    0     0.00165   0.000591 
#> engine_b 0.00165     0      -0.000929
#> engine_c 0.000591 -0.000929    0
```

The three sources are mutually concordant, so the verdict is
*concordant*. The test is more than a single yes-or-no, though: it
returns the full pairwise discrepancy matrix, so when one source departs
the rejection can be read down to the offending pair.

``` r

draws$engine_c <- draws$engine_c + 1   # engine_c now disagrees
fit <- concordance_test(draws, n_permutations = 199, seed = 1)
fit$pairwise
#>             engine_a    engine_b  engine_c
#> engine_a 0.000000000 0.001600717 0.1787739
#> engine_b 0.001600717 0.000000000 0.2067248
#> engine_c 0.178773878 0.206724754 0.0000000
```

The pairwise matrix localises the problem: the `engine_c` row and column
carry the large discrepancies, while `engine_a` and `engine_b` remain
close. The single joint-permutation null keeps the overall p-value
calibrated across all pairs, so this is one test with one verdict, not a
multiple-comparison sweep.

## The validation pattern

Used together the two tests express a simple discipline. Concordance
across independent sources is corroborating evidence – agreement is hard
to fake. Calibration against a target density is the absolute check that
the agreed-upon answer is also the correct one. Divergence on either
test is informative:
[`concordance_test()`](https://max578.github.io/kernR/reference/concordance_test.md)
localises *which* source departs, and
[`ksd_test()`](https://max578.github.io/kernR/reference/ksd_test.md)
identifies *whether* a sample departs from its claimed target. A
pipeline whose ensembles pass both has earned more trust than one
validated on the posterior mean alone.
