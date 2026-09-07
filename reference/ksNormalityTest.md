# Kolmogorov-Smirnov Test for Normality

Tests whether a one-dimensional signal is normally distributed by
comparing its empirical cumulative distribution function (ECDF) to that
of a normal distribution fitted to the sample (mean and SD). The
statistic is the Kolmogorov-Smirnov distance \\D = \sup_x \|F_n(x) -
\Phi((x-\mu)/\sigma)\|\\, the largest vertical gap between the two CDFs;
the two-sided asymptotic (Kolmogorov) p-value follows. The \\D\\
statistic reproduces
[`stats::ks.test`](https://rdrr.io/r/stats/ks.test.html) and
`scipy.stats.kstest` bit-for-bit.

## Usage

``` r
ksNormalityTest(x)
```

## Arguments

- x:

  A numeric vector (the time series).

## Value

A list with `statistic` (the KS distance \\D\\), `p_value` (two-sided
asymptotic), and `n` (sample size).

## Details

A distribution-shape (empirical-CDF) normality test, complementary to
the moment-based
[`jarqueBeraTest`](https://x-biosignal.github.io/PhysioAnalysis/reference/jarqueBeraTest.md):
KS responds to the overall shape of the distribution, Jarque-Bera to its
skewness and kurtosis (the tails). The two can disagree – on a long,
mildly non-Gaussian signal Jarque-Bera's power often rejects while KS
does not.

## Caveat

Because the normal parameters are estimated from the same sample, the
naive p-value is anti-conservative (a Lilliefors correction gives an
exact test); the \\D\\ statistic itself is exact and is what this
function certifies.

## See also

[`jarqueBeraTest`](https://x-biosignal.github.io/PhysioAnalysis/reference/jarqueBeraTest.md)
for the moment-based normality test,
[`signalMoments`](https://x-biosignal.github.io/PhysioAnalysis/reference/signalMoments.md)
for the skewness and kurtosis.

## Examples

``` r
set.seed(1); ksNormalityTest(rnorm(500))$p_value    # Gaussian -> large p
#> [1] 0.5003172
set.seed(1); ksNormalityTest(rexp(500))$statistic   # skewed -> large D
#> [1] 0.1432139
```
