# Ljung-Box test for autocorrelation (white-noise / portmanteau test)

Applies the Ljung-Box portmanteau test (Ljung & Box, 1978) to a
one-dimensional signal: the statistic \\Q = N(N+2)\sum\_{k=1}^{h}
\hat\rho_k^2 / (N-k)\\ aggregates the first `lag` squared
autocorrelations (from
[`autocorrelation`](https://x-biosignal.github.io/PhysioAnalysis/reference/autocorrelation.md)),
and under the null hypothesis that the series is white noise \\Q\\
follows a chi-squared distribution with `lag` degrees of freedom. It is
the standard formal test for the presence of autocorrelation (and,
applied to model residuals, for model adequacy). Same definition as
`statsmodels.stats.diagnostic.acorr_ljungbox` and
`stats::Box.test(type = "Ljung-Box")`.

## Usage

``` r
ljungBoxTest(x, lag = 10L)
```

## Arguments

- x:

  A numeric vector (the time series).

- lag:

  Number of lags \\h\\ to include (default 10), and the degrees of
  freedom of the reference chi-squared.

## Value

A list with `statistic` (the Ljung-Box \\Q\\), `p_value` (the
chi-squared upper-tail probability), and `df` (`= lag`).

## References

Ljung, G.M., Box, G.E.P. (1978). On a measure of lack of fit in time
series models. *Biometrika*, 65(2), 297-303.

## See also

[`autocorrelation`](https://x-biosignal.github.io/PhysioAnalysis/reference/autocorrelation.md),
[`Box.test`](https://rdrr.io/r/stats/box.test.html)

## Examples

``` r
set.seed(1)
ljungBoxTest(rnorm(500), lag = 10)$p_value   # white noise -> large p
#> [1] 0.3980809
```
