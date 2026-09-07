# KPSS test for stationarity

Computes the KPSS statistic (Kwiatkowski, Phillips, Schmidt & Shin,
1992) for level stationarity of a one-dimensional signal. The residuals
from regressing the signal on a constant are cumulatively summed, and
the statistic is \\\eta = N^{-2}\sum_t S_t^2 / \hat\lambda^2\\, where
\\S_t\\ are the partial sums and \\\hat\lambda^2\\ is a Bartlett
(Newey-West) long-run variance estimate with truncation `lag`. Unlike
the Augmented Dickey-Fuller test, the KPSS *null* is stationarity: a
*small* statistic (below the critical value) means the series is
consistent with stationarity. The statistic is the same deterministic
quantity as
`statsmodels.tsa.stattools.kpss(regression = "c", nlags = lag)`.

## Usage

``` r
kpssTest(x, lag = 10L)
```

## Arguments

- x:

  A numeric vector (the time series).

- lag:

  Bartlett truncation lag \\L\\ for the long-run variance (default 10).

## Value

A list with `statistic` (the KPSS \\\eta\\), `n` (the signal length),
and `lag`.

## Details

Note: only the *statistic* is computed here (an exact quantity). The
KPSS critical values and p-value come from the reference tables (an
approximation) and are not returned.

## References

Kwiatkowski, D., Phillips, P.C.B., Schmidt, P., Shin, Y. (1992). Testing
the null hypothesis of stationarity against the alternative of a unit
root. *Journal of Econometrics*, 54(1-3), 159-178.

## See also

[`adfTest`](https://x-biosignal.github.io/PhysioAnalysis/reference/adfTest.md),
[`autocorrelation`](https://x-biosignal.github.io/PhysioAnalysis/reference/autocorrelation.md)

## Examples

``` r
set.seed(1)
kpssTest(rnorm(500), lag = 10)$statistic   # stationary -> small statistic
#> [1] 0.07720669
```
