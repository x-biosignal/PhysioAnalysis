# Partial autocorrelation function (PACF)

Computes the sample partial autocorrelation function of a
one-dimensional signal at lags `0, 1, ..., lag_max`. The partial
autocorrelation at lag \\k\\ is the correlation between \\x_t\\ and
\\x\_{t+k}\\ with the linear effect of the intervening lags removed,
obtained here by the Durbin-Levinson recursion on the biased, demeaned
autocorrelation sequence (via
[`autocorrelation`](https://x-biosignal.github.io/PhysioAnalysis/reference/autocorrelation.md))
– the same definition as
[`stats::pacf`](https://rdrr.io/r/stats/acf.html) and
`statsmodels.tsa.pacf(method = "ldb")`. Where the ACF *decays*, the PACF
*cuts off* after the autoregressive order, so it is the standard tool
for choosing AR model order.

## Usage

``` r
partialAutocorrelation(x, lag_max = 20L)
```

## Arguments

- x:

  A numeric vector (the time series).

- lag_max:

  Maximum lag (default 20); the returned vector has `lag_max + 1` values
  (lags 0 through `lag_max`).

## Value

A numeric vector of partial autocorrelation coefficients at lags
`0:lag_max`; the lag-0 value is 1 by convention and the lag-1 value
equals the lag-1 autocorrelation.

## References

Box, G.E.P., Jenkins, G.M. (1976). Time Series Analysis: Forecasting and
Control. Holden-Day.

## See also

[`autocorrelation`](https://x-biosignal.github.io/PhysioAnalysis/reference/autocorrelation.md),
[`pacf`](https://rdrr.io/r/stats/acf.html)

## Examples

``` r
set.seed(1)
partialAutocorrelation(as.numeric(arima.sim(list(ar = 0.6), 500)), lag_max = 10)
#>  [1]  1.000000000  0.570267215 -0.014477934 -0.032401946  0.004478176
#>  [6]  0.013318710 -0.056459201  0.027871328 -0.030965381  0.071692262
#> [11] -0.093560384
```
