# Autocorrelation function (ACF)

Computes the sample autocorrelation function of a one-dimensional signal
at lags `0, 1, ..., lag_max` using the biased estimator (dividing each
autocovariance by \\N\\), with the mean removed by default – the same
definition as [`stats::acf()`](https://rdrr.io/r/stats/acf.html) and
`statsmodels.tsa.acf(adjusted = FALSE)`. The autocorrelation measures
how similar a signal is to a time-shifted copy of itself, and is the
basic tool for detecting periodicity, rhythmicity, and the decorrelation
time of a physiological signal.

## Usage

``` r
autocorrelation(x, lag_max = 30L, demean = TRUE)
```

## Arguments

- x:

  A numeric vector (the time series).

- lag_max:

  Maximum lag (default 30); the returned vector has `lag_max + 1` values
  (lags 0 through `lag_max`).

- demean:

  If `TRUE` (default), subtract the mean before computing the
  autocovariance (the standard definition).

## Value

A numeric vector of autocorrelation coefficients at lags `0:lag_max`;
the lag-0 value is always 1.

## References

Box, G.E.P., Jenkins, G.M. (1976). Time Series Analysis: Forecasting and
Control. Holden-Day.

## See also

[`acf`](https://rdrr.io/r/stats/acf.html)

## Examples

``` r
set.seed(1)
autocorrelation(sin(seq(0, 20 * pi, length.out = 500)), lag_max = 20)
#>  [1]  1.00000000  0.99208311  0.96852099  0.92974895  0.87644110  0.80949880
#>  [7]  0.73003544  0.63935797  0.53894543  0.43042487  0.31554498  0.19614791
#> [13]  0.07413976 -0.04853987 -0.16994869 -0.28817248 -0.40135531 -0.50772859
#> [19] -0.60563878 -0.69357304 -0.77018262
```
