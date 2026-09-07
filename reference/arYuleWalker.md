# Autoregressive model fit by Yule-Walker (AR(p))

Fits an autoregressive model of order `order` to a one-dimensional
signal by the Yule-Walker method: the AR coefficients solve the
Yule-Walker equations formed from the signal's biased demeaned
autocovariance (obtained via the Durbin-Levinson recursion on
[`autocorrelation`](https://x-biosignal.github.io/PhysioAnalysis/reference/autocorrelation.md)),
and the innovation variance is the residual (one-step prediction)
variance. This is the same definition as
`statsmodels.regression.linear_model.yule_walker(method = "mle")`. AR
modelling underlies forecasting, parametric (AR) spectra, and – with the
order chosen from the
[`partialAutocorrelation`](https://x-biosignal.github.io/PhysioAnalysis/reference/partialAutocorrelation.md)
cut-off – compact descriptions of oscillatory dynamics.

## Usage

``` r
arYuleWalker(x, order = 4L)
```

## Arguments

- x:

  A numeric vector (the time series).

- order:

  Autoregressive order \\p\\ (default 4).

## Value

A list with `ar` (the length-`order` AR coefficient vector \\\phi_1,
\dots, \phi_p\\), `var_pred` (the innovation / one-step prediction
variance), and `order`.

## References

Box, G.E.P., Jenkins, G.M. (1976). Time Series Analysis: Forecasting and
Control. Holden-Day.

## See also

[`autocorrelation`](https://x-biosignal.github.io/PhysioAnalysis/reference/autocorrelation.md),
[`partialAutocorrelation`](https://x-biosignal.github.io/PhysioAnalysis/reference/partialAutocorrelation.md),
[`ar.yw`](https://rdrr.io/r/stats/ar.html)

## Examples

``` r
set.seed(1)
arYuleWalker(as.numeric(arima.sim(list(ar = c(0.5, -0.3)), 500)), order = 2)
#> $ar
#> [1]  0.4890456 -0.3153433
#> 
#> $var_pred
#> [1] 1.042764
#> 
#> $order
#> [1] 2
#> 
```
