# Augmented Dickey-Fuller (ADF) unit-root test statistic

Computes the Augmented Dickey-Fuller test statistic for a unit root
(non-stationarity) in a one-dimensional signal. The test regresses the
first difference on a constant, the lagged level, and `lag` lagged first
differences, \\\Delta y_t = \alpha + \beta y\_{t-1} + \sum\_{i=1}^{L}
\gamma_i \Delta y\_{t-i} + \varepsilon_t\\, and returns the t-statistic
on \\\beta\\. A large negative statistic (below the Dickey-Fuller
critical value) rejects the unit-root null in favour of stationarity.
The statistic is the same deterministic OLS quantity as
`statsmodels.tsa.stattools.adfuller(..., autolag = None)`.

## Usage

``` r
adfTest(x, lag = 1L)
```

## Arguments

- x:

  A numeric vector (the time series).

- lag:

  Number of augmenting lagged differences \\L\\ (default 1).

## Value

A list with `statistic` (the ADF t-statistic on the lagged level), `n`
(the regression sample size), and `lag`.

## Details

Note: only the *statistic* is computed here (it is an exact regression
quantity). The critical values and p-value require the Dickey-Fuller /
MacKinnon reference tables (an approximation), so they are not returned;
compare the statistic against the tabulated critical value for the
decision.

## References

Dickey, D.A., Fuller, W.A. (1979). Distribution of the estimators for
autoregressive time series with a unit root. *JASA*, 74, 427-431.

## See also

[`autocorrelation`](https://x-biosignal.github.io/PhysioAnalysis/reference/autocorrelation.md),
[`ljungBoxTest`](https://x-biosignal.github.io/PhysioAnalysis/reference/ljungBoxTest.md)

## Examples

``` r
set.seed(1)
adfTest(cumsum(rnorm(500)), lag = 1)$statistic   # random walk -> not far below 0
#> [1] -2.866323
```
