# Autoregressive (parametric) spectral density

Computes the parametric AR spectral density of a one-dimensional signal
from an autoregressive model of order `order` fitted by Yule-Walker (via
[`arYuleWalker`](https://x-biosignal.github.io/PhysioAnalysis/reference/arYuleWalker.md)):
\\S(f) = \sigma^2 / \|1 - \sum_k \phi_k e^{-i 2\pi f k}\|^2\\, evaluated
on `n_freq` normalized frequencies in \\\[0, 0.5\]\\. Unlike the
periodogram / Welch estimate, the AR spectrum is smooth and
low-variance, with sharp peaks at resonant frequencies, using only
`order` parameters. It matches
[`stats::spec.ar`](https://rdrr.io/r/stats/spec.ar.html) (including its
small-sample innovation-variance convention \\N/(N-p-1)\\).

## Usage

``` r
arSpectrum(x, order = 8L, n_freq = 500L)
```

## Arguments

- x:

  A numeric vector (the time series).

- order:

  Autoregressive order \\p\\ (default 8).

- n_freq:

  Number of frequencies in \\\[0, 0.5\]\\ (default 500).

## Value

A list with `freq` (normalized frequencies, cycles/sample), `spec` (the
AR spectral density), and `order`.

## References

Percival, D.B., Walden, A.T. (1993). Spectral Analysis for Physical
Applications. Cambridge University Press.

## See also

[`arYuleWalker`](https://x-biosignal.github.io/PhysioAnalysis/reference/arYuleWalker.md),
[`spec.ar`](https://rdrr.io/r/stats/spec.ar.html)

## Examples

``` r
set.seed(1)
s <- arSpectrum(as.numeric(arima.sim(list(ar = c(0.6, -0.3)), 500)), order = 2)
s$freq[which.max(s$spec)]   # peak (normalized) frequency
#> [1] 0.1452906
```
