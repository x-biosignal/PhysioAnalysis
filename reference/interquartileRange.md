# Interquartile Range (robust spread)

Computes the interquartile range (IQR) of a one-dimensional signal – the
spread of the middle 50%, \\Q_3 - Q_1\\ – a robust, quantile-based
measure of dispersion that ignores the tails entirely. The quartiles use
R's default type-7 quantile (the same linear interpolation as
`scipy.stats.iqr` and `numpy.percentile`), so the result matches both
bit-for-bit.

## Usage

``` r
interquartileRange(x, type = 7L)
```

## Arguments

- x:

  A numeric vector (the time series).

- type:

  Quantile algorithm passed to
  [`quantile`](https://rdrr.io/r/stats/quantile.html) (default `7`,
  matching base R, scipy and numpy).

## Value

A single numeric value, the interquartile range.

## Details

A quantile-based robust spread that complements
[`medianAbsDev`](https://x-biosignal.github.io/PhysioAnalysis/reference/medianAbsDev.md)
(a deviation-based robust scale) and
[`signalMoments`](https://x-biosignal.github.io/PhysioAnalysis/reference/signalMoments.md)
(the classical SD). For approximately Gaussian data, \\\mathrm{IQR} /
(2\\\Phi^{-1}(0.75))\\ (\\\approx \mathrm{IQR}/1.349\\) is a robust
estimate of the SD; a value well below the SD flags a skewed or
heavy-tailed distribution.

## See also

[`medianAbsDev`](https://x-biosignal.github.io/PhysioAnalysis/reference/medianAbsDev.md)
for the deviation-based robust scale,
[`signalMoments`](https://x-biosignal.github.io/PhysioAnalysis/reference/signalMoments.md)
for the classical moments.

## Examples

``` r
interquartileRange(c(1, 2, 3, 4, 100))   # robust to the outlier
#> [1] 2
interquartileRange(rnorm(1000)) / 1.349  # ~1 (robust SD estimate)
#> [1] 1.006017
```
