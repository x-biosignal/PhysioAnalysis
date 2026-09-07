# Trimmed Mean (robust location)

Computes the trimmed mean of a one-dimensional signal – the arithmetic
mean after discarding the most extreme `trim` fraction of the values
from each end. A robust measure of central tendency that, unlike the
ordinary mean, is not pulled by a heavy tail or a few outliers, while
(unlike the median) it still uses all of the retained data. It matches
base R `mean(x, trim = ...)` and `scipy.stats.trim_mean` bit-for-bit,
both of which discard \\\lfloor n\\\mathrm{trim}\rfloor\\ values from
each end.

## Usage

``` r
trimmedMean(x, trim = 0.1)
```

## Arguments

- x:

  A numeric vector (the time series).

- trim:

  Fraction (0 to 0.5) of values trimmed from EACH end (default 0.1).

## Value

A single numeric value, the trimmed mean.

## Details

A robust location estimator that completes the robust toolkit alongside
[`medianAbsDev`](https://x-biosignal.github.io/PhysioAnalysis/reference/medianAbsDev.md)
and
[`interquartileRange`](https://x-biosignal.github.io/PhysioAnalysis/reference/interquartileRange.md)
(robust dispersion). On a skewed distribution the trimmed mean lies
between the ordinary mean and the median, and moves toward the median as
`trim` increases – a measure of how much the tail inflates the mean.

## See also

[`signalMoments`](https://x-biosignal.github.io/PhysioAnalysis/reference/signalMoments.md)
for the ordinary mean and SD,
[`medianAbsDev`](https://x-biosignal.github.io/PhysioAnalysis/reference/medianAbsDev.md)
and
[`interquartileRange`](https://x-biosignal.github.io/PhysioAnalysis/reference/interquartileRange.md)
for robust dispersion.

## Examples

``` r
trimmedMean(c(1, 2, 3, 4, 100), trim = 0.2)  # 3 (robust); the mean is 22
#> [1] 3
set.seed(1); trimmedMean(rexp(1000), trim = 0.2)  # below the tail-inflated mean
#> [1] 0.8067975
```
