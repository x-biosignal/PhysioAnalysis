# Median Absolute Deviation (robust dispersion)

Computes the median absolute deviation (MAD) of a one-dimensional signal
– \\\mathrm{MAD} = c \cdot \mathrm{median}(\|x - \mathrm{median}(x)\|)\\
– a robust measure of dispersion that, unlike the standard deviation, is
not inflated by a small number of extreme values (spikes, movement
artifacts). With the default `constant = 1` it returns the raw MAD; set
`constant = 1 / stats::qnorm(0.75)` (approximately 1.4826) to obtain the
normal-consistent robust estimate of the standard deviation (the scale
MAD and SD share for Gaussian data).

## Usage

``` r
medianAbsDev(x, constant = 1)
```

## Arguments

- x:

  A numeric vector (the time series).

- constant:

  Scale factor (default `1`, the raw MAD). Use `1 / stats::qnorm(0.75)`
  for the normal-consistent robust SD estimate.

## Value

A single numeric value, the (scaled) median absolute deviation.

## Details

The robust dispersion companion of
[`signalMoments`](https://x-biosignal.github.io/PhysioAnalysis/reference/signalMoments.md)
(whose SD is the classical, non-robust dispersion). Comparing the
normal-consistent MAD with the SD is a quick artifact check: they agree
for clean, near-Gaussian data and diverge (MAD-scale below SD) when
heavy-tailed artifacts inflate the SD.

## See also

[`signalMoments`](https://x-biosignal.github.io/PhysioAnalysis/reference/signalMoments.md)
for the classical moments (mean, SD, skewness, kurtosis).

## Examples

``` r
medianAbsDev(c(1, 2, 3, 4, 100))                         # robust to the outlier
#> [1] 1
set.seed(1); medianAbsDev(rnorm(1000), 1 / stats::qnorm(0.75))  # ~1 (robust SD)
#> [1] 1.031688
```
