# Signal distribution moments (mean, SD, skewness, kurtosis)

Computes the first four moments of a one-dimensional signal's amplitude
distribution: the mean, the (population) standard deviation, the
skewness, and the excess kurtosis. Skewness measures asymmetry (negative
= a longer left tail); excess kurtosis measures tailedness relative to a
Gaussian (0 = Gaussian, positive = heavier tails / more peaked). These
are the standard shape descriptors used, for example, to flag EEG
artifacts (high kurtosis) or characterize amplitude asymmetry.
Definitions match `scipy.stats.skew` (biased) and `scipy.stats.kurtosis`
(Fisher / excess, biased), with the population standard deviation
(`ddof = 0`).

## Usage

``` r
signalMoments(x)
```

## Arguments

- x:

  A numeric vector (the time series).

## Value

A list with `mean`, `sd` (population), `skewness`, and `kurtosis`
(excess).

## References

Joanes, D.N., Gill, C.A. (1998). Comparing measures of sample skewness
and kurtosis. *The Statistician*, 47(1), 183-189.

## See also

[`sd`](https://rdrr.io/r/stats/sd.html)

## Examples

``` r
set.seed(1)
signalMoments(rnorm(1000))   # ~ 0 skewness, ~ 0 excess kurtosis for Gaussian
#> $mean
#> [1] -0.01164814
#> 
#> $sd
#> [1] 1.034398
#> 
#> $skewness
#> [1] -0.0191671
#> 
#> $kurtosis
#> [1] -0.001775465
#> 
```
