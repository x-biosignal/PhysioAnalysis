# Jarque-Bera test for normality

Applies the Jarque-Bera goodness-of-fit test for normality (Jarque &
Bera, 1980) to a one-dimensional signal: the statistic \\JB =
\frac{N}{6}\left(S^2 + \frac{K^2}{4}\right)\\, where \\S\\ is the
skewness and \\K\\ the excess kurtosis (from
[`signalMoments`](https://x-biosignal.github.io/PhysioAnalysis/reference/signalMoments.md)),
is compared to a chi-squared distribution with 2 degrees of freedom. It
is the formal, moment-based counterpart of the descriptive near-Gaussian
check: a large statistic rejects normality. Same definition as
`scipy.stats.jarque_bera`.

## Usage

``` r
jarqueBeraTest(x)
```

## Arguments

- x:

  A numeric vector (the time series).

## Value

A list with `statistic` (the JB statistic), `p_value` (the chi-squared
upper-tail probability), and `df` (2).

## References

Jarque, C.M., Bera, A.K. (1980). Efficient tests for normality,
homoscedasticity and serial independence of regression residuals.
*Economics Letters*, 6(3), 255-259.

## See also

[`signalMoments`](https://x-biosignal.github.io/PhysioAnalysis/reference/signalMoments.md)

## Examples

``` r
set.seed(1)
jarqueBeraTest(rnorm(1000))$p_value   # Gaussian -> large p (fail to reject)
#> [1] 0.9697854
```
