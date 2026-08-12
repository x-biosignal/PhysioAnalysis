# Cliff's delta for two independent samples

Cliff's delta, `delta = P(X > Y) - P(X < Y)`, a non-parametric effect
size in `[-1, 1]` whose sign follows the direction of the difference.
The confidence interval uses Cliff's (1993) consistent variance with his
asymmetric in-bounds transform.

## Usage

``` r
cliffsDelta(x, y, conf_level = 0.95)
```

## Arguments

- x, y:

  Numeric vectors for the two samples.

- conf_level:

  Confidence level (default 0.95).

## Value

A list with `delta`, `ci_lower`, `ci_upper`, the variance estimate and
the sample sizes.

## References

Cliff, N. (1993). Psychol Bull 114(3):494-509.

## See also

[`rankBiserial()`](https://x-biosignal.github.io/PhysioAnalysis/reference/rankBiserial.md),
[`effectSize()`](https://x-biosignal.github.io/PhysioAnalysis/reference/effectSize.md)

## Examples

``` r
cliffsDelta(c(5, 6, 7, 8), c(1, 2, 3, 9))
#> $delta
#> [1] 0.5
#> 
#> $ci_lower
#> [1] -0.4878146
#> 
#> $ci_upper
#> [1] 0.9263176
#> 
#> $variance
#> [1] 0.25
#> 
#> $n1
#> [1] 4
#> 
#> $n2
#> [1] 4
#> 
#> $conf_level
#> [1] 0.95
#> 
```
