# Find significant time windows

Identifies contiguous time periods with significant effects.

## Usage

``` r
findSignificantWindows(p_values, times = NULL, alpha = 0.05, min_duration = 0)
```

## Arguments

- p_values:

  Vector of p-values across time.

- times:

  Vector of time points.

- alpha:

  Significance threshold (default: 0.05).

- min_duration:

  Minimum duration of significant window in time units.

## Value

A data.frame with columns: start, end, duration, min_p.

## References

Maris, E. & Oostenveld, R. (2007). "Nonparametric statistical testing of
EEG- and MEG-data." *Journal of Neuroscience Methods*, 164(1), 177-190.
[doi:10.1016/j.jneumeth.2007.03.024](https://doi.org/10.1016/j.jneumeth.2007.03.024)

## See also

[`correctPValues()`](https://x-biosignal.github.io/PhysioAnalysis/reference/correctPValues.md)
for multiple comparison correction,
[`clusterPermutationTest()`](https://x-biosignal.github.io/PhysioAnalysis/reference/clusterPermutationTest.md)
for cluster-based permutation testing,
[`tTestEpochs()`](https://x-biosignal.github.io/PhysioAnalysis/reference/tTestEpochs.md)
for pointwise t-tests.

## Examples

``` r
# Example p-values
times <- seq(-0.2, 0.8, length.out = 100)
p <- c(rep(0.5, 30), rep(0.01, 20), rep(0.5, 50))
windows <- findSignificantWindows(p, times)
```
