# Multiple comparison correction

Applies multiple comparison correction to p-values.

## Usage

``` r
correctPValues(p_values, method = c("fdr", "bonferroni", "holm", "bh", "none"))
```

## Arguments

- p_values:

  Matrix or vector of p-values.

- method:

  Correction method: "bonferroni", "holm", "fdr" (Benjamini-Hochberg),
  "bh" (alias for fdr), or "none".

## Value

Corrected p-values in the same format as input.

## References

Maris, E. & Oostenveld, R. (2007). "Nonparametric statistical testing of
EEG- and MEG-data." *Journal of Neuroscience Methods*, 164(1), 177-190.
[doi:10.1016/j.jneumeth.2007.03.024](https://doi.org/10.1016/j.jneumeth.2007.03.024)

## See also

[`clusterPermutationTest()`](https://x-biosignal.github.io/PhysioAnalysis/reference/clusterPermutationTest.md)
for cluster-based correction,
[`tTestEpochs()`](https://x-biosignal.github.io/PhysioAnalysis/reference/tTestEpochs.md)
for pointwise t-tests,
[`findSignificantWindows()`](https://x-biosignal.github.io/PhysioAnalysis/reference/findSignificantWindows.md)
for identifying significant time windows.

## Examples

``` r
# Example p-values
p <- matrix(runif(100), nrow = 10)
# Apply FDR correction
p_corrected <- correctPValues(p, method = "fdr")
```
