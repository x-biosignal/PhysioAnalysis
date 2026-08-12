# Average epochs

Computes the average across epochs, optionally by condition.

## Usage

``` r
averageEpochs(x, by = NULL)
```

## Arguments

- x:

  An epoched PhysioExperiment object.

- by:

  Optional column name in epoch_info to group by.

## Value

A `PhysioExperiment` object with an `"averaged"` assay containing the
mean signal across epochs. If `by` is specified, the third dimension
corresponds to conditions.

## References

Luck, S.J. (2014). "An Introduction to the Event-Related Potential
Technique." 2nd ed. MIT Press.

## See also

[`epochData()`](https://x-biosignal.github.io/PhysioAnalysis/reference/epochData.md)
to create epoched data,
[`grandAverage()`](https://x-biosignal.github.io/PhysioAnalysis/reference/grandAverage.md)
for averaging across subjects,
[`plotERP()`](https://x-biosignal.github.io/PhysioAnalysis/reference/plotERP.md)
for ERP visualization.
