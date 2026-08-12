# Grand average across subjects/samples

Computes grand average across multiple PhysioExperiment objects.

## Usage

``` r
grandAverage(...)
```

## Arguments

- ...:

  PhysioExperiment objects or a list of them.

## Value

A `PhysioExperiment` object with a `"grand_average"` assay containing
the mean across all input objects. Metadata includes `n_subjects`
indicating the number of objects averaged.

## References

Luck, S.J. (2014). "An Introduction to the Event-Related Potential
Technique." 2nd ed. MIT Press.

## See also

[`averageEpochs()`](https://x-biosignal.github.io/PhysioAnalysis/reference/averageEpochs.md)
for within-subject averaging,
[`epochData()`](https://x-biosignal.github.io/PhysioAnalysis/reference/epochData.md)
to create epoched data,
[`plotERP()`](https://x-biosignal.github.io/PhysioAnalysis/reference/plotERP.md)
for ERP visualization.
