# Get epoch time vector

Returns the time vector for epoched data relative to event onset.

## Usage

``` r
epochTimes(x)
```

## Arguments

- x:

  An epoched PhysioExperiment object.

## Value

Numeric vector of times in seconds relative to event onset (negative
values represent pre-stimulus time).

## References

Luck, S.J. (2014). "An Introduction to the Event-Related Potential
Technique." 2nd ed. MIT Press.

## See also

[`epochData()`](https://x-biosignal.github.io/PhysioAnalysis/reference/epochData.md)
which creates the epoched data,
[`averageEpochs()`](https://x-biosignal.github.io/PhysioAnalysis/reference/averageEpochs.md)
for epoch averaging,
[`plotERP()`](https://x-biosignal.github.io/PhysioAnalysis/reference/plotERP.md)
for ERP visualization.
