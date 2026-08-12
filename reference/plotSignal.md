# Plot a signal channel

Generates a simple line plot for the selected channel and sample from
the default assay. Supports both 2D (time x channel) and 3D (time x
channel x sample) assay arrays.

## Usage

``` r
plotSignal(x, channel = 1L, sample = 1L, assay_name = NULL)
```

## Arguments

- x:

  A `PhysioExperiment` object.

- channel:

  Integer index for the channel.

- sample:

  Integer index for the sample (only used for 3D arrays).

- assay_name:

  Optional assay name. If NULL, uses the default assay.

## Value

A `ggplot` object.

## References

Wickham, H. (2016). *ggplot2: Elegant Graphics for Data Analysis*.
Springer-Verlag New York.
[doi:10.1007/978-3-319-24277-4](https://doi.org/10.1007/978-3-319-24277-4)

## See also

[`plotMultiChannel()`](https://x-biosignal.github.io/PhysioAnalysis/reference/plotMultiChannel.md)
for multi-channel visualization,
[`plotPSD()`](https://x-biosignal.github.io/PhysioAnalysis/reference/plotPSD.md)
for power spectral density plots,
[`plotERP()`](https://x-biosignal.github.io/PhysioAnalysis/reference/plotERP.md)
for event-related potential plots.
