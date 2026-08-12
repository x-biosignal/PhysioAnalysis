# Plot event-related potential (ERP) waveform

Generates an ERP plot from epoched data.

## Usage

``` r
plotERP(x, channel = 1L, ci = 0.95, show_epochs = FALSE, epoch_alpha = 0.2)
```

## Arguments

- x:

  An epoched PhysioExperiment object.

- channel:

  Integer index or character name of the channel.

- ci:

  Confidence interval level (0-1). NULL for no CI.

- show_epochs:

  Logical. If TRUE, shows individual epoch traces.

- epoch_alpha:

  Alpha value for individual epoch traces.

## Value

A ggplot object.

## References

Wickham, H. (2016). *ggplot2: Elegant Graphics for Data Analysis*.
Springer-Verlag New York.
[doi:10.1007/978-3-319-24277-4](https://doi.org/10.1007/978-3-319-24277-4)

Delorme, A. & Makeig, S. (2004). "EEGLAB: an open source toolbox for
analysis of single-trial EEG dynamics including independent component
analysis." *Journal of Neuroscience Methods*, 134(1), 9-21.
[doi:10.1016/j.jneumeth.2003.10.009](https://doi.org/10.1016/j.jneumeth.2003.10.009)

## See also

[`plotMultiChannel()`](https://x-biosignal.github.io/PhysioAnalysis/reference/plotMultiChannel.md)
for multi-channel visualization,
[`plotTopomap()`](https://x-biosignal.github.io/PhysioAnalysis/reference/plotTopomap.md)
for topographic maps,
[`epochData()`](https://x-biosignal.github.io/PhysioAnalysis/reference/epochData.md)
for creating epoched data,
[`bootstrapCI()`](https://x-biosignal.github.io/PhysioAnalysis/reference/bootstrapCI.md)
for bootstrap confidence intervals.
