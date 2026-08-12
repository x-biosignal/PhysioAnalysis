# Plot spectrogram

Creates a visualization of a spectrogram result.

## Usage

``` r
plotSpectrogram(spec, freq_range = NULL, log_power = TRUE)
```

## Arguments

- spec:

  Spectrogram result from spectrogram().

- freq_range:

  Optional frequency range to display.

- log_power:

  If TRUE, displays log power (dB scale).

## Value

A `ggplot` object showing the spectrogram as a filled raster plot with
time on the x-axis and frequency on the y-axis.

## References

Oppenheim, A.V. & Willsky, A.S. (1997). "Signals and Systems." 2nd ed.
Prentice Hall.

## See also

[`spectrogram()`](https://x-biosignal.github.io/PhysioAnalysis/reference/spectrogram.md)
to compute the spectrogram data,
[`plotPSD()`](https://x-biosignal.github.io/PhysioAnalysis/reference/plotPSD.md)
for power spectral density plots,
[`plotSignal()`](https://x-biosignal.github.io/PhysioAnalysis/reference/plotSignal.md)
for time-domain visualization.

## Examples

``` r
pe <- PhysioExperiment(
  assays = list(raw = matrix(rnorm(1000 * 4), nrow = 1000)),
  samplingRate = 250
)

# Compute spectrogram
spec <- spectrogram(pe, channel = 1)

# Plot with frequency range filter
plotSpectrogram(spec, freq_range = c(1, 50))
```
