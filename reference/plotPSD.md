# Plot power spectral density

Generates a PSD plot for specified channels.

## Usage

``` r
plotPSD(x, channels = NULL, sample = 1L, log_scale = TRUE, freq_range = NULL)
```

## Arguments

- x:

  A PhysioExperiment object.

- channels:

  Integer vector of channel indices. If NULL, plots all.

- sample:

  Integer index for the sample (for 3D data).

- log_scale:

  Logical. If TRUE, uses log scale for power.

- freq_range:

  Numeric vector of length 2 specifying frequency range.

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

[`plotSignal()`](https://x-biosignal.github.io/PhysioAnalysis/reference/plotSignal.md)
for time-domain signal plots,
[`fftSignals()`](https://x-biosignal.github.io/PhysioAnalysis/reference/fftSignals.md)
for FFT computation,
[`bandPower()`](https://x-biosignal.github.io/PhysioAnalysis/reference/bandPower.md)
for frequency band power extraction,
[`plotSpectrogram()`](https://x-biosignal.github.io/PhysioAnalysis/reference/plotSpectrogram.md)
for time-frequency visualization.
