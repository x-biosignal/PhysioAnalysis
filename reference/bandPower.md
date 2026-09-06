# Compute band power

Extracts power in specified frequency bands.

## Usage

``` r
bandPower(x, bands = NULL, method = c("welch", "wavelet"), relative = FALSE)
```

## Arguments

- x:

  A PhysioExperiment object.

- bands:

  Named list of frequency bands. Each element should be c(low, high).
  Default includes standard EEG bands.

- method:

  Method: "welch" (PSD) or "wavelet".

- relative:

  If TRUE, returns relative power (proportion of total).

## Value

A `data.frame` with one row per channel and columns for the channel name
and power in each frequency band. If `relative = TRUE`, values represent
the proportion of total power.

## References

Oppenheim, A.V. & Willsky, A.S. (1997). "Signals and Systems." 2nd ed.
Prentice Hall.

## See also

[`spectrogram()`](https://x-biosignal.github.io/PhysioAnalysis/reference/spectrogram.md)
for full time-frequency decomposition,
[`waveletTransform()`](https://x-biosignal.github.io/PhysioAnalysis/reference/waveletTransform.md)
for wavelet-based power,
[`fftSignals()`](https://x-biosignal.github.io/PhysioAnalysis/reference/fftSignals.md)
for raw FFT,
[`plotPSD()`](https://x-biosignal.github.io/PhysioAnalysis/reference/plotPSD.md)
for power spectral density visualization.

## Examples

``` r
pe <- PhysioExperiment(
  assays = list(raw = matrix(rnorm(2560), nrow = 256, ncol = 10)),
  colData = S4Vectors::DataFrame(label = paste0("Ch", 1:10)),
  samplingRate = 256
)

# Compute band power for standard EEG bands
bp <- bandPower(pe)
head(bp)
#>   channel      delta      theta      alpha       beta     gamma
#> 1     Ch1 0.05834715 0.04624460 0.02555233 0.11101998 0.4392125
#> 2     Ch2 0.02116564 0.01409393 0.05685884 0.20376474 0.5593013
#> 3     Ch3 0.04033103 0.01007837 0.01538094 0.18172863 0.6396738
#> 4     Ch4 0.01202950 0.05605040 0.04575043 0.10880683 0.5550231
#> 5     Ch5 0.03958326 0.09531850 0.02874965 0.07207095 0.6544251
#> 6     Ch6 0.06042356 0.07204888 0.14383959 0.17269528 0.4446339

# Compute relative band power
bp_rel <- bandPower(pe, relative = TRUE)
```
