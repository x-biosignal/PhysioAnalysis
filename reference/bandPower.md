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
#> 1     Ch1 0.05300117 0.02949889 0.10228801 0.13206227 0.5155690
#> 2     Ch2 0.05452699 0.05939071 0.02381740 0.11298972 0.5845233
#> 3     Ch3 0.02799270 0.04337279 0.05831042 0.09003666 0.4975200
#> 4     Ch4 0.06736209 0.11126386 0.08849311 0.15262704 0.4814074
#> 5     Ch5 0.04508777 0.04211203 0.04159759 0.12080708 0.5253201
#> 6     Ch6 0.06136755 0.02283829 0.01946117 0.08501015 0.4921390

# Compute relative band power
bp_rel <- bandPower(pe, relative = TRUE)
```
