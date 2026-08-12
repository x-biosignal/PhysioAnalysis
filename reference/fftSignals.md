# Fast Fourier transform helper

Computes the discrete Fourier transform along the time axis of the
default assay and stores the magnitude spectrum in a new assay named
`"fft"`.

## Usage

``` r
fftSignals(x)
```

## Arguments

- x:

  A `PhysioExperiment` object.

## Value

A `PhysioExperiment` object with an additional `"fft"` assay containing
the magnitude spectrum (same dimensions as the input assay).

## References

Oppenheim, A.V. & Willsky, A.S. (1997). "Signals and Systems." 2nd ed.
Prentice Hall.

## See also

[`spectrogram()`](https://x-biosignal.github.io/PhysioAnalysis/reference/spectrogram.md)
for time-frequency analysis,
[`bandPower()`](https://x-biosignal.github.io/PhysioAnalysis/reference/bandPower.md)
for frequency band power extraction,
[`plotPSD()`](https://x-biosignal.github.io/PhysioAnalysis/reference/plotPSD.md)
for power spectral density visualization.
