# Spectral parameterization (specparam / FOOOF)

Parameterizes each channel's power spectrum into an aperiodic 1/f
component (offset, exponent, and an optional knee) plus a set of
Gaussian periodic peaks, following the specparam / FOOOF algorithm
(Donoghue et al. 2020). The aperiodic exponent is a marker of
excitation/inhibition balance and cortical state; the peaks capture
genuine oscillations separated from the 1/f background. The fit runs on
the Welch PSD computed by the same machinery as
[`bandPower()`](https://x-biosignal.github.io/PhysioAnalysis/reference/bandPower.md).

## Usage

``` r
specparam(
  x,
  freq_range = c(1, 45),
  peak_width_limits = c(1, 12),
  max_n_peaks = 6L,
  aperiodic_mode = c("fixed", "knee"),
  min_peak_height = 0.05,
  peak_threshold = 2,
  nperseg = NULL,
  assay_name = NULL
)

# S3 method for class 'specparam_result'
print(x, ...)
```

## Arguments

- x:

  A `specparam_result`.

- freq_range:

  Numeric length-2 frequency range in Hz to fit (default: `c(1, 45)`).

- peak_width_limits:

  Numeric length-2 minimum and maximum peak bandwidth in Hz (default:
  `c(1, 12)`).

- max_n_peaks:

  Maximum number of periodic peaks to fit (default: 6).

- aperiodic_mode:

  `"fixed"` for a pure 1/f (offset, exponent) or `"knee"` to also fit a
  knee frequency.

- min_peak_height:

  Minimum peak height above the aperiodic fit, in `log10` power
  (default: 0.05).

- peak_threshold:

  Minimum peak height in units of the standard deviation of the
  flattened spectrum (default: 2).

- nperseg:

  Welch segment length in samples (default: about a 2-second segment for
  adequate low-frequency resolution).

- assay_name:

  Input assay name (default: the default assay).

- ...:

  Ignored.

## Value

An object of class `"specparam_result"`: a list with `aperiodic` (a
data.frame of per-channel offset, exponent, and knee), `peaks` (a
data.frame of per-channel center frequency CF, power PW, and bandwidth
BW), `fit` (per-channel R-squared and error), `spectra` (per-channel log
spectra and fitted model for plotting), and the settings used.

## References

Donoghue, T., et al. (2020). Parameterizing neural power spectra into
periodic and aperiodic components. Nature Neuroscience, 23(12),
1655-1665.

## See also

[`plotSpecparam()`](https://x-biosignal.github.io/PhysioAnalysis/reference/plotSpecparam.md),
[`specparamBiomarker()`](https://x-biosignal.github.io/PhysioAnalysis/reference/specparamBiomarker.md),
[`bandPower()`](https://x-biosignal.github.io/PhysioAnalysis/reference/bandPower.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pe <- make_pe_2d(n_time = 4000, n_channels = 4, sr = 250)
sp <- specparam(pe, freq_range = c(1, 40))
sp$aperiodic
} # }
```
