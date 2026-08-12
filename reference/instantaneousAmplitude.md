# Extract instantaneous amplitude (envelope)

Extracts the instantaneous amplitude (envelope) from the analytic
signal.

## Usage

``` r
instantaneousAmplitude(x, assay_name = "analytic", output_assay = "amplitude")
```

## Arguments

- x:

  A PhysioExperiment object with analytic signal.

- assay_name:

  Name of the analytic signal assay.

- output_assay:

  Name for the output assay.

## Value

A `PhysioExperiment` object with an additional assay (default
`"amplitude"`) containing the real-valued instantaneous amplitude
(envelope) of the signal.

## References

Oppenheim, A.V. & Willsky, A.S. (1997). "Signals and Systems." 2nd ed.
Prentice Hall.

## See also

[`hilbertTransform()`](https://x-biosignal.github.io/PhysioAnalysis/reference/hilbertTransform.md)
which must be called first,
[`instantaneousPhase()`](https://x-biosignal.github.io/PhysioAnalysis/reference/instantaneousPhase.md)
for the companion phase extraction,
[`waveletTransform()`](https://x-biosignal.github.io/PhysioAnalysis/reference/waveletTransform.md)
for alternative time-frequency decomposition.

## Examples

``` r
pe <- PhysioExperiment(
  assays = list(raw = matrix(rnorm(500 * 4), nrow = 500)),
  samplingRate = 100
)

# First compute Hilbert transform, then extract amplitude
pe <- hilbertTransform(pe)
pe <- instantaneousAmplitude(pe)
```
