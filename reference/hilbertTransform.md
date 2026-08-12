# Hilbert transform for instantaneous amplitude/phase

Computes the analytic signal using the Hilbert transform. The analytic
signal can be used to extract instantaneous amplitude and phase.

## Usage

``` r
hilbertTransform(x, output_assay = "analytic")
```

## Arguments

- x:

  A PhysioExperiment object.

- output_assay:

  Name for the output assay.

## Value

A `PhysioExperiment` object with an additional assay (default
`"analytic"`) containing the complex-valued analytic signal.

## References

Oppenheim, A.V. & Willsky, A.S. (1997). "Signals and Systems." 2nd ed.
Prentice Hall.

## See also

[`instantaneousAmplitude()`](https://x-biosignal.github.io/PhysioAnalysis/reference/instantaneousAmplitude.md)
to extract the signal envelope,
[`instantaneousPhase()`](https://x-biosignal.github.io/PhysioAnalysis/reference/instantaneousPhase.md)
to extract the instantaneous phase,
[`plv()`](https://x-biosignal.github.io/PhysioAnalysis/reference/plv.md)
for phase-based connectivity.

## Examples

``` r
pe <- PhysioExperiment(
  assays = list(raw = matrix(rnorm(500 * 4), nrow = 500)),
  samplingRate = 100
)

# Compute Hilbert transform
pe <- hilbertTransform(pe)

# Extract amplitude and phase
pe <- instantaneousAmplitude(pe)
pe <- instantaneousPhase(pe)
```
