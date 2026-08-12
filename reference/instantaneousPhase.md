# Extract instantaneous phase

Extracts the instantaneous phase from the analytic signal.

## Usage

``` r
instantaneousPhase(x, assay_name = "analytic", output_assay = "phase")
```

## Arguments

- x:

  A PhysioExperiment object with analytic signal.

- assay_name:

  Name of the analytic signal assay.

- output_assay:

  Name for the output assay.

## Value

A `PhysioExperiment` object with an additional assay (default `"phase"`)
containing instantaneous phase values in radians (range \\\[-\pi,
\pi\]\\).

## References

Oppenheim, A.V. & Willsky, A.S. (1997). "Signals and Systems." 2nd ed.
Prentice Hall.

## See also

[`hilbertTransform()`](https://x-biosignal.github.io/PhysioAnalysis/reference/hilbertTransform.md)
which must be called first,
[`instantaneousAmplitude()`](https://x-biosignal.github.io/PhysioAnalysis/reference/instantaneousAmplitude.md)
for the companion amplitude extraction,
[`plv()`](https://x-biosignal.github.io/PhysioAnalysis/reference/plv.md)
for phase-based connectivity analysis.

## Examples

``` r
pe <- PhysioExperiment(
  assays = list(raw = matrix(rnorm(500 * 4), nrow = 500)),
  samplingRate = 100
)

# First compute Hilbert transform, then extract phase
pe <- hilbertTransform(pe)
pe <- instantaneousPhase(pe)
```
