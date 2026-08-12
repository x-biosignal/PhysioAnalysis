# Plot topographic map animation

Creates a series of topographic maps across time.

## Usage

``` r
plotTopomapSeries(x, times, ...)
```

## Arguments

- x:

  A PhysioExperiment object with electrode positions.

- times:

  Numeric vector of time points to plot.

- ...:

  Additional arguments passed to
  [`plotTopomap()`](https://x-biosignal.github.io/PhysioAnalysis/reference/plotTopomap.md),
  including `interpolation` and the spherical-spline controls.

## Value

A list of ggplot objects.

## Details

The interpolation method and controls are forwarded unchanged to
[`plotTopomap()`](https://x-biosignal.github.io/PhysioAnalysis/reference/plotTopomap.md).

## References

Shepard, D. (1968). "A two-dimensional interpolation function for
irregularly-spaced data." *Proceedings of the 1968 23rd ACM National
Conference*, 517-524.
[doi:10.1145/800186.810616](https://doi.org/10.1145/800186.810616)

Perrin F, Pernier J, Bertrand O, Echallier J. (1989). Spherical splines
for scalp potential and current density mapping. *Electroencephalography
and Clinical Neurophysiology*, 72(2), 184-187.
[doi:10.1016/0013-4694(89)90180-6](https://doi.org/10.1016/0013-4694%2889%2990180-6)

## See also

[`plotTopomap()`](https://x-biosignal.github.io/PhysioAnalysis/reference/plotTopomap.md)
for a single topographic map,
[`plotERP()`](https://x-biosignal.github.io/PhysioAnalysis/reference/plotERP.md)
for event-related potential plots,
[`plotMultiChannel()`](https://x-biosignal.github.io/PhysioAnalysis/reference/plotMultiChannel.md)
for multi-channel signal visualization.

## Examples

``` r
# \donttest{
pe <- PhysioExperiment(
  assays = list(raw = matrix(rnorm(400), nrow = 100, ncol = 4)),
  colData = S4Vectors::DataFrame(label = c("Fz", "Cz", "Pz", "Oz")),
  samplingRate = 100
)
pe <- applyMontage(pe, "10-20")

# Create topomaps at multiple time points
plots <- plotTopomapSeries(pe, times = c(0.1, 0.2, 0.3, 0.4))
# }
```
