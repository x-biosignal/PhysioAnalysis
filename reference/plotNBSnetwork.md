# Plot a significant NBS subnetwork

Draws the edges belonging to significant NBS components on a circular
node layout.

## Usage

``` r
plotNBSnetwork(
  result,
  alpha = NULL,
  node_labels = NULL,
  title = "NBS subnetwork"
)
```

## Arguments

- result:

  A result from
  [`networkBasedStatistic()`](https://x-biosignal.github.io/PhysioAnalysis/reference/networkBasedStatistic.md).

- alpha:

  Significance level for the components to draw (default: uses the value
  stored in `result`).

- node_labels:

  Optional character node labels.

- title:

  Plot title.

## Value

A `ggplot` object.

## See also

[`networkBasedStatistic()`](https://x-biosignal.github.io/PhysioAnalysis/reference/networkBasedStatistic.md)

## Examples

``` r
if (FALSE) { # \dontrun{
plotNBSnetwork(networkBasedStatistic(g1, g2))
} # }
```
