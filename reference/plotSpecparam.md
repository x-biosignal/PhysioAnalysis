# Plot a specparam model fit

Plots the observed log power spectrum for one channel together with the
full specparam model and the aperiodic-only fit.

## Usage

``` r
plotSpecparam(result, channel = 1)
```

## Arguments

- result:

  A `specparam_result` from
  [`specparam()`](https://x-biosignal.github.io/PhysioAnalysis/reference/specparam.md).

- channel:

  Channel index or label to plot (default: 1).

## Value

A `ggplot` object.

## See also

[`specparam()`](https://x-biosignal.github.io/PhysioAnalysis/reference/specparam.md)

## Examples

``` r
if (FALSE) { # \dontrun{
sp <- specparam(make_pe_2d(n_time = 4000, sr = 250))
plotSpecparam(sp, channel = 1)
} # }
```
