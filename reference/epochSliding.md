# Create epochs using sliding window

Segments continuous data into overlapping epochs using a sliding window
approach. This is useful for time-frequency analysis or when events are
not available.

## Usage

``` r
epochSliding(x, window, step = NULL, baseline = NULL)
```

## Arguments

- x:

  A PhysioExperiment object

- window:

  Window size in seconds

- step:

  Step size in seconds (default: window/2 for 50% overlap)

- baseline:

  Optional baseline correction window as c(start, end) in seconds
  relative to epoch start. NULL for no correction.

## Value

PhysioExperiment with epoched data

## Examples

``` r
pe <- PhysioExperiment(
  assays = list(raw = matrix(rnorm(1000 * 4), nrow = 1000, ncol = 4)),
  samplingRate = 100)
# Create 0.5 second windows with 0.1 second step
epoched <- epochSliding(pe, window = 0.5, step = 0.1)
```
