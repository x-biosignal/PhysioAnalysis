# Sliding window connectivity

Computes connectivity matrices over sliding time windows.

## Usage

``` r
slidingWindowConnectivity(
  x,
  window_size = 256L,
  step = 64L,
  method = c("correlation", "coherence", "plv"),
  ...
)
```

## Arguments

- x:

  A PhysioExperiment object.

- window_size:

  Window size in samples.

- step:

  Step size in samples.

- method:

  Connectivity method: "correlation", "coherence", or "plv".

- ...:

  Additional arguments passed to connectivity function.

## Value

A list with time-varying connectivity matrices.

## Examples

``` r
set.seed(123)
pe <- PhysioExperiment(
  assays = list(raw = matrix(rnorm(1000 * 4), nrow = 1000, ncol = 4)),
  samplingRate = 100
)
dyn_conn <- slidingWindowConnectivity(pe, window_size = 200, step = 50)
```
