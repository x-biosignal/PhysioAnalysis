# Compute temporal stability of network

Calculates how stable the network topology is over time.

## Usage

``` r
temporalStability(dyn_conn, metric = c("correlation", "distance"))
```

## Arguments

- dyn_conn:

  Result from slidingWindowConnectivity().

- metric:

  Stability metric: "correlation" or "distance".

## Value

A list with stability metrics.

## Examples

``` r
set.seed(123)
pe <- PhysioExperiment(
  assays = list(raw = matrix(rnorm(1000 * 4), nrow = 1000, ncol = 4)),
  samplingRate = 100
)
dyn_conn <- slidingWindowConnectivity(pe, window_size = 200, step = 50)
stability <- temporalStability(dyn_conn)
```
