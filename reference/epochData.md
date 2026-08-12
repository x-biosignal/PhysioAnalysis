# Epoching operations for PhysioExperiment

Functions for segmenting continuous data into epochs/trials based on
events. Epoch data around events

## Usage

``` r
epochData(
  x,
  tmin = -0.2,
  tmax = 0.8,
  event_type = NULL,
  baseline = NULL,
  reject = NULL,
  events = NULL,
  min_length = NULL
)
```

## Arguments

- x:

  A PhysioExperiment object.

- tmin:

  Time before event onset in seconds (negative for pre-stimulus).

- tmax:

  Time after event onset in seconds.

- event_type:

  Character vector of event types to epoch around. If NULL, uses all
  events. Ignored if `events` is provided.

- baseline:

  Numeric vector of length 2 specifying baseline period (tmin, tmax) for
  baseline correction. NULL for no correction.

- reject:

  Amplitude threshold for epoch rejection. NULL to keep all.

- events:

  An EventQuery object for advanced event filtering. If provided,
  overrides `event_type`.

- min_length:

  Minimum epoch length in seconds when using variable-length epochs
  (tmax as event name). Epochs shorter than this are excluded.

## Value

A new `PhysioExperiment` object with a 4D assay (time x channel x epoch
x sample) named `"epoched"`. Metadata includes `epoch_tmin`,
`epoch_tmax`, `epoch_info` (a DataFrame with epoch_id, event_type,
event_value, event_onset), and `n_epochs`.

## Details

Extracts epochs (segments) of data around specified events.

## References

Luck, S.J. (2014). "An Introduction to the Event-Related Potential
Technique." 2nd ed. MIT Press.

## See also

[`averageEpochs()`](https://x-biosignal.github.io/PhysioAnalysis/reference/averageEpochs.md)
to average across epochs,
[`epochTimes()`](https://x-biosignal.github.io/PhysioAnalysis/reference/epochTimes.md)
to get the time vector,
[`tTestEpochs()`](https://x-biosignal.github.io/PhysioAnalysis/reference/tTestEpochs.md)
for statistical testing on epoched data,
[`plotERP()`](https://x-biosignal.github.io/PhysioAnalysis/reference/plotERP.md)
for ERP visualization.

## Examples

``` r
# Create continuous data with events
pe <- PhysioExperiment(
  assays = list(raw = matrix(rnorm(1000 * 4), nrow = 1000)),
  samplingRate = 100
)
pe <- addEvents(pe, onset = c(1, 2, 3, 4, 5), type = "stimulus")

# Extract epochs: 200ms before to 800ms after stimulus
epochs <- epochData(pe, tmin = -0.2, tmax = 0.8)

# With baseline correction
epochs_bl <- epochData(pe, tmin = -0.2, tmax = 0.8,
                       baseline = c(-0.2, 0))

# With artifact rejection
epochs_clean <- epochData(pe, tmin = -0.2, tmax = 0.8,
                          reject = 100)
```
