# Aperiodic parameter as a PhysioBiomarker

Wraps a specparam aperiodic parameter (by default the exponent) for one
channel as a reliability-characterised
[`PhysioBiomarker`](https://x-biosignal.github.io/PhysioCore//reference/PhysioBiomarker.html).

## Usage

``` r
specparamBiomarker(
  result,
  channel = 1,
  param = c("exponent", "offset", "knee"),
  reliability = list(icc = NA_real_, sem = NA_real_)
)
```

## Arguments

- result:

  A `specparam_result` from
  [`specparam()`](https://x-biosignal.github.io/PhysioAnalysis/reference/specparam.md).

- channel:

  Channel index or label (default: 1).

- param:

  Which aperiodic parameter to export: `"exponent"`, `"offset"`, or
  `"knee"`.

- reliability:

  Named list of reliability indices to attach (default: `icc` and `sem`
  placeholders).

## Value

A
[`PhysioBiomarker`](https://x-biosignal.github.io/PhysioCore//reference/PhysioBiomarker.html).

## See also

[`specparam()`](https://x-biosignal.github.io/PhysioAnalysis/reference/specparam.md),
[`PhysioCore::physioBiomarker()`](https://x-biosignal.github.io/PhysioCore//reference/physioBiomarker-constructor.html)

## Examples

``` r
if (FALSE) { # \dontrun{
sp <- specparam(make_pe_2d(n_time = 4000, sr = 250))
specparamBiomarker(sp, channel = 1)
} # }
```
