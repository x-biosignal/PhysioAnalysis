# Functional random-intercept (multilevel) model for waveforms

Fits the subject random-intercept model to a set of curves nested in
subjects. It splits the pointwise variance into between- and
within-subject components (a functional ICC curve) and, when `group` is
given, tests the between-subject fixed effect on the subject-mean curves
so strides are not treated as independent.

## Usage

``` r
functionalMixedModel(
  curves,
  subject,
  group = NULL,
  n_perm = 1000L,
  seed = NULL
)
```

## Arguments

- curves:

  An `N x P` matrix of curves: `N` observations (e.g. strides pooled
  over subjects), `P` domain points.

- subject:

  Length-`N` grouping identifying the random-intercept level (the
  subject each curve belongs to).

- group:

  Optional between-subject fixed effect (length `N`), constant within
  each subject (e.g. a patient/control label).

- n_perm:

  Permutations for the curve-wide fixed-effect p-value (default 1000; 0
  to skip).

- seed:

  Optional RNG seed for the permutation test.

## Value

a `functional_mixed` object: `mean_curve`, `var_between`, `var_within`,
`icc` (functional ICC curve), `subject_means` (one row per subject),
`n_subjects`, and – when `group` is given – `f_curve`, `p_pointwise`,
`p_global` (subject-level fixed-effect test) plus `group_diff` (the
signed difference curve for two groups).

## References

Morris JS (2015) Functional regression, Annu Rev Stat Appl; Pini &
Vantini (2017) interval-wise testing.

## See also

[`functionalRegression()`](https://x-biosignal.github.io/PhysioAnalysis/reference/functionalRegression.md),
[`spmTTest()`](https://x-biosignal.github.io/PhysioAnalysis/reference/spmTTest.md)

## Examples

``` r
set.seed(1)
P <- 101; t <- seq(0, 1, length.out = P)
make <- function(off) t(sapply(1:12, function(k) off + sin(2 * pi * t) + rnorm(P, 0, 0.1)))
ctrl <- do.call(rbind, lapply(rnorm(6, 0, 0.5), make))
pat  <- do.call(rbind, lapply(rnorm(6, 1.0, 0.5), make))   # shifted group
Y <- rbind(ctrl, pat)
subj <- rep(1:12, each = 12); grp <- rep(c("ctrl", "pat"), each = 72)
fm <- functionalMixedModel(Y, subj, grp, n_perm = 200)
fm$p_global                         # group difference, subject-level correct
#> [1] 0.0199005
```
