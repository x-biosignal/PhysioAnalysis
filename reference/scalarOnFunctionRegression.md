# Scalar-on-function regression (functional coefficient)

Predicts a scalar outcome from a predictor CURVE (e.g. a clinical score
from a joint-moment waveform). The curve is reduced to functional
principal components, the outcome is regressed on the scores, and the
coefficients are mapped back to a coefficient function beta(t): where
along the curve the predictor matters.

## Usage

``` r
scalarOnFunctionRegression(y, curves, npc = NULL, ve = 0.95)
```

## Arguments

- y:

  Length-`N` numeric outcome.

- curves:

  `N x P` matrix of predictor curves.

- npc:

  Number of functional principal components (default: enough to explain
  `ve` variance).

- ve:

  Variance-explained target for choosing `npc` when `npc` is `NULL`
  (default 0.95).

## Value

a `sofr_result`: `beta` (length-`P` coefficient function), `intercept`,
`fitted`, `r_squared`, `npc`, `scores`, `pve` (proportion variance
explained by the retained PCs).

## References

Reiss PT, Ogden RT (2007) scalar-on-function regression; Ramsay &
Silverman (2005).

## See also

[`functionalRegression()`](https://x-biosignal.github.io/PhysioAnalysis/reference/functionalRegression.md)

## Examples

``` r
set.seed(2)
P <- 80; t <- seq(0, 1, length.out = P)
X <- matrix(rnorm(50 * P), 50, P) + outer(rnorm(50), sin(2 * pi * t))
beta_true <- dnorm(t, 0.3, 0.05)
y <- as.numeric(scale(X %*% beta_true)) + rnorm(50, 0, 0.3)
m <- scalarOnFunctionRegression(y, X, npc = 5)
m$r_squared
#> [1] 0.8919058
```
