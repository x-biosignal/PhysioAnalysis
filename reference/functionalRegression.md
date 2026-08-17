# Function-on-scalar regression (coefficient curves)

Models a set of response CURVES as a linear function of one or more
scalar predictors, giving a coefficient curve for each predictor (how
the waveform changes per unit predictor) with pointwise inference and a
domain-wide permutation test.

## Usage

``` r
functionalRegression(curves, predictor, n_perm = 1000L, seed = NULL)
```

## Arguments

- curves:

  An `N x P` matrix: `N` observations (curves), `P` domain points (e.g.
  % gait cycle).

- predictor:

  A length-`N` vector, `N x q` matrix, or data frame of scalar
  predictors. An intercept is added automatically.

- n_perm:

  Permutations for the curve-wide (family-wise) p-value per coefficient
  (default 1000); set 0 to skip.

- seed:

  Optional integer seed for the permutation test (reproducibility).

## Value

a `fosr_result`: `coefficients` (`k x P`, incl. intercept), `se`,
`tval`, `p_pointwise` (`k x P`), `p_global` (length `k`, family-wise),
the fitted curves and `terms`.

## References

Ramsay JO, Silverman BW (2005) Functional Data Analysis; Reiss PT, et
al. (2010) function-on-scalar regression.

## See also

[`scalarOnFunctionRegression()`](https://x-biosignal.github.io/PhysioAnalysis/reference/scalarOnFunctionRegression.md),
[`spmRegression()`](https://x-biosignal.github.io/PhysioAnalysis/reference/spmRegression.md)

## Examples

``` r
set.seed(1)
P <- 101; t <- seq(0, 1, length.out = P)
speed <- runif(40, 0.8, 1.6)
curves <- outer(speed, rep(1, P)) * matrix(sin(2 * pi * t), 40, P, byrow = TRUE) +
  matrix(rnorm(40 * P, 0, 0.05), 40, P)
fit <- functionalRegression(curves, speed, n_perm = 200)
fit$p_global                       # speed coefficient curve is significant
#> (Intercept)          x1 
#>          NA 0.004975124 
```
