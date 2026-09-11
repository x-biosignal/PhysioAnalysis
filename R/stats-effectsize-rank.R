# Rank-based effect sizes: the rank-biserial correlation (from the
# Mann-Whitney U) and Cliff's delta, which coincide for two independent
# samples, with Cliff's (1993) consistent-variance confidence intervals.

#' Dominance statistic (Cliff's delta) with a consistent-variance CI
#'
#' delta = P(X > Y) - P(X < Y), estimated as the mean of sign(x_i - y_j) over
#' all pairs. The variance uses Cliff's (1993) consistent estimator from the
#' per-row and per-column dominance means, and the interval is his asymmetric
#' (in-bounds) transform so the limits stay within [-1, 1].
#' @keywords internal
#' @noRd
.dominanceStats <- function(x, y, conf_level = 0.95) {
  stopifnot(is.numeric(conf_level), length(conf_level) == 1L,
            is.finite(conf_level), conf_level > 0, conf_level < 1)
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  n1 <- length(x); n2 <- length(y)
  if (n1 < 1L || n2 < 1L) {
    stop("both samples must contain at least one finite value.", call. = FALSE)
  }
  D <- sign(outer(x, y, "-"))                      # n1 x n2 dominance matrix
  delta <- mean(D)
  di <- rowMeans(D); dj <- colMeans(D)
  # Cliff (1993) consistent variance: the numerator uses SUMS of squared
  # deviations about delta (mean(di) = mean(dj) = mean(D) = delta), i.e.
  # SS = (k - 1) * var(.), not the variances themselves.
  ss_di <- if (n1 > 1L) (n1 - 1) * stats::var(di) else 0
  ss_dj <- if (n2 > 1L) (n2 - 1) * stats::var(dj) else 0
  ss_dij <- if (n1 * n2 > 1L) (n1 * n2 - 1) * stats::var(as.vector(D)) else 0
  var_delta <- if (n1 > 1L && n2 > 1L) {
    (n2^2 * ss_di + n1^2 * ss_dj - ss_dij) / (n1 * n2 * (n1 - 1) * (n2 - 1))
  } else {
    NA_real_
  }
  var_delta <- if (is.finite(var_delta)) max(var_delta, 0) else NA_real_

  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  if (is.na(var_delta) || (abs(delta) >= 1 && var_delta == 0)) {
    ci <- c(lower = delta, upper = delta)
  } else {
    se <- sqrt(var_delta)
    denom <- 1 - delta^2 + z^2 * se^2
    num1 <- delta - delta^3
    num2 <- z * se * sqrt((1 - delta^2)^2 + z^2 * se^2)
    ci <- c(lower = max(-1, (num1 - num2) / denom),
            upper = min(1, (num1 + num2) / denom))
  }
  # Mann-Whitney U for group x (favourable pairs, ties count half)
  u1 <- sum(D > 0) + 0.5 * sum(D == 0)
  list(delta = delta, ci = ci, var = var_delta, n1 = n1, n2 = n2,
       u = u1, conf_level = conf_level)
}

#' Rank-biserial correlation for two independent samples
#'
#' The rank-biserial correlation derived from the Mann-Whitney U statistic,
#' equal to `2*U1/(n1*n2) - 1` and to `1 - 2*U2/(n1*n2)` (Kerby 2014), where
#' `U1` counts the pairs in which the first sample exceeds the second. For two
#' independent samples this equals Cliff's delta; the confidence interval is
#' Cliff's (1993) consistent-variance interval.
#'
#' @param x,y Numeric vectors for the two samples.
#' @param conf_level Confidence level (default 0.95).
#' @return A list with `r` (rank-biserial correlation), `ci_lower`, `ci_upper`,
#'   the Mann-Whitney `u` for `x`, and the sample sizes.
#' @references Kerby, D.S. (2014). Compr Psychol 3:11.IT.3.1. Cliff, N. (1993).
#'   Psychol Bull 114(3):494-509.
#' @seealso [cliffsDelta()], [effectSize()]
#' @export
#' @examples
#' rankBiserial(c(5, 6, 7, 8), c(1, 2, 3, 9))
rankBiserial <- function(x, y, conf_level = 0.95) {
  st <- .dominanceStats(as.numeric(x), as.numeric(y), conf_level)
  list(r = st$delta, ci_lower = unname(st$ci["lower"]),
       ci_upper = unname(st$ci["upper"]), u = st$u,
       n1 = st$n1, n2 = st$n2, conf_level = conf_level)
}

#' Cliff's delta for two independent samples
#'
#' Cliff's delta, `delta = P(X > Y) - P(X < Y)`, a non-parametric effect size in
#' `[-1, 1]` whose sign follows the direction of the difference. The confidence
#' interval uses Cliff's (1993) consistent variance with his asymmetric
#' in-bounds transform.
#'
#' @param x,y Numeric vectors for the two samples.
#' @param conf_level Confidence level (default 0.95).
#' @return A list with `delta`, `ci_lower`, `ci_upper`, the variance estimate
#'   and the sample sizes.
#' @references Cliff, N. (1993). Psychol Bull 114(3):494-509.
#' @seealso [rankBiserial()], [effectSize()]
#' @export
#' @examples
#' cliffsDelta(c(5, 6, 7, 8), c(1, 2, 3, 9))
cliffsDelta <- function(x, y, conf_level = 0.95) {
  st <- .dominanceStats(as.numeric(x), as.numeric(y), conf_level)
  list(delta = st$delta, ci_lower = unname(st$ci["lower"]),
       ci_upper = unname(st$ci["upper"]), variance = st$var,
       n1 = st$n1, n2 = st$n2, conf_level = conf_level)
}
