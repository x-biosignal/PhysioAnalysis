# Functional regression for movement/biomechanics waveforms.
#
# SPM tests differences between groups of curves; functional REGRESSION models a
# curve as a function of covariates (or a scalar outcome as a function of a
# curve). Two directions, both dependency-free base R:
#   * function-on-scalar regression -- the response is a curve, predictors are
#     scalars: fit a pointwise linear model and get COEFFICIENT CURVES with
#     pointwise SE/t and a curve-wide permutation p-value (family-wise across the
#     domain, in the SPM spirit).
#   * scalar-on-function regression -- the response is a scalar, the predictor is
#     a curve: reduce the curve to functional principal components, regress on
#     the scores and map the coefficients back to a COEFFICIENT FUNCTION beta(t).

# Pointwise OLS of Y (N x P) on design D (N x k): returns k x P coefficients etc.
.fosr_fit <- function(Y, D) {
  XtX_inv <- solve(crossprod(D))
  B <- XtX_inv %*% crossprod(D, Y)                   # k x P
  resid <- Y - D %*% B
  dfr <- nrow(D) - ncol(D)
  sigma2 <- colSums(resid^2) / dfr                   # length P
  se <- sqrt(outer(diag(XtX_inv), sigma2))           # k x P
  list(B = B, se = se, tval = B / se, dfr = dfr, sigma2 = sigma2)
}

#' Function-on-scalar regression (coefficient curves)
#'
#' Models a set of response CURVES as a linear function of one or more scalar
#' predictors, giving a coefficient curve for each predictor (how the waveform
#' changes per unit predictor) with pointwise inference and a domain-wide
#' permutation test.
#'
#' @param curves An `N x P` matrix: `N` observations (curves), `P` domain points
#'   (e.g. % gait cycle).
#' @param predictor A length-`N` vector, `N x q` matrix, or data frame of scalar
#'   predictors. An intercept is added automatically.
#' @param n_perm Permutations for the curve-wide (family-wise) p-value per
#'   coefficient (default 1000); set 0 to skip.
#' @param seed Optional integer seed for the permutation test (reproducibility).
#' @return a `fosr_result`: `coefficients` (`k x P`, incl. intercept), `se`,
#'   `tval`, `p_pointwise` (`k x P`), `p_global` (length `k`, family-wise), the
#'   fitted curves and `terms`.
#' @references Ramsay JO, Silverman BW (2005) Functional Data Analysis; Reiss PT,
#'   et al. (2010) function-on-scalar regression.
#' @seealso [scalarOnFunctionRegression()], [spmRegression()]
#' @export
#' @examples
#' set.seed(1)
#' P <- 101; t <- seq(0, 1, length.out = P)
#' speed <- runif(40, 0.8, 1.6)
#' curves <- outer(speed, rep(1, P)) * matrix(sin(2 * pi * t), 40, P, byrow = TRUE) +
#'   matrix(rnorm(40 * P, 0, 0.05), 40, P)
#' fit <- functionalRegression(curves, speed, n_perm = 200)
#' fit$p_global                       # speed coefficient curve is significant
functionalRegression <- function(curves, predictor, n_perm = 1000L, seed = NULL) {
  Y <- as.matrix(curves); N <- nrow(Y); P <- ncol(Y)
  X <- if (is.data.frame(predictor)) as.matrix(predictor) else as.matrix(predictor)
  if (nrow(X) != N) stop("`predictor` must have N rows matching `curves`.", call. = FALSE)
  D <- cbind(`(Intercept)` = 1, X)
  if (is.null(colnames(X))) colnames(D) <- c("(Intercept)", paste0("x", seq_len(ncol(X))))
  k <- ncol(D)
  fit <- .fosr_fit(Y, D)
  p_point <- 2 * stats::pt(-abs(fit$tval), df = fit$dfr)
  # curve-wide family-wise p per coefficient: null of max|t| under row permutation
  p_global <- rep(NA_real_, k)
  if (n_perm > 0) {
    if (!is.null(seed)) set.seed(seed)
    obs_max <- apply(abs(fit$tval), 1, max)
    ge <- integer(k)
    for (b in seq_len(n_perm)) {
      ord <- sample.int(N)                           # random row permutation
      Dp <- D; Dp[, -1] <- D[ord, -1, drop = FALSE]  # permute predictors, keep intercept
      tb <- .fosr_fit(Y, Dp)$tval
      ge <- ge + (apply(abs(tb), 1, max) >= obs_max)
    }
    p_global <- (ge + 1) / (n_perm + 1)
    p_global[1] <- NA_real_                          # intercept not permutation-tested
  }
  names(p_global) <- colnames(D)
  structure(list(
    coefficients = fit$B, se = fit$se, tval = fit$tval,
    p_pointwise = p_point, p_global = p_global,
    fitted = D %*% fit$B, terms = colnames(D), df_residual = fit$dfr,
    N = N, P = P), class = "fosr_result")
}

#' @export
print.fosr_result <- function(x, ...) {
  cat(sprintf("Function-on-scalar regression -- %d curves x %d points, %d terms\n",
              x$N, x$P, length(x$terms)))
  for (j in seq_along(x$terms)) {
    pg <- x$p_global[j]
    cat(sprintf("  %-14s peak|t| = %5.2f  global p = %s\n", x$terms[j],
                max(abs(x$tval[j, ])), if (is.na(pg)) "-" else sprintf("%.3f", pg)))
  }
  invisible(x)
}

#' Scalar-on-function regression (functional coefficient)
#'
#' Predicts a scalar outcome from a predictor CURVE (e.g. a clinical score from a
#' joint-moment waveform). The curve is reduced to functional principal
#' components, the outcome is regressed on the scores, and the coefficients are
#' mapped back to a coefficient function beta(t): where along the curve the
#' predictor matters.
#'
#' @param y Length-`N` numeric outcome.
#' @param curves `N x P` matrix of predictor curves.
#' @param npc Number of functional principal components (default: enough to
#'   explain `ve` variance).
#' @param ve Variance-explained target for choosing `npc` when `npc` is `NULL`
#'   (default 0.95).
#' @return a `sofr_result`: `beta` (length-`P` coefficient function), `intercept`,
#'   `fitted`, `r_squared`, `npc`, `scores`, `pve` (proportion variance explained
#'   by the retained PCs).
#' @references Reiss PT, Ogden RT (2007) scalar-on-function regression; Ramsay &
#'   Silverman (2005).
#' @seealso [functionalRegression()]
#' @export
#' @examples
#' set.seed(2)
#' P <- 80; t <- seq(0, 1, length.out = P)
#' X <- matrix(rnorm(50 * P), 50, P) + outer(rnorm(50), sin(2 * pi * t))
#' beta_true <- dnorm(t, 0.3, 0.05)
#' y <- as.numeric(scale(X %*% beta_true)) + rnorm(50, 0, 0.3)
#' m <- scalarOnFunctionRegression(y, X, npc = 5)
#' m$r_squared
scalarOnFunctionRegression <- function(y, curves, npc = NULL, ve = 0.95) {
  y <- as.numeric(y); X <- as.matrix(curves); N <- nrow(X); P <- ncol(X)
  if (length(y) != N) stop("`y` must have length N (rows of `curves`).", call. = FALSE)
  cm <- colMeans(X); Xc <- sweep(X, 2L, cm)
  sv <- svd(Xc)
  ev <- sv$d^2; pve_cum <- cumsum(ev) / sum(ev)
  if (is.null(npc)) npc <- max(1L, which(pve_cum >= ve)[1])
  npc <- min(npc, length(sv$d))
  V <- sv$v[, seq_len(npc), drop = FALSE]            # P x npc loadings
  scores <- Xc %*% V                                 # N x npc
  df <- data.frame(y = y, scores)
  fit <- stats::lm(y ~ ., data = df)
  b <- stats::coef(fit)[-1]                          # npc coefficients
  b[is.na(b)] <- 0
  beta <- as.numeric(V %*% b)                        # P coefficient function
  fitted <- as.numeric(stats::predict(fit))
  structure(list(
    beta = beta, intercept = mean(y) - sum(cm * beta), fitted = fitted,
    r_squared = summary(fit)$r.squared, npc = npc, scores = scores,
    pve = pve_cum[npc], N = N, P = P), class = "sofr_result")
}

#' @export
print.sofr_result <- function(x, ...) {
  cat(sprintf("Scalar-on-function regression -- %d obs x %d points\n", x$N, x$P))
  cat(sprintf("  %d functional PCs (%.1f%% variance)   R-squared = %.3f\n",
              x$npc, 100 * x$pve, x$r_squared))
  invisible(x)
}
