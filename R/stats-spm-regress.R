# Regression SPM: an SPM{t} field for the slope of a continuous covariate at
# each node (spm1d.stats.regress). Reuses the RFT machinery in stats-spm.R
# (.estimateFWHMresiduals / .reselCount / .rftCriticalT / cluster inference).

#' SPM linear-regression t-field
#'
#' Computes an SPM\{t\} field testing, at every node, the slope of a linear
#' regression of the waveform on a continuous predictor. The node-wise statistic
#' is the ordinary least-squares slope t-value (identical to
#' \code{summary(lm(y ~ predictor))} at that node), and field-level significance
#' uses the same random-field-theory threshold as [spmTTest()] with
#' \code{df = n - 2}.
#'
#' @param x A \code{PhysioExperiment} or a numeric matrix (time x observations).
#' @param predictor Numeric covariate, one value per observation (column).
#' @param alpha Significance level (default 0.05).
#' @param two_tailed Logical; two-tailed slope test (default \code{TRUE}).
#' @return A list of class \code{"spm_result"} (\code{test_type = "regression"})
#'   with the \code{t} field, RFT \code{threshold}, significant \code{clusters},
#'   pointwise \code{p_values}, \code{df}, \code{fwhm}, \code{resel_count}, and
#'   the fitted \code{slope}/\code{intercept} fields.
#' @references Pataky 2016; Friston et al. 2007 (RFT). spm1d.stats.regress.
#' @seealso [spmTTest()], [spmMANOVA()], [spmSnPM()]
#' @export
#' @examples
#' set.seed(1)
#' pred <- rnorm(20)
#' data <- matrix(rnorm(100 * 20), nrow = 100)
#' data[40:60, ] <- data[40:60, ] + outer(rep(1.5, 21), pred)  # slope 30% window
#' spmRegression(data, pred)
spmRegression <- function(x, predictor, alpha = 0.05, two_tailed = TRUE) {
  if (inherits(x, "PhysioExperiment")) {
    data <- SummarizedExperiment::assay(x, defaultAssay(x))
  } else if (is.matrix(x)) {
    data <- x
  } else {
    stop("Input must be a PhysioExperiment or matrix", call. = FALSE)
  }
  n_time <- nrow(data); n_obs <- ncol(data)
  predictor <- as.numeric(predictor)
  if (length(predictor) != n_obs) {
    stop("'predictor' must have one value per observation (", n_obs, ").",
         call. = FALSE)
  }
  if (n_obs < 3L) stop("regression SPM needs at least 3 observations.",
                       call. = FALSE)
  xc <- predictor - mean(predictor)
  sxx <- sum(xc^2)
  if (!is.finite(sxx) || sxx <= 0) {
    stop("'predictor' has zero variance.", call. = FALSE)
  }

  # node-wise OLS slope, vectorised over time
  ybar <- rowMeans(data, na.rm = TRUE)
  slope <- as.numeric(data %*% xc) / sxx                 # sum(xc * y) / Sxx
  intercept <- ybar - slope * mean(predictor)
  fitted <- outer(ybar, rep(1, n_obs)) + outer(slope, xc)
  residuals <- data - fitted                             # n_time x n_obs
  sse <- rowSums(residuals^2)
  df <- n_obs - 2L
  mse <- sse / df
  se_slope <- sqrt(mse / sxx)
  t_stat <- slope / se_slope
  t_stat[!is.finite(t_stat)] <- 0

  fwhm <- .estimateFWHMresiduals(t(residuals))
  resel_count <- .reselCount(n_time, fwhm)
  threshold <- .rftCriticalT(alpha, df, resel_count, two_tailed)

  clusters <- .findSPMClusters(t_stat, threshold, two_tailed)
  if (length(clusters) > 0) {
    ec <- .ecDensityT(threshold, df)
    for (i in seq_along(clusters)) {
      clusters[[i]]$p_cluster <- .rftClusterPValue(
        clusters[[i]]$extent, resel_count, fwhm, ec["rho0"], ec["rho1"],
        tails = if (two_tailed) 2 else 1)
    }
  }
  p_values <- if (two_tailed) 2 * stats::pt(-abs(t_stat), df) else
    stats::pt(-t_stat, df)

  structure(list(
    t = t_stat, threshold = threshold, clusters = clusters,
    p_values = p_values, df = df, fwhm = fwhm, resel_count = resel_count,
    alpha = alpha, two_tailed = two_tailed, n_time = n_time,
    slope = slope, intercept = intercept, test_type = "regression"),
    class = c("spm_result", "list"))
}
