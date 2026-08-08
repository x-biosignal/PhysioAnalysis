#' Statistical Parametric Mapping (SPM1D) for Biomechanics
#'
#' Functions for statistical analysis of continuous waveform data using
#' Statistical Parametric Mapping (SPM) methodology adapted from neuroimaging.
#' These methods test hypotheses over entire waveforms rather than discrete points.

#' SPM t-test for waveform comparison
#'
#' Performs a t-test at each time point and computes an SPM t-statistic map
#' with Random Field Theory (RFT) correction for multiple comparisons.
#'
#' @param x A PhysioExperiment object or matrix (time x observations).
#' @param group1 Indices for first group (for two-sample test).
#' @param group2 Indices for second group. If NULL, performs one-sample test.
#' @param alpha Significance level (default: 0.05).
#' @param two_tailed Logical; if TRUE, performs two-tailed test.
#'
#' @return A list of class "spm_result" containing:
#'   \item{t}{T-statistic at each time point}
#'   \item{threshold}{Critical threshold from RFT}
#'   \item{clusters}{Significant clusters (start, end, extent, p-value)}
#'   \item{p_values}{Pointwise p-values}
#'   \item{alpha}{Significance level used}
#'
#' @details
#' SPM analyzes continuous biomechanical waveforms (e.g., joint angles, moments)
#' by computing t-statistics at each time point and using Random Field Theory
#' to control family-wise error rate across the entire waveform.
#'
#' @references
#' Pataky TC (2012). One-dimensional statistical parametric mapping in Python.
#' Computer Methods in Biomechanics and Biomedical Engineering.
#'
#' @export
#' @examples
#' # Create example gait data (100 time points x 20 subjects)
#' set.seed(123)
#' # Group 1: normal gait
#' g1 <- matrix(rnorm(100 * 10), nrow = 100, ncol = 10)
#' # Group 2: altered gait (effect at 40-60% of cycle)
#' g2 <- matrix(rnorm(100 * 10), nrow = 100, ncol = 10)
#' g2[40:60, ] <- g2[40:60, ] + 1.5
#'
#' data <- cbind(g1, g2)
#' pe <- PhysioExperiment(
#'   assays = list(values = data),
#'   samplingRate = 100
#' )
#'
#' # Two-sample SPM t-test
#' result <- spmTTest(pe, group1 = 1:10, group2 = 11:20)
#' print(result)
spmTTest <- function(x, group1 = NULL, group2 = NULL,
                     alpha = 0.05, two_tailed = TRUE) {

  # Extract data matrix
 if (inherits(x, "PhysioExperiment")) {
    assay_name <- defaultAssay(x)
    data <- SummarizedExperiment::assay(x, assay_name)
  } else if (is.matrix(x)) {
    data <- x
  } else {
    stop("Input must be a PhysioExperiment or matrix", call. = FALSE)
  }

  n_time <- nrow(data)
  n_obs <- ncol(data)

  # Determine test type
  if (is.null(group1)) {
    group1 <- seq_len(n_obs)
  }

  two_sample <- !is.null(group2)

  if (two_sample) {
    # Two-sample t-test
    data1 <- data[, group1, drop = FALSE]
    data2 <- data[, group2, drop = FALSE]
    n1 <- ncol(data1)
    n2 <- ncol(data2)

    # Compute t-statistic at each time point
    m1 <- rowMeans(data1, na.rm = TRUE)
    m2 <- rowMeans(data2, na.rm = TRUE)
    v1 <- apply(data1, 1, var, na.rm = TRUE)
    v2 <- apply(data2, 1, var, na.rm = TRUE)

    # Pooled variance (equal variance assumption for SPM)
    sp <- sqrt(((n1 - 1) * v1 + (n2 - 1) * v2) / (n1 + n2 - 2))
    se <- sp * sqrt(1/n1 + 1/n2)

    t_stat <- (m1 - m2) / se
    df <- n1 + n2 - 2
    # model residuals: each observation minus its group mean
    residuals <- cbind(data1 - m1, data2 - m2)

  } else {
    # One-sample t-test against zero
    data1 <- data[, group1, drop = FALSE]
    n1 <- ncol(data1)

    m <- rowMeans(data1, na.rm = TRUE)
    s <- apply(data1, 1, sd, na.rm = TRUE)
    se <- s / sqrt(n1)

    t_stat <- m / se
    df <- n1 - 1
    residuals <- data1 - m
  }

  # Handle NA/Inf
  t_stat[!is.finite(t_stat)] <- 0

  # Estimate smoothness (FWHM) from the residual field (not the statistic)
  fwhm <- .estimateFWHMresiduals(t(residuals))
  resel_count <- .reselCount(n_time, fwhm)

  # RFT critical height via the expected Euler characteristic
  threshold <- .rftCriticalT(alpha, df, resel_count, two_tailed)

  # Find significant clusters
  clusters <- .findSPMClusters(t_stat, threshold, two_tailed)

  # Corrected cluster-extent p-values (Worsley/Friston RFT at the cluster-defining
  # threshold). The .rftClusterPValue formula now matches rft1d.p_cluster to
  # machine precision (VAL-07). spm1d's reported cluster P uses an idiosyncratic
  # internal cluster height rather than the threshold, so full-pipeline P differs.
  if (length(clusters) > 0) {
    ec <- .ecDensityT(threshold, df)
    for (i in seq_along(clusters)) {
      clusters[[i]]$p_cluster <- .rftClusterPValue(
        clusters[[i]]$extent, resel_count, fwhm, ec["rho0"], ec["rho1"],
        tails = if (two_tailed) 2 else 1
      )
    }
  }

  # Pointwise p-values (uncorrected)
  if (two_tailed) {
    p_values <- 2 * stats::pt(-abs(t_stat), df)
  } else {
    p_values <- stats::pt(-t_stat, df)
  }

  result <- list(
    t = t_stat,
    threshold = threshold,
    clusters = clusters,
    p_values = p_values,
    df = df,
    fwhm = fwhm,
    resel_count = resel_count,
    alpha = alpha,
    two_tailed = two_tailed,
    n_time = n_time,
    test_type = if (two_sample) "two-sample" else "one-sample"
  )

  class(result) <- c("spm_result", "list")
  result
}

#' SPM paired t-test
#'
#' Performs SPM analysis for paired/repeated measures data.
#'
#' @param x A PhysioExperiment object or matrix (time x observations).
#' @param condition1 Indices for first condition.
#' @param condition2 Indices for second condition.
#' @param alpha Significance level.
#' @param two_tailed Logical; if TRUE, performs two-tailed test.
#'
#' @return A list of class "spm_result".
#' @export
#' @examples
#' # Pre-post intervention comparison
#' set.seed(123)
#' pre <- matrix(rnorm(100 * 15), nrow = 100)
#' post <- pre + 0.8  # Effect across all time points
#' post[30:50, ] <- post[30:50, ] + 0.5  # Additional effect
#'
#' data <- cbind(pre, post)
#' pe <- PhysioExperiment(assays = list(values = data), samplingRate = 100)
#'
#' result <- spmPairedTTest(pe, condition1 = 1:15, condition2 = 16:30)
spmPairedTTest <- function(x, condition1, condition2,
                            alpha = 0.05, two_tailed = TRUE) {

  # Extract data
  if (inherits(x, "PhysioExperiment")) {
    assay_name <- defaultAssay(x)
    data <- SummarizedExperiment::assay(x, assay_name)
  } else if (is.matrix(x)) {
    data <- x
  } else {
    stop("Input must be a PhysioExperiment or matrix", call. = FALSE)
  }

  data1 <- data[, condition1, drop = FALSE]
  data2 <- data[, condition2, drop = FALSE]

  if (ncol(data1) != ncol(data2)) {
    stop("Paired test requires equal number of observations in each condition",
         call. = FALSE)
  }

  # Compute differences
  diffs <- data1 - data2

  # One-sample t-test on differences
  result <- spmTTest(diffs, alpha = alpha, two_tailed = two_tailed)
  result$test_type <- "paired"

  result
}

#' SPM ANOVA (F-test)
#'
#' Performs SPM F-test for comparing multiple groups.
#'
#' @param x A PhysioExperiment object or matrix (time x observations).
#' @param groups Factor or list indicating group membership.
#' @param alpha Significance level.
#'
#' @return A list of class "spm_result" containing F-statistics and clusters.
#' @export
#' @examples
#' # Three-group comparison
#' set.seed(123)
#' g1 <- matrix(rnorm(100 * 8), nrow = 100)
#' g2 <- matrix(rnorm(100 * 8), nrow = 100) + 0.5
#' g3 <- matrix(rnorm(100 * 8), nrow = 100) + 1.0
#'
#' data <- cbind(g1, g2, g3)
#' groups <- factor(rep(c("A", "B", "C"), each = 8))
#'
#' pe <- PhysioExperiment(assays = list(values = data), samplingRate = 100)
#' result <- spmAnova(pe, groups = groups)
spmAnova <- function(x, groups, alpha = 0.05) {

  # Extract data
  if (inherits(x, "PhysioExperiment")) {
    assay_name <- defaultAssay(x)
    data <- SummarizedExperiment::assay(x, assay_name)
  } else if (is.matrix(x)) {
    data <- x
  } else {
    stop("Input must be a PhysioExperiment or matrix", call. = FALSE)
  }

  n_time <- nrow(data)
  n_obs <- ncol(data)

  groups <- droplevels(as.factor(groups))           # drop unused/empty levels
  group_levels <- levels(groups)
  n_groups <- length(group_levels)

  if (n_groups < 2) {
    stop("At least 2 groups required for ANOVA", call. = FALSE)
  }
  df1 <- n_groups - 1
  df2 <- n_obs - n_groups
  if (df2 < 1L) {
    stop("ANOVA requires more observations than groups (residual df >= 1).",
         call. = FALSE)
  }

  # Compute F-statistic at each time point
  f_stat <- numeric(n_time)
  residuals <- matrix(0, n_time, n_obs)             # data minus group means

  for (t_idx in seq_len(n_time)) {
    vals <- data[t_idx, ]

    # Between-group variance
    grand_mean <- mean(vals, na.rm = TRUE)
    ss_between <- 0
    ss_within <- 0

    for (g in group_levels) {
      gi <- which(groups == g)
      group_vals <- vals[gi]
      n_g <- sum(!is.na(group_vals))
      group_mean <- mean(group_vals, na.rm = TRUE)

      ss_between <- ss_between + n_g * (group_mean - grand_mean)^2
      ss_within <- ss_within + sum((group_vals - group_mean)^2, na.rm = TRUE)
      residuals[t_idx, gi] <- group_vals - group_mean
    }

    ms_between <- ss_between / df1
    ms_within <- ss_within / df2

    if (ms_within > 0) {
      f_stat[t_idx] <- ms_between / ms_within
    }
  }

  f_stat[!is.finite(f_stat)] <- 0                   # guard degenerate nodes

  # Estimate smoothness from the residual field, then the RFT F-threshold
  fwhm <- .estimateFWHMresiduals(t(residuals))
  resel_count <- .reselCount(n_time, fwhm)
  threshold <- .rftCriticalF(alpha, df1, df2, resel_count)

  # Find clusters + corrected cluster-extent p-values
  clusters <- .findSPMClustersF(f_stat, threshold)
  # Corrected cluster-extent p-values (RFT at the threshold); the formula matches
  # rft1d.p_cluster to machine precision (VAL-07).
  if (length(clusters) > 0) {
    ec <- .ecDensityF(threshold, df1, df2)
    for (i in seq_along(clusters)) {
      clusters[[i]]$p_cluster <- .rftClusterPValue(
        clusters[[i]]$extent, resel_count, fwhm, ec["rho0"], ec["rho1"], tails = 1
      )
    }
  }

  # Pointwise p-values
  p_values <- stats::pf(f_stat, df1, df2, lower.tail = FALSE)

  result <- list(
    f = f_stat,
    threshold = threshold,
    clusters = clusters,
    p_values = p_values,
    df1 = df1,
    df2 = df2,
    fwhm = fwhm,
    resel_count = resel_count,
    alpha = alpha,
    n_time = n_time,
    groups = group_levels,
    test_type = "anova"
  )

  class(result) <- c("spm_result", "list")
  result
}

#' Estimate FWHM (Full Width at Half Maximum) for RFT
#' @noRd
.estimateFWHM <- function(x) {
  # Legacy statistic-based estimator; superseded by .estimateFWHMresiduals().
  # Retained only as a fallback when residuals are unavailable.
  n <- length(x)
  if (n < 3) return(1)
  var_dx <- var(diff(x), na.rm = TRUE)
  var_x <- var(x, na.rm = TRUE)
  if (var_x == 0 || !is.finite(var_dx) || var_dx == 0) return(n)
  max(1, min(sqrt(4 * log(2) * var_x / var_dx), n))
}

#' Estimate field smoothness (FWHM) from standardized model residuals
#'
#' Residual-based smoothness estimator (Kiebel et al. 1999): for each node the
#' normalized variance of the residual field's gradient gives the resels per
#' node; the FWHM is the reciprocal of their mean. This is the smoothness that
#' enters RFT, and must be estimated from the model RESIDUALS (data minus group
#' means), not from the test statistic itself.
#' @param residuals An `n_obs x n_time` matrix of model residuals.
#' @noRd
.estimateFWHMresiduals <- function(residuals) {
  R <- as.matrix(residuals)
  R[!is.finite(R)] <- 0
  n_obs <- nrow(R); Q <- ncol(R)
  if (n_obs < 2L || Q < 3L) return(Q)
  ssq <- colSums(R^2)                               # per-node sum of squares
  grad <- matrix(0, n_obs, Q)                       # node-wise gradient per obs
  grad[, 1L] <- R[, 2L] - R[, 1L]
  grad[, Q] <- R[, Q] - R[, Q - 1L]
  if (Q > 2L) grad[, 2:(Q - 1L)] <- (R[, 3:Q] - R[, 1:(Q - 2L)]) / 2
  v <- colSums(grad^2)
  keep <- ssq > 0
  if (!any(keep)) return(Q)
  lambda <- v[keep] / ssq[keep]                     # normalized gradient var
  rpn <- sqrt(lambda / (4 * log(2)))                # resels per node
  rpn <- rpn[is.finite(rpn) & rpn > 0]
  if (!length(rpn)) return(Q)
  max(1, min(1 / mean(rpn), Q))
}

#' Euler-characteristic densities of a 1D Student-t field (Worsley 1996)
#' @noRd
.ecDensityT <- function(u, df) {
  rho0 <- stats::pt(u, df, lower.tail = FALSE)
  rho1 <- sqrt(4 * log(2)) / (2 * pi) * (1 + u^2 / df)^(-(df - 1) / 2)
  c(rho0 = rho0, rho1 = rho1)
}

#' Euler-characteristic densities of a 1D F field (Worsley 1996)
#'
#' The 1D (dimension-1) F-field EC density is bracket-free -- the linear-in-`cc`
#' factor belongs to the 2D density. At `df1 = 1` this reduces to twice the
#' one-sided t-field density (the F(1, df2) = t(df2)^2 identity), which pins the
#' normalising constant to `sqrt(4 log 2) / sqrt(2 pi)`.
#' @noRd
.ecDensityF <- function(u, df1, df2) {
  rho0 <- stats::pf(u, df1, df2, lower.tail = FALSE)
  cc <- df1 * u / df2
  logk <- lgamma((df1 + df2 - 1) / 2) + 0.5 * log(2) -
    lgamma(df1 / 2) - lgamma(df2 / 2)
  rho1 <- sqrt(4 * log(2)) / sqrt(2 * pi) * exp(logk) *
    cc^((df1 - 1) / 2) * (1 + cc)^(-(df1 + df2 - 2) / 2)
  c(rho0 = rho0, rho1 = rho1)
}

#' RFT critical threshold for a 1D t-field
#'
#' Solves the expected-Euler-characteristic equation
#' `E[EC(u)] = rho0(u) + R1 * rho1(u) = target` for the critical height `u`,
#' where `R1 = resel_count` and `target = alpha` (one-tailed) or `alpha/2`
#' (two-tailed, by field symmetry). This is the RFT survival threshold used by
#' spm1d / rft1d, replacing the previous resel-Bonferroni approximation.
#' @noRd
.rftCriticalT <- function(alpha, df, resel_count, two_tailed = TRUE) {
  target <- if (two_tailed) alpha / 2 else alpha
  u0 <- stats::qt(1 - target, df)                   # single-comparison height
  if (!is.finite(resel_count) || resel_count <= .Machine$double.eps) return(u0)
  f <- function(u) unname(.ecDensityT(u, df)["rho0"] +
                            resel_count * .ecDensityT(u, df)["rho1"] - target)
  .rftSolve(f, u0, multiplicative = FALSE)
}

#' RFT critical threshold for a 1D F-field
#' @noRd
.rftCriticalF <- function(alpha, df1, df2, resel_count) {
  u0 <- stats::qf(1 - alpha, df1, df2)
  if (!is.finite(resel_count) || resel_count <= .Machine$double.eps) return(u0)
  f <- function(u) unname(.ecDensityF(u, df1, df2)["rho0"] +
                            resel_count * .ecDensityF(u, df1, df2)["rho1"] - alpha)
  .rftSolve(f, u0, multiplicative = TRUE)
}

#' Solve E[EC](u) = target (via uniroot) for the RFT critical height
#'
#' Robust to fields whose 1D EC density does not decay in `u` -- a Student-t
#' field with `df = 1` or an F-field with `df2 <= 3` -- for which the excursion
#' probability never falls to the target and no finite threshold exists; `Inf`
#' is returned in that case (nothing can reach significance) rather than letting
#' `uniroot` fail on a non-bracketing interval.
#' @noRd
.rftSolve <- function(f, u0, multiplicative = FALSE, max_it = 500L) {
  if (f(u0) <= 0) return(u0)                         # correction negligible
  hi <- u0 + 1; it <- 0L
  while (f(hi) > 0 && it < max_it) {
    hi <- if (multiplicative) hi * 1.5 else hi + 1
    it <- it + 1L
  }
  if (f(hi) > 0) return(Inf)                         # no root: EC never decays
  stats::uniroot(f, c(u0, hi), tol = 1e-9)$root
}

#' Number of resolution elements (resels) of a 1D field
#' @noRd
.reselCount <- function(n_time, fwhm) (n_time - 1) / fwhm

#' Find SPM clusters for t-test
#' @noRd
.findSPMClusters <- function(t_stat, threshold, two_tailed) {
  n <- length(t_stat)
  clusters <- list()

  if (two_tailed) {
    # Positive clusters
    pos_sig <- t_stat > threshold
    pos_clusters <- .extractClusters(pos_sig, t_stat, "positive")
    clusters <- c(clusters, pos_clusters)

    # Negative clusters
    neg_sig <- t_stat < -threshold
    neg_clusters <- .extractClusters(neg_sig, t_stat, "negative")
    clusters <- c(clusters, neg_clusters)
  } else {
    sig <- t_stat > threshold
    clusters <- .extractClusters(sig, t_stat, "positive")
  }

  clusters
}

#' Find SPM clusters for F-test
#' @noRd
.findSPMClustersF <- function(f_stat, threshold) {
  sig <- f_stat > threshold
  .extractClusters(sig, f_stat, "positive")
}

#' Extract clusters from significance mask
#' @noRd
.extractClusters <- function(sig_mask, stat_values, direction) {
  if (!any(sig_mask)) return(list())

  rle_sig <- rle(sig_mask)
  clusters <- list()

  pos <- 1
  for (i in seq_along(rle_sig$lengths)) {
    if (rle_sig$values[i]) {
      start <- pos
      end <- pos + rle_sig$lengths[i] - 1
      extent <- end - start + 1

      cluster_stat <- sum(abs(stat_values[start:end]))
      peak_stat <- if (direction == "positive") {
        max(stat_values[start:end])
      } else {
        min(stat_values[start:end])
      }

      clusters <- c(clusters, list(list(
        start = start,
        end = end,
        extent = extent,
        cluster_stat = cluster_stat,
        peak_stat = peak_stat,
        direction = direction
      )))
    }
    pos <- pos + rle_sig$lengths[i]
  }

  clusters
}

#' RFT cluster-extent p-value (Friston 1994 set/cluster-level inference)
#'
#' Corrected probability of a supra-threshold cluster of at least `extent`
#' nodes, for a 1D field with the given resel count and EC densities at the
#' cluster-defining height. Uses the expected number of clusters `E_m` and the
#' expected supra-threshold resel count `E_n` to set the cluster-size scale.
#' @param extent Cluster extent in nodes.
#' @param resel_count Field resel count `R1`.
#' @param fwhm Field FWHM in nodes.
#' @param rho0,rho1 EC densities at the cluster-defining threshold.
#' @param tails 2 for a two-tailed t-field (both signs), 1 otherwise.
#' @noRd
.rftClusterPValue <- function(extent, resel_count, fwhm, rho0, rho1, tails = 1) {
  if (!is.finite(rho0) || !is.finite(rho1) || rho0 <= 0 || rho1 <= 0) return(1)
  k <- extent / fwhm                              # cluster extent in resels
  # Expected number of clusters (per tail) and expected cluster extent in resels
  # follow Friston et al. 1994 / Worsley et al. 1996 as implemented by spm1d/rft1d
  # (spm12 convention). The expected cluster extent is E[k] = rho0 / rho1 (the
  # resel count cancels); an earlier form used (E[m]/E[n])^2 for beta, which adds
  # a spurious 1/resel term and diverges from spm1d/rft1d (VAL-07).
  ec_clusters <- resel_count * rho1 + rho0        # E[m], expected clusters (one tail)
  ek <- rho0 / rho1                               # E[k], expected cluster extent (resels)
  if (!is.finite(ek) || ek <= 0) return(1)
  beta <- (gamma(1.5) / ek)^2                     # 1D (D = 1) cluster-size scale
  p_one <- exp(-beta * k^2)                       # P(a single cluster's extent >= k)
  p1 <- 1 - exp(-ec_clusters * p_one)             # one-tailed cluster probability
  min(1, tails * p1)                              # spm1d doubles for two-tailed
}

#' Plot SPM result
#'
#' Visualizes SPM analysis results with significance thresholds.
#'
#' @param x An spm_result object.
#' @param time_axis Optional time axis values.
#' @param show_threshold Logical; if TRUE, shows significance threshold.
#' @param show_clusters Logical; if TRUE, highlights significant clusters.
#' @param title Plot title.
#'
#' @return A ggplot object.
#' @export
#' @examples
#' # Create and plot SPM result
#' set.seed(123)
#' g1 <- matrix(rnorm(100 * 10), nrow = 100)
#' g2 <- matrix(rnorm(100 * 10), nrow = 100)
#' g2[40:60, ] <- g2[40:60, ] + 1.5
#'
#' pe <- PhysioExperiment(
#'   assays = list(values = cbind(g1, g2)),
#'   samplingRate = 100
#' )
#' result <- spmTTest(pe, group1 = 1:10, group2 = 11:20)
#' plotSPM(result)
plotSPM <- function(x, time_axis = NULL, show_threshold = TRUE,
                    show_clusters = TRUE, title = NULL) {

  if (!inherits(x, "spm_result")) {
    stop("Input must be an spm_result object", call. = FALSE)
  }

  n_time <- x$n_time

  if (is.null(time_axis)) {
    time_axis <- seq(0, 100, length.out = n_time)
  }

  # Determine statistic type
  if (!is.null(x$t)) {
    stat <- x$t
    stat_label <- "SPM{t}"
  } else if (!is.null(x$f)) {
    stat <- x$f
    stat_label <- "SPM{F}"
  } else if (!is.null(x$T2)) {
    stat <- x$T2
    stat_label <- "SPM{T2}"
  } else if (!is.null(x$X2)) {
    stat <- x$X2
    stat_label <- "SPM{X2}"
  } else {
    stop("No statistic found in result", call. = FALSE)
  }

  # Create data frame
  df <- data.frame(
    time = time_axis,
    stat = stat
  )

  # Base plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$time, y = .data$stat)) +
    ggplot2::geom_line(linewidth = 1, color = "black") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    ggplot2::labs(x = "Time (%)", y = stat_label) +
    ggplot2::theme_minimal()

  # Add threshold lines
  if (show_threshold) {
    p <- p +
      ggplot2::geom_hline(yintercept = x$threshold,
                          linetype = "dashed", color = "red") +
      ggplot2::annotate("text", x = max(time_axis) * 0.95, y = x$threshold,
                        label = sprintf("p = %.3f", x$alpha),
                        color = "red", vjust = -0.5, size = 3)

    if (!is.null(x$t) && x$two_tailed) {
      p <- p +
        ggplot2::geom_hline(yintercept = -x$threshold,
                            linetype = "dashed", color = "red")
    }
  }

  # Highlight significant clusters
  if (show_clusters && length(x$clusters) > 0) {
    for (cl in x$clusters) {
      start_time <- time_axis[cl$start]
      end_time <- time_axis[cl$end]

      p <- p +
        ggplot2::annotate("rect",
                          xmin = start_time, xmax = end_time,
                          ymin = -Inf, ymax = Inf,
                          fill = "red", alpha = 0.2)
    }
  }

  # Add title
  if (!is.null(title)) {
    p <- p + ggplot2::ggtitle(title)
  } else {
    test_desc <- switch(x$test_type,
      "two-sample" = "Two-sample t-test",
      "one-sample" = "One-sample t-test",
      "paired" = "Paired t-test",
      "anova" = sprintf("ANOVA (%d groups)", length(x$groups)),
      "regression" = "Linear regression",
      "manova" = sprintf("MANOVA (%d groups, %d components)",
                         length(x$groups), x$n_comp),
      "snpm" = sprintf("Permutation SPM{%s}", x$statistic),
      "SPM Analysis"
    )
    p <- p + ggplot2::ggtitle(test_desc)
  }

  p
}

#' Print SPM result
#' @param x An `spm_result` object from [spmTTest()], [spmPairedTTest()] or
#'   [spmAnova()].
#' @param ... Ignored; present for S3 `print` method consistency.
#' @return `x`, invisibly.
#' @keywords internal
#' @export
print.spm_result <- function(x, ...) {
  cat("SPM Analysis Result\n")
  cat("==================\n")
  cat(sprintf("Test type: %s\n", x$test_type))
  cat(sprintf("Time points: %d\n", x$n_time))
  cat(sprintf("Alpha: %.3f\n", x$alpha))
  cat(sprintf("Threshold: %.3f\n", x$threshold))
  if (!is.null(x$fwhm)) {
    cat(sprintf("FWHM (smoothness): %.2f\n", x$fwhm))
  }
  if (!is.null(x$resel_count)) {
    cat(sprintf("Resel count: %.2f\n", x$resel_count))
  }
  if (!is.null(x$n_permutations)) {
    cat(sprintf("Permutations: %d\n", x$n_permutations))
  }

  if (length(x$clusters) > 0) {
    cat(sprintf("\nSignificant clusters: %d\n", length(x$clusters)))
    for (i in seq_along(x$clusters)) {
      cl <- x$clusters[[i]]
      cat(sprintf("  Cluster %d: [%d-%d] extent=%d, p=%.4f\n",
                  i, cl$start, cl$end, cl$extent,
                  if (!is.null(cl$p_cluster)) cl$p_cluster else NA))
    }
  } else {
    cat("\nNo significant clusters found\n")
  }

  invisible(x)
}
