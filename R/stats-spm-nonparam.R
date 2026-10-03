# Non-parametric SPM (SnPM): the field-maximum permutation distribution gives a
# distribution-free critical threshold and FWER control (Nichols & Holmes 2002),
# a validity cross-check on the parametric RFT thresholds in stats-spm.R. Supports
# t (one/two-sample), F (one-way ANOVA), and T^2 (two-group MANOVA) statistics.

# --- fast node-wise statistic fields (vectorised where possible) ---

.snpm_t_2samp <- function(data, g1, g2) {
  d1 <- data[, g1, drop = FALSE]; d2 <- data[, g2, drop = FALSE]
  n1 <- length(g1); n2 <- length(g2)
  m1 <- rowMeans(d1); m2 <- rowMeans(d2)
  v1 <- rowSums((d1 - m1)^2) / (n1 - 1); v2 <- rowSums((d2 - m2)^2) / (n2 - 1)
  sp <- sqrt(((n1 - 1) * v1 + (n2 - 1) * v2) / (n1 + n2 - 2))
  t <- (m1 - m2) / (sp * sqrt(1 / n1 + 1 / n2))
  t[!is.finite(t)] <- 0; t
}

.snpm_t_1samp <- function(data, signs) {
  d <- sweep(data, 2, signs, `*`)
  n <- ncol(d); m <- rowMeans(d)
  s <- sqrt(rowSums((d - m)^2) / (n - 1))
  t <- m / (s / sqrt(n)); t[!is.finite(t)] <- 0; t
}

.snpm_f <- function(data, labels, glevs) {
  n_obs <- ncol(data); k <- length(glevs)
  grand <- rowMeans(data)
  ssb <- numeric(nrow(data)); ssw <- numeric(nrow(data))
  for (l in glevs) {
    idx <- which(labels == l); ng <- length(idx)
    mj <- rowMeans(data[, idx, drop = FALSE])
    ssb <- ssb + ng * (mj - grand)^2
    ssw <- ssw + rowSums((data[, idx, drop = FALSE] - mj)^2)
  }
  fstat <- (ssb / (k - 1)) / (ssw / (n_obs - k))
  fstat[!is.finite(fstat)] <- 0; fstat
}

.snpm_t2 <- function(arr, labels, glevs) {
  n_time <- dim(arr)[1]; n_obs <- dim(arr)[2]; n_comp <- dim(arr)[3]
  i1 <- which(labels == glevs[1]); i2 <- which(labels == glevs[2])
  n1 <- length(i1); n2 <- length(i2); df_err <- n_obs - 2
  out <- numeric(n_time)
  for (t_idx in seq_len(n_time)) {
    Y <- matrix(arr[t_idx, , ], n_obs, n_comp)
    m1 <- colMeans(Y[i1, , drop = FALSE]); m2 <- colMeans(Y[i2, , drop = FALSE])
    W <- crossprod(sweep(Y[i1, , drop = FALSE], 2, m1)) +
      crossprod(sweep(Y[i2, , drop = FALSE], 2, m2))
    d <- m1 - m2
    Sinv <- tryCatch(solve(W / df_err), error = function(e) NULL)
    if (!is.null(Sinv)) {
      out[t_idx] <- (n1 * n2 / n_obs) * as.numeric(t(d) %*% Sinv %*% d)
    }
  }
  out
}

#' Non-parametric (permutation) SPM
#'
#' Builds the permutation distribution of the field-maximum statistic to obtain a
#' distribution-free critical threshold with strong FWER control (Nichols &
#' Holmes 2002), plus permutation cluster- and set-level p-values. This is the
#' validity cross-check for the parametric RFT thresholds: as the field smoothness
#' grows the permutation threshold converges to the RFT threshold.
#'
#' @param x A \code{PhysioExperiment}/matrix (time x obs) for \code{t}/\code{F},
#'   or a 3D array / list of component matrices for \code{T2}.
#' @param groups Grouping factor (one per observation) for \code{F} or \code{T2},
#'   or a two-sample \code{t}.
#' @param group1,group2 Observation indices for a two-sample \code{t}
#'   (\code{group2 = NULL} gives a one-sample sign-flip test).
#' @param statistic \code{"t"}, \code{"F"}, or \code{"T2"}.
#' @param n_permutations Number of random permutations (default 1000).
#' @param alpha Significance level (default 0.05).
#' @param two_tailed For \code{t}: use the two-tailed field maximum \code{max|t|}.
#' @param vector_components Passed to the \code{T2} array coercion.
#' @param seed Optional RNG seed for reproducibility.
#' @return A list of class \code{"spm_result"} (\code{test_type = "snpm"}) with
#'   the observed statistic field, the permutation \code{threshold},
#'   \code{clusters} with permutation p-values, pointwise permutation
#'   \code{p_values}, and the \code{perm_max} distribution.
#' @references Nichols & Holmes 2002; Pataky 2016. spm1d.stats.nonparam.
#' @seealso [spmTTest()], [spmAnova()], [spmMANOVA()]
#' @export
#' @examples
#' set.seed(1)
#' g1 <- matrix(rnorm(50 * 10), 50); g2 <- matrix(rnorm(50 * 10), 50)
#' g2[20:30, ] <- g2[20:30, ] + 1.5
#' spmSnPM(cbind(g1, g2), group1 = 1:10, group2 = 11:20,
#'         n_permutations = 200, seed = 1)
spmSnPM <- function(x, groups = NULL, group1 = NULL, group2 = NULL,
                    statistic = c("t", "F", "T2"), n_permutations = 1000,
                    alpha = 0.05, two_tailed = TRUE, vector_components = NULL,
                    seed = NULL) {
  statistic <- match.arg(statistic)
  if (!is.null(seed)) set.seed(seed)

  # --- observed statistic field + a permutation generator ---
  if (statistic == "t") {
    if (inherits(x, "PhysioExperiment")) {
      data <- SummarizedExperiment::assay(x, defaultAssay(x))
    } else if (is.matrix(x)) data <- x else
      stop("t/F SPM needs a PhysioExperiment or matrix.", call. = FALSE)
    n_obs <- ncol(data); n_time <- nrow(data)
    one_sample <- is.null(group2)
    if (one_sample) {
      g1 <- if (is.null(group1)) seq_len(n_obs) else group1
      obs <- .snpm_t_1samp(data[, g1, drop = FALSE], rep(1, length(g1)))
      perm_fun <- function() .snpm_t_1samp(
        data[, g1, drop = FALSE],
        sample(c(-1, 1), length(g1), replace = TRUE))
      df <- length(g1) - 1
    } else {
      obs <- .snpm_t_2samp(data, group1, group2)
      pool <- c(group1, group2); n1 <- length(group1)
      perm_fun <- function() {
        p <- sample(pool); .snpm_t_2samp(data, p[seq_len(n1)],
                                         p[(n1 + 1):length(pool)])
      }
      df <- length(group1) + length(group2) - 2
    }
    field_max <- function(v) if (two_tailed) max(abs(v)) else max(v)
  } else if (statistic == "F") {
    if (inherits(x, "PhysioExperiment")) {
      data <- SummarizedExperiment::assay(x, defaultAssay(x))
    } else if (is.matrix(x)) data <- x else
      stop("F SPM needs a PhysioExperiment or matrix.", call. = FALSE)
    labels <- droplevels(as.factor(groups)); glevs <- levels(labels)
    n_time <- nrow(data)
    obs <- .snpm_f(data, labels, glevs)
    perm_fun <- function() .snpm_f(data, sample(labels), glevs)
    df <- c(length(glevs) - 1L, ncol(data) - length(glevs))
    two_tailed <- FALSE
    field_max <- function(v) max(v)
  } else {                                               # T2
    arr <- .spmVectorArray(x, vector_components)
    labels <- droplevels(as.factor(groups)); glevs <- levels(labels)
    if (length(glevs) != 2L) {
      stop("statistic = 'T2' supports two groups.", call. = FALSE)
    }
    n_time <- dim(arr)[1]
    obs <- .snpm_t2(arr, labels, glevs)
    perm_fun <- function() .snpm_t2(arr, sample(labels), glevs)
    df <- c(dim(arr)[3], dim(arr)[2] - dim(arr)[3] - 1L)
    two_tailed <- FALSE
    field_max <- function(v) max(v)
  }

  # --- permutation field-maximum distribution (identity permutation included) ---
  n_permutations <- as.integer(n_permutations)
  if (n_permutations < 1L) stop("'n_permutations' must be >= 1.", call. = FALSE)
  perm_max <- numeric(n_permutations)
  perm_max[1] <- field_max(obs)
  # seq_len(n)[-1] is empty at n = 1 (avoids the 2:1 = c(2,1) antipattern)
  for (b in seq_len(n_permutations)[-1L]) perm_max[b] <- field_max(perm_fun())

  # critical threshold: the (1 - alpha) permutation quantile of the field max
  threshold <- as.numeric(stats::quantile(perm_max, 1 - alpha, type = 1,
                                          names = FALSE))
  # pointwise permutation p-values (proportion of perm maxima >= observed)
  stat_for_p <- if (statistic == "t" && two_tailed) abs(obs) else obs
  p_values <- vapply(stat_for_p, function(s) mean(perm_max >= s), numeric(1))

  if (statistic == "t" && two_tailed) {
    clusters <- .findSPMClusters(obs, threshold, TRUE)
  } else {
    clusters <- .findSPMClustersF(obs, threshold)
  }
  for (i in seq_along(clusters)) {
    clusters[[i]]$p_cluster <- mean(perm_max >= abs(clusters[[i]]$peak_stat))
  }

  fname <- switch(statistic, t = "t", F = "f", T2 = "T2")
  res <- list(
    threshold = threshold, clusters = clusters, p_values = p_values,
    df = df, alpha = alpha, two_tailed = two_tailed, n_time = n_time,
    n_permutations = n_permutations, perm_max = perm_max,
    statistic = statistic, test_type = "snpm")
  res[[fname]] <- obs
  structure(res, class = c("spm_result", "list"))
}
