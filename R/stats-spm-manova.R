# Vector-field (multivariate) SPM: a one-way MANOVA at each node. For two groups
# the field is Hotelling's T^2 (which reduces to the squared two-sample t-field
# when there is a single component); for more than two groups it is the
# chi-square field from the Bartlett-transformed Wilks' Lambda. Reuses the RFT
# machinery in stats-spm.R and adds the chi-square field EC densities.

#' Euler-characteristic densities of a 1D chi-square field (Worsley 1994)
#' @noRd
.ecDensityChi2 <- function(u, df) {
  rho0 <- stats::pchisq(u, df, lower.tail = FALSE)
  logrho1 <- 0.5 * log(4 * log(2)) - 0.5 * log(2 * pi) +
    (df - 1) / 2 * log(u) - u / 2 - ((df - 2) / 2) * log(2) - lgamma(df / 2)
  rho1 <- exp(logrho1)
  c(rho0 = rho0, rho1 = rho1)
}

#' RFT critical threshold for a 1D chi-square field
#' @noRd
.rftCriticalChi2 <- function(alpha, df, resel_count) {
  u0 <- stats::qchisq(1 - alpha, df)
  if (!is.finite(resel_count) || resel_count <= .Machine$double.eps) return(u0)
  f <- function(u) unname(.ecDensityChi2(u, df)["rho0"] +
                            resel_count * .ecDensityChi2(u, df)["rho1"] - alpha)
  .rftSolve(f, u0, multiplicative = TRUE)
}

# Coerce the multivariate input to an (n_time x n_obs x n_comp) array.
.spmVectorArray <- function(x, vector_components = NULL) {
  if (is.array(x) && length(dim(x)) == 3L) return(x)
  if (is.list(x) && !is.data.frame(x)) {
    comps <- lapply(x, as.matrix)
    d <- dim(comps[[1]])
    if (!all(vapply(comps, function(m) identical(dim(m), d), logical(1)))) {
      stop("all component matrices must have the same dimensions.",
           call. = FALSE)
    }
    arr <- array(0, c(d[1], d[2], length(comps)))
    for (k in seq_along(comps)) arr[, , k] <- comps[[k]]
    return(arr)
  }
  if (is.matrix(x)) {
    if (is.null(vector_components) || vector_components == 1L) {
      return(array(x, c(nrow(x), ncol(x), 1L)))
    }
    stop("for a matrix with >1 component, pass a 3D array or a list of ",
         "component matrices.", call. = FALSE)
  }
  stop("Input must be a 3D array (time x obs x component) or a list of ",
       "component matrices.", call. = FALSE)
}

#' SPM one-way MANOVA (Hotelling's T^2 / chi-square vector field)
#'
#' Computes a multivariate SPM field for a one-way MANOVA at every node, for
#' vector-valued waveforms (e.g. 3D joint angles). With two groups the field is
#' Hotelling's \eqn{T^2} (reducing to the squared two-sample t-field for a single
#' component); with more than two groups it is the chi-square field of the
#' Bartlett-transformed Wilks' \eqn{\Lambda}. Field-level significance uses
#' random-field theory (an F-field threshold for \eqn{T^2}, a chi-square-field
#' threshold for the \eqn{X^2} field).
#'
#' @param x A 3D array \code{time x obs x component}, or a list of component
#'   matrices (each \code{time x obs}), or a single \code{time x obs} matrix
#'   (one component).
#' @param groups A grouping factor, one value per observation.
#' @param vector_components Optional integer; only used to disambiguate a 2D
#'   matrix input (defaults to a single component).
#' @param alpha Significance level (default 0.05).
#' @return A list of class \code{"spm_result"} (\code{test_type = "manova"}) with
#'   the statistic field (\code{T2} or \code{X2}), RFT \code{threshold},
#'   significant \code{clusters}, \code{p_values}, degrees of freedom,
#'   \code{fwhm}, and \code{resel_count}.
#' @references Pataky 2016 (vector-field 1D SPM); Worsley 1994 (chi-square RFT).
#'   spm1d.stats.manova1 / hotellings2.
#' @seealso [spmTTest()], [spmRegression()], [spmSnPM()]
#' @export
#' @examples
#' set.seed(1)
#' # two groups, 3 components (e.g. hip flexion/abduction/rotation)
#' arr <- array(rnorm(50 * 16 * 3), c(50, 16, 3))
#' arr[20:30, 1:8, ] <- arr[20:30, 1:8, ] + 1.2
#' spmMANOVA(arr, groups = rep(c("A", "B"), each = 8))
spmMANOVA <- function(x, groups, vector_components = NULL, alpha = 0.05) {
  arr <- .spmVectorArray(x, vector_components)
  n_time <- dim(arr)[1]; n_obs <- dim(arr)[2]; n_comp <- dim(arr)[3]
  groups <- droplevels(as.factor(groups))
  if (length(groups) != n_obs) {
    stop("'groups' must have one value per observation (", n_obs, ").",
         call. = FALSE)
  }
  glevs <- levels(groups); g <- length(glevs)
  if (g < 2L) stop("MANOVA requires at least 2 groups.", call. = FALSE)
  df_err <- n_obs - g
  if (df_err < n_comp + 1L) {
    stop("MANOVA needs residual df (n - groups = ", df_err,
         ") of at least components + 1 (", n_comp + 1L, ").", call. = FALSE)
  }

  stat <- numeric(n_time)
  resid_arr <- array(0, c(n_time, n_obs, n_comp))       # residuals for smoothness
  two_group <- (g == 2L)
  gidx <- lapply(glevs, function(l) which(groups == l))

  for (t_idx in seq_len(n_time)) {
    Y <- matrix(arr[t_idx, , ], n_obs, n_comp)          # obs x component
    grand <- colMeans(Y)
    W <- matrix(0, n_comp, n_comp)                      # within SSCP
    B <- matrix(0, n_comp, n_comp)                      # between SSCP
    res_t <- matrix(0, n_obs, n_comp)
    means <- vector("list", g)
    for (j in seq_len(g)) {
      idx <- gidx[[j]]; ng <- length(idx)
      mj <- colMeans(Y[idx, , drop = FALSE]); means[[j]] <- mj
      dev <- sweep(Y[idx, , drop = FALSE], 2, mj)
      W <- W + crossprod(dev)
      B <- B + ng * tcrossprod(mj - grand)
      res_t[idx, ] <- dev
    }
    resid_arr[t_idx, , ] <- res_t
    if (two_group) {
      # Hotelling's two-sample T^2 = a' Sp^{-1} a scaled by n1 n2 / N
      n1 <- length(gidx[[1]]); n2 <- length(gidx[[2]])
      Sp <- W / df_err                                  # pooled covariance
      d <- means[[1]] - means[[2]]
      Spinv <- tryCatch(solve(Sp), error = function(e) NULL)
      if (!is.null(Spinv)) {
        stat[t_idx] <- (n1 * n2 / (n1 + n2)) *
          as.numeric(t(d) %*% Spinv %*% d)
      }
    } else {
      # Wilks' Lambda -> Bartlett chi-square
      WB <- W + B
      lambda <- tryCatch(det(W) / det(WB), error = function(e) NA_real_)
      if (is.finite(lambda) && lambda > 0) {
        stat[t_idx] <- -(n_obs - 1 - (n_comp + g) / 2) * log(lambda)
      }
    }
  }
  stat[!is.finite(stat) | stat < 0] <- 0

  # pooled residual smoothness: stack every (obs, component) residual series as a
  # row (n_obs * n_comp rows x n_time). For a single component this is exactly
  # the signed-residual field spmTTest uses, so T^2 reduces to the squared t-field.
  R_stack <- matrix(aperm(resid_arr, c(2, 3, 1)), n_obs * n_comp, n_time)
  fwhm <- .estimateFWHMresiduals(R_stack)
  resel_count <- .reselCount(n_time, fwhm)

  if (two_group) {
    # T^2 = (df_err * n_comp / (df_err - n_comp + 1)) * F,  F ~ F(n_comp, df_err - n_comp + 1)
    df1 <- n_comp; df2 <- df_err - n_comp + 1L
    scale <- df_err * n_comp / df2
    f_crit <- .rftCriticalF(alpha, df1, df2, resel_count)
    threshold <- scale * f_crit
    f_field <- stat / scale
    p_values <- stats::pf(f_field, df1, df2, lower.tail = FALSE)
    clusters <- .findSPMClustersF(stat, threshold)
    if (length(clusters) > 0) {
      ec <- .ecDensityF(f_crit, df1, df2)
      for (i in seq_along(clusters)) {
        clusters[[i]]$p_cluster <- .rftClusterPValue(
          clusters[[i]]$extent, resel_count, fwhm, ec["rho0"], ec["rho1"], 1)
      }
    }
    df_out <- c(df1 = df1, df2 = df2)
    stat_name <- "T2"
  } else {
    df_chi <- n_comp * (g - 1L)
    threshold <- .rftCriticalChi2(alpha, df_chi, resel_count)
    p_values <- stats::pchisq(stat, df_chi, lower.tail = FALSE)
    clusters <- .findSPMClustersF(stat, threshold)
    if (length(clusters) > 0) {
      ec <- .ecDensityChi2(threshold, df_chi)
      for (i in seq_along(clusters)) {
        clusters[[i]]$p_cluster <- .rftClusterPValue(
          clusters[[i]]$extent, resel_count, fwhm, ec["rho0"], ec["rho1"], 1)
      }
    }
    df_out <- c(df = df_chi)
    stat_name <- "X2"
  }

  res <- list(
    threshold = threshold, clusters = clusters, p_values = p_values,
    df = df_out, fwhm = fwhm, resel_count = resel_count, alpha = alpha,
    n_time = n_time, n_comp = n_comp, groups = glevs,
    stat_name = stat_name, test_type = "manova")
  res[[stat_name]] <- stat
  structure(res, class = c("spm_result", "list"))
}
