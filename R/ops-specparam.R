utils::globalVariables(c("frequency", "power", "component"))

# Spectral parameterization (specparam / FOOOF; Donoghue et al. 2020): separate
# a power spectrum into an aperiodic 1/f component (offset, exponent, optional
# knee) and a set of Gaussian periodic peaks, on the Welch PSD produced by the
# ops-timefreq machinery.

## Aperiodic model in log10 power.
.specparam_gen_ap <- function(freqs, params, mode) {
  if (mode == "fixed") {
    params[1] - params[2] * log10(freqs)
  } else {
    params[1] - log10(params[2] + freqs^params[3])   # offset, knee, exponent
  }
}

## A single Gaussian peak (in log10 power above the aperiodic component).
.specparam_gaussian <- function(freqs, ctr, hgt, wid) {
  hgt * exp(-(freqs - ctr)^2 / (2 * wid^2))
}

## Least-squares aperiodic fit: closed form for "fixed", L-BFGS-B for "knee".
.specparam_simple_ap <- function(freqs, lp, mode) {
  lf <- log10(freqs)
  fit0 <- stats::lm(lp ~ lf)
  off0 <- unname(stats::coef(fit0)[1])
  exp0 <- unname(-stats::coef(fit0)[2])
  if (mode == "fixed") return(c(offset = off0, exponent = exp0))
  obj <- function(p) sum((lp - (p[1] - log10(p[2] + freqs^p[3])))^2)
  opt <- tryCatch(
    stats::optim(c(off0, 1, max(0.1, exp0)), obj, method = "L-BFGS-B",
                 lower = c(-Inf, 1e-8, 0), upper = c(Inf, 1e4, 8)),
    error = function(e) list(par = c(off0, 1e-8, max(0.1, exp0))))
  stats::setNames(opt$par, c("offset", "knee", "exponent"))
}

## Robust aperiodic fit: fit, drop points sitting above the fit (peaks), refit.
.specparam_robust_ap <- function(freqs, lp, mode) {
  p0 <- .specparam_simple_ap(freqs, lp, mode)
  flat <- lp - .specparam_gen_ap(freqs, p0, mode)
  flat[flat < 0] <- 0
  thr <- stats::quantile(flat, 0.025)
  mask <- flat <= thr
  if (sum(mask) < (if (mode == "fixed") 2L else 3L)) mask <- rep(TRUE, length(freqs))
  .specparam_simple_ap(freqs[mask], lp[mask], mode)
}

## Iterative Gaussian peak fitting on the flattened (aperiodic-removed) spectrum,
## followed by a joint refit. Returns a matrix of (center, height, std).
.specparam_fit_peaks <- function(freqs, flat, max_n_peaks, min_peak_height,
                                 peak_threshold, pw_limits) {
  bw_min <- pw_limits[1] / 2; bw_max <- pw_limits[2] / 2
  fstd <- stats::sd(flat); fl <- flat; guesses <- list()
  for (i in seq_len(max_n_peaks)) {
    mi <- which.max(fl); mh <- fl[mi]
    if (!is.finite(mh) || mh <= peak_threshold * fstd || mh <= min_peak_height) break
    cf <- freqs[mi]; half <- mh / 2
    li <- mi; while (li > 1 && fl[li] > half) li <- li - 1
    ri <- mi; while (ri < length(fl) && fl[ri] > half) ri <- ri + 1
    fwhm <- freqs[ri] - freqs[li]
    gstd <- (fwhm / 2) / sqrt(2 * log(2))
    if (!is.finite(gstd) || gstd <= 0) gstd <- stats::median(diff(freqs)) * 2
    gstd <- min(max(gstd, bw_min), bw_max)
    guesses[[length(guesses) + 1L]] <- c(cf, mh, gstd)
    fl <- fl - .specparam_gaussian(freqs, cf, mh, gstd)
  }
  if (!length(guesses)) return(matrix(numeric(0), 0, 3))
  g <- do.call(rbind, guesses)
  obj <- function(p) {
    m <- matrix(p, ncol = 3, byrow = TRUE)
    pred <- Reduce(`+`, lapply(seq_len(nrow(m)),
                   function(k) .specparam_gaussian(freqs, m[k, 1], m[k, 2], m[k, 3])))
    sum((flat - pred)^2)
  }
  lower <- as.vector(t(cbind(g[, 1] - 2, 0, bw_min)))
  upper <- as.vector(t(cbind(g[, 1] + 2, Inf, bw_max)))
  opt <- tryCatch(
    stats::optim(as.vector(t(g)), obj, method = "L-BFGS-B",
                 lower = lower, upper = upper),
    error = function(e) list(par = as.vector(t(g))))
  matrix(opt$par, ncol = 3, byrow = TRUE)
}

## Fit one channel's PSD. Returns aperiodic params, peak table, fit stats, and
## the log spectra needed for plotting.
.specparam_fit_one <- function(freqs, power, freq_range, peak_width_limits,
                               max_n_peaks, aperiodic_mode, min_peak_height,
                               peak_threshold) {
  sel <- freqs >= freq_range[1] & freqs <= freq_range[2] & freqs > 0 & power > 0
  f <- freqs[sel]; lp <- log10(power[sel])
  ap <- .specparam_robust_ap(f, lp, aperiodic_mode)
  flat <- lp - .specparam_gen_ap(f, ap, aperiodic_mode)
  g <- .specparam_fit_peaks(f, flat, max_n_peaks, min_peak_height,
                            peak_threshold, peak_width_limits)
  peak_fit <- if (nrow(g)) {
    Reduce(`+`, lapply(seq_len(nrow(g)),
           function(k) .specparam_gaussian(f, g[k, 1], g[k, 2], g[k, 3])))
  } else numeric(length(f))
  ap <- .specparam_simple_ap(f, lp - peak_fit, aperiodic_mode)
  ap_fit <- .specparam_gen_ap(f, ap, aperiodic_mode)
  model <- ap_fit + peak_fit
  ss_tot <- sum((lp - mean(lp))^2)
  r2 <- if (ss_tot > 0) 1 - sum((lp - model)^2) / ss_tot else NA_real_
  peaks <- if (nrow(g)) {
    ord <- order(g[, 1])
    data.frame(CF = g[ord, 1], PW = g[ord, 2], BW = 2 * g[ord, 3])
  } else data.frame(CF = numeric(0), PW = numeric(0), BW = numeric(0))
  list(aperiodic = ap, peaks = peaks, r_squared = r2,
       error = sqrt(mean((lp - model)^2)),
       freqs = f, log_psd = lp, model = model, ap_fit = ap_fit)
}

#' Spectral parameterization (specparam / FOOOF)
#'
#' Parameterizes each channel's power spectrum into an aperiodic 1/f component
#' (offset, exponent, and an optional knee) plus a set of Gaussian periodic
#' peaks, following the specparam / FOOOF algorithm (Donoghue et al. 2020). The
#' aperiodic exponent is a marker of excitation/inhibition balance and cortical
#' state; the peaks capture genuine oscillations separated from the 1/f
#' background. The fit runs on the Welch PSD computed by the same machinery as
#' [bandPower()].
#'
#' @param x A PhysioExperiment object. 3D (epoched) assays are averaged over
#'   trials before fitting.
#' @param freq_range Numeric length-2 frequency range in Hz to fit
#'   (default: \code{c(1, 45)}).
#' @param peak_width_limits Numeric length-2 minimum and maximum peak bandwidth
#'   in Hz (default: \code{c(1, 12)}).
#' @param max_n_peaks Maximum number of periodic peaks to fit (default: 6).
#' @param aperiodic_mode \code{"fixed"} for a pure 1/f (offset, exponent) or
#'   \code{"knee"} to also fit a knee frequency.
#' @param min_peak_height Minimum peak height above the aperiodic fit, in
#'   \code{log10} power (default: 0.05).
#' @param peak_threshold Minimum peak height in units of the standard deviation
#'   of the flattened spectrum (default: 2).
#' @param nperseg Welch segment length in samples (default: about a 2-second
#'   segment for adequate low-frequency resolution).
#' @param assay_name Input assay name (default: the default assay).
#' @return An object of class \code{"specparam_result"}: a list with
#'   \code{aperiodic} (a data.frame of per-channel offset, exponent, and knee),
#'   \code{peaks} (a data.frame of per-channel center frequency CF, power PW, and
#'   bandwidth BW), \code{fit} (per-channel R-squared and error), \code{spectra}
#'   (per-channel log spectra and fitted model for plotting), and the settings
#'   used.
#' @references
#' Donoghue, T., et al. (2020). Parameterizing neural power spectra into
#' periodic and aperiodic components. Nature Neuroscience, 23(12), 1655-1665.
#' @seealso [plotSpecparam()], [specparamBiomarker()], [bandPower()]
#' @export
#' @examples
#' \dontrun{
#' pe <- make_pe_2d(n_time = 4000, n_channels = 4, sr = 250)
#' sp <- specparam(pe, freq_range = c(1, 40))
#' sp$aperiodic
#' }
specparam <- function(x, freq_range = c(1, 45), peak_width_limits = c(1, 12),
                      max_n_peaks = 6L, aperiodic_mode = c("fixed", "knee"),
                      min_peak_height = 0.05, peak_threshold = 2,
                      nperseg = NULL, assay_name = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))
  aperiodic_mode <- match.arg(aperiodic_mode)
  stopifnot(is.numeric(freq_range), length(freq_range) == 2,
            freq_range[1] < freq_range[2])
  stopifnot(is.numeric(peak_width_limits), length(peak_width_limits) == 2)

  sr <- samplingRate(x)
  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  if (length(dim(data)) == 3) data <- apply(data, c(1, 2), mean)
  data <- as.matrix(data)
  n_ch <- ncol(data)
  if (is.null(nperseg)) nperseg <- min(nrow(data), max(256L, 2L * as.integer(sr)))

  cd <- SummarizedExperiment::colData(x)
  ch_labels <- if ("label" %in% colnames(cd)) as.character(cd$label) else
    paste0("Ch", seq_len(n_ch))

  fits <- lapply(seq_len(n_ch), function(ch) {
    ps <- .welchPSD(data[, ch], sr, nperseg = nperseg)
    .specparam_fit_one(ps$frequencies, ps$power, freq_range, peak_width_limits,
                       max_n_peaks, aperiodic_mode, min_peak_height, peak_threshold)
  })

  aperiodic <- data.frame(
    channel = ch_labels,
    offset = vapply(fits, function(z) z$aperiodic[["offset"]], numeric(1)),
    exponent = vapply(fits, function(z) z$aperiodic[["exponent"]], numeric(1)),
    stringsAsFactors = FALSE)
  if (aperiodic_mode == "knee") {
    aperiodic$knee <- vapply(fits, function(z) z$aperiodic[["knee"]], numeric(1))
  }

  peaks <- do.call(rbind, lapply(seq_len(n_ch), function(ch) {
    pk <- fits[[ch]]$peaks
    if (!nrow(pk)) return(NULL)
    cbind(channel = ch_labels[ch], pk, stringsAsFactors = FALSE)
  }))
  if (is.null(peaks)) {
    peaks <- data.frame(channel = character(0), CF = numeric(0),
                        PW = numeric(0), BW = numeric(0), stringsAsFactors = FALSE)
  }

  fit_df <- data.frame(
    channel = ch_labels,
    r_squared = vapply(fits, function(z) z$r_squared, numeric(1)),
    error = vapply(fits, function(z) z$error, numeric(1)),
    stringsAsFactors = FALSE)

  spectra <- stats::setNames(lapply(fits, function(z)
    list(freqs = z$freqs, log_psd = z$log_psd, model = z$model,
         ap_fit = z$ap_fit)), ch_labels)

  structure(list(aperiodic = aperiodic, peaks = peaks, fit = fit_df,
                 spectra = spectra, freq_range = freq_range,
                 aperiodic_mode = aperiodic_mode, channels = ch_labels),
            class = "specparam_result")
}

#' @param x A \code{specparam_result}.
#' @param ... Ignored.
#' @rdname specparam
#' @exportS3Method print specparam_result
print.specparam_result <- function(x, ...) {
  cat(sprintf("<specparam_result> %d channel(s), mode = %s, %g-%g Hz\n",
              length(x$channels), x$aperiodic_mode,
              x$freq_range[1], x$freq_range[2]))
  cat(sprintf("  aperiodic exponent: %.3f - %.3f (mean %.3f)\n",
              min(x$aperiodic$exponent), max(x$aperiodic$exponent),
              mean(x$aperiodic$exponent)))
  cat(sprintf("  peaks: %d total; model R-squared mean %.3f\n",
              nrow(x$peaks), mean(x$fit$r_squared)))
  invisible(x)
}

#' Plot a specparam model fit
#'
#' Plots the observed log power spectrum for one channel together with the full
#' specparam model and the aperiodic-only fit.
#'
#' @param result A \code{specparam_result} from [specparam()].
#' @param channel Channel index or label to plot (default: 1).
#' @return A \code{ggplot} object.
#' @seealso [specparam()]
#' @export
#' @examples
#' \dontrun{
#' sp <- specparam(make_pe_2d(n_time = 4000, sr = 250))
#' plotSpecparam(sp, channel = 1)
#' }
plotSpecparam <- function(result, channel = 1) {
  stopifnot(inherits(result, "specparam_result"))
  ch <- if (is.numeric(channel)) result$channels[channel] else channel
  sp <- result$spectra[[ch]]
  if (is.null(sp)) stop("channel not found: ", channel, call. = FALSE)
  df <- data.frame(
    frequency = rep(sp$freqs, 3),
    power = c(sp$log_psd, sp$model, sp$ap_fit),
    component = rep(c("observed", "full model", "aperiodic"),
                    each = length(sp$freqs)))
  df$component <- factor(df$component,
                         levels = c("observed", "full model", "aperiodic"))
  ggplot2::ggplot(df, ggplot2::aes(x = frequency, y = power,
                                   colour = component, linetype = component)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::scale_linetype_manual(values = c("solid", "solid", "dashed")) +
    ggplot2::labs(x = "Frequency (Hz)", y = "log10 power",
                  title = sprintf("specparam fit: %s", ch),
                  colour = NULL, linetype = NULL) +
    ggplot2::theme_minimal()
}

#' Aperiodic parameter as a PhysioBiomarker
#'
#' Wraps a specparam aperiodic parameter (by default the exponent) for one
#' channel as a reliability-characterised \code{\link[PhysioCore]{PhysioBiomarker}}.
#'
#' @param result A \code{specparam_result} from [specparam()].
#' @param channel Channel index or label (default: 1).
#' @param param Which aperiodic parameter to export: \code{"exponent"},
#'   \code{"offset"}, or \code{"knee"}.
#' @param reliability Named list of reliability indices to attach (default:
#'   \code{icc} and \code{sem} placeholders).
#' @return A \code{\link[PhysioCore]{PhysioBiomarker}}.
#' @seealso [specparam()], [PhysioCore::physioBiomarker()]
#' @export
#' @examples
#' \dontrun{
#' sp <- specparam(make_pe_2d(n_time = 4000, sr = 250))
#' specparamBiomarker(sp, channel = 1)
#' }
specparamBiomarker <- function(result, channel = 1,
                               param = c("exponent", "offset", "knee"),
                               reliability = list(icc = NA_real_, sem = NA_real_)) {
  stopifnot(inherits(result, "specparam_result"))
  param <- match.arg(param)
  if (param == "knee" && !"knee" %in% names(result$aperiodic)) {
    stop("this result was fit in 'fixed' mode and has no knee.", call. = FALSE)
  }
  ch <- if (is.numeric(channel)) result$channels[channel] else channel
  row <- which(result$aperiodic$channel == ch)
  if (!length(row)) stop("channel not found: ", channel, call. = FALSE)
  unit <- switch(param, exponent = "a.u.", offset = "log10(power)",
                 knee = "Hz")
  ver <- tryCatch(as.character(utils::packageVersion("PhysioAnalysis")),
                  error = function(e) NA_character_)
  physioBiomarker(
    value = result$aperiodic[[param]][row],
    name = paste0("aperiodic_", param), unit = unit, reliability = reliability,
    provenance = list(assay = "specparam",
                      band = sprintf("%g-%gHz", result$freq_range[1],
                                     result$freq_range[2]),
                      method = paste0("specparam:", result$aperiodic_mode),
                      software_version = ver),
    interpretation = sprintf("channel %s", ch))
}
