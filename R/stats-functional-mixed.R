# Functional (waveform) mixed-effects models for MULTILEVEL curves.
#
# Movement data is hierarchical: many strides per subject, many subjects per
# group. SPM (and ordinary curve statistics) treat every curve as independent,
# which pseudo-replicates -- a subject with 30 strides counts 30x, inflating the
# apparent evidence. This module fits the functional random-intercept model:
#   * decompose the waveform variance, point by point, into a BETWEEN-subject and
#     a WITHIN-subject (stride-to-stride) component, giving a functional ICC
#     curve (how much of the variance is stable subject identity), and
#   * test a between-subject fixed effect (e.g. patients vs controls) at the
#     SUBJECT level -- on subject-mean curves -- so the inference respects the
#     multilevel structure instead of counting strides as independent.
# Dependency-free base R (unbalanced one-way random-effects components +
# subject-level permutation inference).

#' Functional random-intercept (multilevel) model for waveforms
#'
#' Fits the subject random-intercept model to a set of curves nested in subjects.
#' It splits the pointwise variance into between- and within-subject components
#' (a functional ICC curve) and, when `group` is given, tests the between-subject
#' fixed effect on the subject-mean curves so strides are not treated as
#' independent.
#'
#' @param curves An `N x P` matrix of curves: `N` observations (e.g. strides
#'   pooled over subjects), `P` domain points.
#' @param subject Length-`N` grouping identifying the random-intercept level
#'   (the subject each curve belongs to).
#' @param group Optional between-subject fixed effect (length `N`), constant
#'   within each subject (e.g. a patient/control label).
#' @param n_perm Permutations for the curve-wide fixed-effect p-value (default
#'   1000; 0 to skip).
#' @param seed Optional RNG seed for the permutation test.
#' @return a `functional_mixed` object: `mean_curve`, `var_between`,
#'   `var_within`, `icc` (functional ICC curve), `subject_means` (one row per
#'   subject), `n_subjects`, and -- when `group` is given -- `f_curve`,
#'   `p_pointwise`, `p_global` (subject-level fixed-effect test) plus
#'   `group_diff` (the signed difference curve for two groups).
#' @references Morris JS (2015) Functional regression, Annu Rev Stat Appl; Pini &
#'   Vantini (2017) interval-wise testing.
#' @seealso [functionalRegression()], [spmTTest()]
#' @export
#' @examples
#' set.seed(1)
#' P <- 101; t <- seq(0, 1, length.out = P)
#' make <- function(off) t(sapply(1:12, function(k) off + sin(2 * pi * t) + rnorm(P, 0, 0.1)))
#' ctrl <- do.call(rbind, lapply(rnorm(6, 0, 0.5), make))
#' pat  <- do.call(rbind, lapply(rnorm(6, 1.0, 0.5), make))   # shifted group
#' Y <- rbind(ctrl, pat)
#' subj <- rep(1:12, each = 12); grp <- rep(c("ctrl", "pat"), each = 72)
#' fm <- functionalMixedModel(Y, subj, grp, n_perm = 200)
#' fm$p_global                         # group difference, subject-level correct
functionalMixedModel <- function(curves, subject, group = NULL,
                                  n_perm = 1000L, seed = NULL) {
  Y <- as.matrix(curves); N <- nrow(Y); P <- ncol(Y)
  subject <- as.factor(subject)
  if (length(subject) != N) stop("`subject` must have length N.", call. = FALSE)
  subs <- levels(subject); ns <- length(subs)
  if (ns < 2L) stop("need >= 2 subjects.", call. = FALSE)
  ni <- as.integer(table(subject))
  grand <- colMeans(Y)
  # subject-mean curves + within/between sums of squares (pointwise)
  Sbar <- matrix(0, ns, P, dimnames = list(subs, NULL))
  ss_within <- numeric(P)
  for (i in seq_len(ns)) {
    Yi <- Y[subject == subs[i], , drop = FALSE]
    mi <- colMeans(Yi); Sbar[i, ] <- mi
    ss_within <- ss_within + colSums(sweep(Yi, 2L, mi)^2)
  }
  ss_between <- colSums(ni * sweep(Sbar, 2L, grand)^2)
  ms_within <- ss_within / (N - ns)
  ms_between <- ss_between / (ns - 1)
  n0 <- (N - sum(ni^2) / N) / (ns - 1)               # unbalanced design constant
  var_within <- ms_within
  var_between <- pmax((ms_between - ms_within) / n0, 0)
  icc <- var_between / pmax(var_between + var_within, .Machine$double.eps)

  out <- list(mean_curve = grand, var_between = var_between, var_within = var_within,
              icc = icc, subject_means = Sbar, n_subjects = ns, n_obs = N, P = P,
              subject_n = ni)
  # between-subject fixed effect, tested at the SUBJECT level
  if (!is.null(group)) {
    group <- as.factor(group)
    if (length(group) != N) stop("`group` must have length N.", call. = FALSE)
    # one group label per subject (must be constant within subject)
    sg <- vapply(subs, function(s) {
      g <- unique(as.character(group[subject == s]))
      if (length(g) != 1L) stop("`group` must be constant within each subject.", call. = FALSE)
      g
    }, character(1))
    sg <- as.factor(sg)
    fcur <- .fmm_group_F(Sbar, sg)
    p_point <- stats::pf(fcur$F, fcur$df1, fcur$df2, lower.tail = FALSE)
    p_global <- NA_real_
    if (n_perm > 0) {
      if (!is.null(seed)) set.seed(seed)
      obs_max <- max(fcur$F)
      ge <- 0L
      for (b in seq_len(n_perm)) {
        fp <- .fmm_group_F(Sbar, sg[sample.int(ns)])$F
        ge <- ge + (max(fp) >= obs_max)
      }
      p_global <- (ge + 1) / (n_perm + 1)
    }
    out$f_curve <- fcur$F; out$df1 <- fcur$df1; out$df2 <- fcur$df2
    out$p_pointwise <- p_point; out$p_global <- p_global
    out$groups <- levels(sg)
    if (nlevels(sg) == 2L)
      out$group_diff <- colMeans(Sbar[sg == levels(sg)[1], , drop = FALSE]) -
                        colMeans(Sbar[sg == levels(sg)[2], , drop = FALSE])
  }
  structure(out, class = "functional_mixed")
}

# Pointwise one-way ANOVA F across groups on subject-mean curves (ns x P).
.fmm_group_F <- function(Sbar, sg) {
  ns <- nrow(Sbar); G <- nlevels(sg); glev <- levels(sg)
  grand <- colMeans(Sbar)
  ss_bg <- numeric(ncol(Sbar)); ss_wg <- numeric(ncol(Sbar))
  for (g in glev) {
    Sg <- Sbar[sg == g, , drop = FALSE]; mg <- colMeans(Sg)
    ss_bg <- ss_bg + nrow(Sg) * (mg - grand)^2
    ss_wg <- ss_wg + colSums(sweep(Sg, 2L, mg)^2)
  }
  df1 <- G - 1; df2 <- ns - G
  Fv <- (ss_bg / df1) / pmax(ss_wg / df2, .Machine$double.eps)
  list(F = Fv, df1 = df1, df2 = df2)
}

#' @export
print.functional_mixed <- function(x, ...) {
  cat(sprintf("Functional mixed model -- %d curves in %d subjects, %d points\n",
              x$n_obs, x$n_subjects, x$P))
  cat(sprintf("  functional ICC: mean %.3f (range %.3f-%.3f)\n",
              mean(x$icc), min(x$icc), max(x$icc)))
  if (!is.null(x$f_curve)) {
    cat(sprintf("  group effect (subject-level): groups {%s}, peak F = %.2f",
                paste(x$groups, collapse = ", "), max(x$f_curve)))
    if (!is.na(x$p_global)) cat(sprintf(", curve-wide p = %.3f", x$p_global))
    cat("\n")
  }
  invisible(x)
}
