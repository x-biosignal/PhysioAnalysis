library(testthat)
library(PhysioAnalysis)

# VAL-01: numerical parity of PhysioAnalysis's SPM/RFT against an ACTUAL spm1d +
# rft1d run (spm1d 0.4.53, rft1d 0.2.5). The fixture bundles the input data and
# spm1d's outputs so spmTTest/spmPairedTTest/spmAnova run on identical inputs.
# Provenance and generation: data-raw/spm1d_reference.{py,R}.
#
# Parity targets and what is checked here (all tight):
#   * RFT critical t/F thresholds  vs rft1d.{t,f}.isf   (rel < 1e-2; observed < 7e-3)
#   * RFT cluster-extent p (.rftClusterPValue) vs rft1d.p_cluster (~1e-15; VAL-07)
#   * t/F statistic fields          vs spm1d            (< 1e-6; observed ~1e-13)
#   * residual FWHM                 vs spm1d            (< 1e-4; observed ~0)
#   * field-level critical height   vs per-case rft1d   (rel < 1e-2; observed < 6e-3)
#   * cluster count and endpoints   vs spm1d            (exact)
#
# KNOWN divergences (NOT asserted here):
#   * spm1d's anova1 Fstar uses df2-1 and its own effective smoothness, so it
#     differs from rft1d.f.isf at the reported FWHM; PhysioAnalysis matches
#     rft1d (the pure RFT theory), which is the correct target.
#   * Full-pipeline cluster P: VAL-07 fixed the .rftClusterPValue formula (it now
#     matches rft1d.p_cluster exactly -- see the grid test below). spm1d's
#     *reported* cluster P additionally evaluates the RFT cluster probability at
#     an idiosyncratic internal cluster height (interpolation-aware, doubled for
#     two-tailed) rather than the cluster-defining threshold, so the end-to-end
#     cluster P still differs; PhysioAnalysis uses the standard threshold-based
#     Worsley/Friston convention.

fx_path <- test_path("fixtures", "spm1d-reference.rds")

test_that("RFT critical t and F thresholds match rft1d across a grid", {
  skip_if(!file.exists(fx_path), "spm1d/rft1d reference fixture not bundled")
  fx <- readRDS(fx_path)
  t_now <- mapply(function(a, df, rc) PhysioAnalysis:::.rftCriticalT(a, df, rc, TRUE),
                  fx$rft_t$alpha, fx$rft_t$df, fx$rft_t$resel)
  expect_lt(max(abs(t_now - fx$rft_t$crit) / fx$rft_t$crit), 1e-2)
  f_now <- mapply(function(a, d1, d2, rc) PhysioAnalysis:::.rftCriticalF(a, d1, d2, rc),
                  fx$rft_f$alpha, fx$rft_f$df1, fx$rft_f$df2, fx$rft_f$resel)
  expect_lt(max(abs(f_now - fx$rft_f$crit) / fx$rft_f$crit), 1e-2)
})

test_that("RFT cluster-extent probabilities match rft1d.p_cluster (VAL-07)", {
  skip_if(!file.exists(fx_path), "spm1d/rft1d reference fixture not bundled")
  fx <- readRDS(fx_path)
  # .rftClusterPValue(k*fwhm, R, fwhm, rho0(u), rho1(u), tails=1) must equal
  # rft1d.{t,f}.p_cluster(k, u, df, Q, fwhm) -- the corrected cluster-extent
  # formula (E[k] = rho0/rho1, not the earlier e_m/e_n form).
  gt <- fx$clusterp_t
  p_t <- mapply(function(u, df, R, k, fw) {
    e <- PhysioAnalysis:::.ecDensityT(u, df)
    PhysioAnalysis:::.rftClusterPValue(k * fw, R, fw, e["rho0"], e["rho1"], 1)
  }, gt$u, gt$df, gt$R, gt$k, gt$fwhm)
  expect_equal(as.numeric(p_t), gt$p, tolerance = 1e-9)
  gf <- fx$clusterp_f
  p_f <- mapply(function(u, d1, d2, R, k, fw) {
    e <- PhysioAnalysis:::.ecDensityF(u, d1, d2)
    PhysioAnalysis:::.rftClusterPValue(k * fw, R, fw, e["rho0"], e["rho1"], 1)
  }, gf$u, gf$df1, gf$df2, gf$R, gf$k, gf$fwhm)
  expect_equal(as.numeric(p_f), gf$p, tolerance = 1e-9)
})

# compare a PhysioAnalysis spm_result to an spm1d case (field/fwhm/threshold/clusters)
# check_clusters = FALSE for anova, where spm1d's Fstar convention gives a
# different threshold (and hence different supra-threshold clusters) than rft1d.
expect_spm_parity <- function(res, case, field, check_clusters = TRUE) {
  # statistic field (t or F): essentially exact
  expect_equal(as.numeric(field), case$field, tolerance = 1e-6)
  # residual FWHM: same Kiebel/Friston estimator
  expect_lt(abs(res$fwhm - case$fwhm) / case$fwhm, 1e-4)
  # field-level critical height vs rft1d at this case's (df, Q, FWHM)
  expect_lt(abs(res$threshold - case$rft_crit) / case$rft_crit, 1e-2)
  if (!check_clusters) return(invisible())
  # cluster count and endpoints (spm1d endpoints are 0-based; R clusters 1-based)
  expect_equal(length(res$clusters), length(case$clusters))
  if (length(res$clusters) && length(case$clusters)) {
    r_starts <- vapply(res$clusters, `[[`, numeric(1), "start")
    ord <- order(r_starts)
    for (i in seq_along(case$clusters)) {
      rc <- res$clusters[[ord[i]]]
      ep <- case$clusters[[i]]$endpoints
      expect_equal(rc$start, ep[1] + 1)
      expect_equal(rc$end,   ep[2] + 1)
    }
  }
}

test_that("two-sample SPM{t} matches spm1d (field, FWHM, threshold, clusters)", {
  skip_if(!file.exists(fx_path), "spm1d/rft1d reference fixture not bundled")
  fx <- readRDS(fx_path)
  for (nm in c("effect", "null")) {
    c <- fx$ttest2[[nm]]
    res <- spmTTest(c$data, group1 = c$group1, group2 = c$group2, alpha = c$alpha)
    expect_spm_parity(res, c, res$t)
  }
})

test_that("paired SPM{t} matches spm1d (field, FWHM, threshold, clusters)", {
  skip_if(!file.exists(fx_path), "spm1d/rft1d reference fixture not bundled")
  fx <- readRDS(fx_path)
  c <- fx$paired$effect
  res <- spmPairedTTest(c$data, condition1 = c$condition1, condition2 = c$condition2,
                        alpha = c$alpha)
  expect_spm_parity(res, c, res$t)
})

test_that("one-way SPM{F} ANOVA matches spm1d field/FWHM and rft1d threshold", {
  skip_if(!file.exists(fx_path), "spm1d/rft1d reference fixture not bundled")
  fx <- readRDS(fx_path)
  c <- fx$anova1$effect
  res <- spmAnova(c$data, c$groups, alpha = c$alpha)
  expect_equal(res$df1, c$df1)
  expect_equal(res$df2, c$df2)          # F-statistic residual df (spm1d reports df2-1)
  # clusters not compared: spm1d's anova Fstar (recorded in c$zstar) diverges
  # from rft1d, so its supra-threshold clusters differ (see VAL-07).
  expect_spm_parity(res, c, res$f, check_clusters = FALSE)
})
