library(testthat)
library(PhysioAnalysis)

# VAL-04: numeric EEG<->MNE parity, replacing the manuscript's documentation-level
# "MNE-equivalent" correspondence with a bundled, runnable reference captured from
# MNE-Python 1.12.1 + mne_connectivity 0.8.1 on a fixed seeded EEG signal (stored
# in the fixture so the R side runs on identical input). Provenance and
# generation: data-raw/eeg_mne_reference.{py,R}.
#
# What is checked:
#   * Welch band power (theta/alpha/beta/gamma) vs MNE psd_array_welch -- tight
#     (cor > 0.999, < 5% per channel). This reproduces the manuscript's headline
#     "alpha power" cross-validation with MNE.
#   * PLV (alpha) / wPLI (theta) vs mne_connectivity -- structural agreement
#     (cor > 0.9); see the note below.
#
# Documented differences:
#   * The delta band (0.5-4 Hz) is NOT asserted: MNE's Welch detrends each segment
#     (removing DC/drift) while PhysioAnalysis's .welchPSD does not, so the DC-
#     dominated delta band diverges. Bands >= theta are unaffected.
#   * plv()/wPLI() are continuous Hilbert-phase estimators; mne_connectivity uses
#     an epoch-spectral estimator. They agree on coupling STRUCTURE (~0.95
#     correlation here) but are not identical, so the manuscript's PLV=0.9996 /
#     wPLI=0.9994 figures (a different estimator pairing) are not reproduced.

fx_path <- test_path("fixtures", "eeg-mne-reference.rds")

.load_pe <- function(fx) {
  PhysioExperiment(assays = list(raw = fx$data), samplingRate = fx$sf)
}

test_that("Welch band power matches MNE psd_array_welch (theta/alpha/beta/gamma)", {
  skip_if(!file.exists(fx_path), "EEG/MNE reference fixture not bundled")
  fx <- readRDS(fx_path)
  bp <- as.data.frame(bandPower(.load_pe(fx)))
  for (b in c("theta", "alpha", "beta", "gamma")) {
    r_col <- as.numeric(bp[[b]])
    m_col <- as.numeric(fx$band_powers[, b])
    expect_gt(cor(r_col, m_col), 0.999, label = sprintf("%s power cor vs MNE", b))
    expect_lt(max(abs(r_col - m_col) / m_col), 0.05,
              label = sprintf("%s power max rel diff vs MNE", b))
  }
})

test_that("PLV and wPLI agree structurally with mne_connectivity", {
  skip_if(!file.exists(fx_path), "EEG/MNE reference fixture not bundled")
  fx <- readRDS(fx_path)
  pe <- .load_pe(fx)
  ut <- upper.tri(matrix(0, fx$nch, fx$nch))
  plv_m <- matrix(as.numeric(fx$plv), fx$nch)
  wpli_m <- matrix(as.numeric(fx$wpli), fx$nch)
  # continuous-Hilbert (PhysioAnalysis) vs epoch-spectral (mne_connectivity):
  # correlated on coupling structure, not identical.
  expect_gt(cor(plv(pe, c(8, 12))[ut], plv_m[ut]), 0.90)
  expect_gt(cor(abs(wPLI(pe, c(4, 8))[ut]), abs(wpli_m[ut])), 0.90)
})
