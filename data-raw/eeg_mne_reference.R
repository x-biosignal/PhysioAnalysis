#!/usr/bin/env Rscript
# Build tests/testthat/fixtures/eeg-mne-reference.rds for the EEG<->MNE parity
# test (VAL-04), from the MNE extraction in data-raw/eeg_mne_reference.py.
#
# Provenance: band powers are MNE (mne.time_frequency.psd_array_welch); PLV/wPLI
# are mne_connectivity.spectral_connectivity_epochs; all on the same seeded
# signal (stored) so the R side runs on identical input. Band power matches MNE
# to < 0.5%; PLV/wPLI agree on coupling structure (~0.95) but use a different
# (epoch-spectral vs continuous-Hilbert) estimator -- documented in the test.
#
# Run:  micromamba run -n mne python data-raw/eeg_mne_reference.py   # -> JSON
#       OUT=/tmp/val04 Rscript data-raw/eeg_mne_reference.R          # -> .rds

WD <- Sys.getenv("OUT", "/tmp/val04")
ref <- jsonlite::fromJSON(file.path(WD, "eeg_mne_reference.json"), simplifyMatrix = TRUE)

data <- t(as.matrix(ref$data))          # (nch x n) -> (n samples x nch)
storage.mode(data) <- "double"
bp <- as.matrix(as.data.frame(ref$band_powers))   # nch x nbands
colnames(bp) <- ref$bands

fixture <- list(
  sf = ref$sf, nch = ref$nch,
  data = data,                          # samples x channels
  band_powers = bp,                     # nch x nbands (MNE Welch)
  bands = ref$bands,
  plv = as.matrix(ref$plv),             # nch x nch (mne_connectivity, alpha)
  wpli = as.matrix(ref$wpli),           # nch x nch (mne_connectivity, theta)
  provenance = list(
    description = "MNE Welch band power + mne_connectivity PLV/wPLI on a fixed seeded EEG signal.",
    mne = ref$mne,
    generator = "data-raw/eeg_mne_reference.{py,R}",
    note = paste("PSD band power matches MNE tightly; plv()/wPLI() are continuous",
                 "Hilbert estimators vs mne_connectivity's epoch-spectral, so they",
                 "agree on coupling structure (~0.95) but are not identical.")
  )
)

out <- file.path("tests", "testthat", "fixtures", "eeg-mne-reference.rds")
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
saveRDS(fixture, out, version = 2, compress = "xz")
message(sprintf("wrote %s (%d bytes; %d ch x %d samples, mne %s)",
                out, file.info(out)$size, ncol(data), nrow(data), ref$mne))
