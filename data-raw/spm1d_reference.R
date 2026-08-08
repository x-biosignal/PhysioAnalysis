#!/usr/bin/env Rscript
# Build tests/testthat/fixtures/spm1d-reference.rds for the spm1d/rft1d parity
# test (VAL-01), from the extraction in data-raw/spm1d_reference.py.
#
# Provenance: fields (t/F), FWHM, critical heights and cluster inference are
# captured from an ACTUAL spm1d 0.4.53 + rft1d 0.2.5 run (see the .py for seeds,
# smoothing, and effect sizes). spm1d stores fields as (Nobs x Q); this builder
# transposes to PhysioAnalysis's (Q x Nobs). `rft_crit` per case is the rft1d
# critical threshold at that case's exact (df, Q, FWHM) -- the pure-RFT parity
# target PhysioAnalysis is checked against (spm1d's own anova Fstar uses a
# different effective-smoothness/df convention; recorded but not the target).
#
# Run:  micromamba run -n spm1d python data-raw/spm1d_reference.py   # -> JSON
#       OUT=/tmp/spm1d_ref Rscript data-raw/spm1d_reference.R        # -> .rds

WD <- Sys.getenv("OUT", "/tmp/spm1d_ref")
ref <- jsonlite::fromJSON(file.path(WD, "spm1d_reference.json"), simplifyMatrix = TRUE)

mat <- function(x) { m <- as.matrix(x); storage.mode(m) <- "double"; m }
cl_list <- function(cl) {
  # normalise spm1d clusters (data.frame or list) to a list of {endpoints, P, extent}
  if (is.null(cl) || (is.data.frame(cl) && nrow(cl) == 0)) return(list())
  if (is.data.frame(cl))
    return(lapply(seq_len(nrow(cl)), function(i)
      list(endpoints = as.numeric(cl$endpoints[[i]]), P = cl$P[i], extent = cl$extent[i])))
  lapply(cl, function(c) list(endpoints = as.numeric(c$endpoints), P = c$P, extent = c$extent))
}

t2 <- function(cs) {
  yA <- mat(cs$yA); yB <- mat(cs$yB)
  list(data = t(rbind(yA, yB)), group1 = seq_len(nrow(yA)),
       group2 = nrow(yA) + seq_len(nrow(yB)),
       alpha = cs$alpha, df = cs$df, fwhm = cs$fwhm, zstar = cs$zstar,
       rft_crit = cs$rft_crit, field = as.numeric(cs$tfield), clusters = cl_list(cs$clusters))
}
pr <- function(cs) {
  yA <- mat(cs$yA); yB <- mat(cs$yB)
  list(data = t(rbind(yA, yB)), condition1 = seq_len(nrow(yA)),
       condition2 = nrow(yA) + seq_len(nrow(yB)),
       alpha = cs$alpha, df = cs$df, fwhm = cs$fwhm, zstar = cs$zstar,
       rft_crit = cs$rft_crit, field = as.numeric(cs$tfield), clusters = cl_list(cs$clusters))
}
an <- function(cs) {
  list(data = t(mat(cs$Y)), groups = factor(cs$labels),
       alpha = cs$alpha, df1 = cs$df1, df2 = cs$df2, fwhm = cs$fwhm,
       zstar = cs$zstar, rft_crit = cs$rft_crit, spm1d_df = cs$spm1d_df,
       field = as.numeric(cs$Ffield), clusters = cl_list(cs$clusters))
}

fixture <- list(
  rft_t = ref$rft_t, rft_f = ref$rft_f,
  clusterp_t = ref$clusterp_t, clusterp_f = ref$clusterp_f,
  ttest2 = list(effect = t2(ref$ttest2$effect), null = t2(ref$ttest2$null)),
  paired = list(effect = pr(ref$paired$effect)),
  anova1 = list(effect = an(ref$anova1$effect)),
  provenance = c(ref$provenance, list(
    parity_targets = paste(
      "thresholds vs rft1d (pure RFT); t/F fields and FWHM vs spm1d.",
      "KNOWN divergences (follow-up VAL-07): spm1d anova1 Fstar uses df2-1 and its",
      "own effective smoothness (recorded as $zstar, not the target); cluster-extent",
      "p-values use PhysioAnalysis's Friston-1994 GRF approx which differs from spm1d.")))
)

out <- file.path("tests", "testthat", "fixtures", "spm1d-reference.rds")
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
saveRDS(fixture, out, version = 2, compress = "xz")
message(sprintf("wrote %s (%d bytes); spm1d %s rft1d %s",
                out, file.info(out)$size, ref$provenance$spm1d, ref$provenance$rft1d))
