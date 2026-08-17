# Golden-reference generator for PhysioAnalysis (WSF-08).
#
# Each golden is an INDEPENDENT reference value (igraph / base R / an analytic
# one-sided FFT PSD formula) -- NOT captured from the package's own output.
# The companion test file (tests/testthat/test-golden.R) rebuilds the SAME
# deterministic input, runs the PACKAGE function, and compares to these goldens.
#
# Run from the repo root, e.g.:
#   _R_CHECK_FORCE_SUGGESTS_=false Rscript physio-ecosystem/PhysioAnalysis/data-raw/golden.R
#
# Determinism: every input is built under an explicit set.seed(), so rerunning
# regenerates byte-identical .rds fixtures.

suppressMessages({
  if (requireNamespace("devtools", quietly = TRUE) &&
      !"PhysioAnalysis" %in% loadedNamespaces()) {
    # Load PhysioCore (dependency) then PhysioAnalysis, if not already loaded.
    core <- "physio-ecosystem/PhysioCore"
    self <- "physio-ecosystem/PhysioAnalysis"
    if (!dir.exists(core)) {                     # allow running from pkg dir
      core <- file.path("..", "PhysioCore")
      self <- "."
    }
    if (!"PhysioCore" %in% loadedNamespaces()) {
      try(devtools::load_all(core, quiet = TRUE), silent = TRUE)
    }
    devtools::load_all(self, quiet = TRUE)
  }
})

stopifnot(requireNamespace("igraph", quietly = TRUE))

# Resolve the golden dir independent of testthat::test_path (so the script runs
# outside a test context). Mirror .golden_dir(): <pkg>/tests/testthat/_golden.
.pkg_root <- local({
  cand <- c("physio-ecosystem/PhysioAnalysis", ".")
  hit <- cand[file.exists(file.path(cand, "DESCRIPTION"))]
  if (!length(hit)) stop("cannot locate PhysioAnalysis package root")
  normalizePath(hit[1])
})
golden_dir <- file.path(.pkg_root, "tests", "testthat", "_golden")

# Minimal writer (mirrors helper-golden.R::write_golden but without testthat),
# with a FIXED capture timestamp so the .dcf is deterministic too.
write_golden <- function(value, key, source, tol = 1e-8, dir = golden_dir) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  saveRDS(value, file.path(dir, paste0(key, ".rds")))
  write.dcf(
    data.frame(
      Key = key, Source = source, Tolerance = format(tol),
      Captured = "fixed (deterministic generator)",
      R = as.character(getRversion()),
      stringsAsFactors = FALSE
    ),
    file.path(dir, paste0(key, ".dcf"))
  )
  message(sprintf("wrote golden '%s' (tol=%s)", key, format(tol)))
  invisible(value)
}

# ---------------------------------------------------------------------------
# 1. Network metrics vs igraph
#    Deterministic random undirected binary graph.
# ---------------------------------------------------------------------------
set.seed(101)
n_net <- 12
A <- matrix(0, n_net, n_net)
ut <- upper.tri(A)
A[ut] <- rbinom(sum(ut), 1, 0.35)
A <- A + t(A)
diag(A) <- 0

g <- igraph::graph_from_adjacency_matrix(A, mode = "undirected", diag = FALSE)

# Local clustering coefficient: igraph transitivity(type="local") with
# isolates="zero" matches the package's (2*triangles)/(k(k-1)) definition.
cc_ref <- igraph::transitivity(g, type = "local", isolates = "zero")
write_golden(
  cc_ref, "net_clustering_igraph",
  source = paste0("igraph::transitivity(type='local', isolates='zero') on a ",
                  "seed-101 rbinom(0.35) 12-node undirected binary graph"),
  tol = 1e-8
)

# Betweenness centrality: igraph betweenness(normalized=TRUE). For an
# undirected graph the package's Brandes sum over ordered pairs / ((n-1)(n-2))
# equals igraph's unordered sum / ((n-1)(n-2)/2); both are matched here with
# normalized=TRUE. (The package's normalized=FALSE differs by a factor of 2 --
# it counts ordered pairs -- so the golden uses normalized=TRUE.)
bc_ref <- igraph::betweenness(g, directed = FALSE, normalized = TRUE)
write_golden(
  bc_ref, "net_betweenness_igraph",
  source = paste0("igraph::betweenness(directed=FALSE, normalized=TRUE) on the ",
                  "same seed-101 12-node undirected binary graph"),
  tol = 1e-8
)

# ---------------------------------------------------------------------------
# 2. Correlation connectivity vs base cor()
# ---------------------------------------------------------------------------
set.seed(202)
n_time_c <- 500L
n_ch_c <- 5L
Xc <- matrix(stats::rnorm(n_time_c * n_ch_c), n_time_c, n_ch_c)
colnames(Xc) <- paste0("Ch", seq_len(n_ch_c))

cor_pearson_ref <- stats::cor(Xc, method = "pearson")
dimnames(cor_pearson_ref) <- list(colnames(Xc), colnames(Xc))
write_golden(
  cor_pearson_ref, "conn_correlation_pearson_base",
  source = "base stats::cor(method='pearson') on seed-202 rnorm 500x5 matrix",
  tol = 1e-10
)

cor_spearman_ref <- stats::cor(Xc, method = "spearman")
dimnames(cor_spearman_ref) <- list(colnames(Xc), colnames(Xc))
write_golden(
  cor_spearman_ref, "conn_correlation_spearman_base",
  source = "base stats::cor(method='spearman') on seed-202 rnorm 500x5 matrix",
  tol = 1e-10
)

# ---------------------------------------------------------------------------
# 3. Spectrogram vs a hand-computed one-sided FFT PSD (single window).
#    Reference derived from the documented PSD formula using base fft():
#      psd = |FFT(seg*w)|^2 / (sr * sum(w^2)); non-DC/non-Nyquist bins doubled.
# ---------------------------------------------------------------------------
set.seed(404)
sr_s <- 256
Nwin <- 256L
sig_s <- stats::rnorm(Nwin)

# hanning window matching .getWindow("hanning", 256): 0.5*(1-cos(2*pi*t)),
# t = seq(0,1,length.out = n)
w_han <- 0.5 * (1 - cos(2 * pi * seq(0, 1, length.out = Nwin)))
wp <- sum(w_han^2)
seg <- sig_s * w_han
ft <- stats::fft(seg)
nf <- floor(Nwin / 2) + 1
psd_ref <- Mod(ft[seq_len(nf)])^2 / (sr_s * wp)
psd_ref[2:(nf - 1)] <- 2 * psd_ref[2:(nf - 1)]   # one-sided
write_golden(
  psd_ref, "spectrogram_onesided_psd_handfft",
  source = paste0("hand-computed one-sided FFT PSD |fft(seg*hanning)|^2/",
                  "(sr*sum(w^2)) (non-DC/Nyquist doubled), seed-404 rnorm(256), ",
                  "sr=256, single 256-sample window"),
  tol = 1e-9
)

# ---------------------------------------------------------------------------
# 4. bandPower (Welch) vs an independent from-scratch Welch PSD.
#    Same documented normalization, re-derived with base fft() (not the package
#    kernel): hanning window, 50% overlap, one-sided, averaged over segments.
# ---------------------------------------------------------------------------
set.seed(505)
sr_b <- 256
Nb <- 512L
sig_b <- stats::rnorm(Nb)
bands <- list(theta = c(4, 8), alpha = c(8, 13), beta = c(13, 30))

welch_ref <- function(x, sr, nperseg = 256L) {
  n <- length(x)
  noverlap <- floor(nperseg / 2)
  step <- nperseg - noverlap
  nseg <- floor((n - noverlap) / step)
  nf <- floor(nperseg / 2) + 1
  w <- 0.5 * (1 - cos(2 * pi * seq(0, 1, length.out = nperseg)))
  wsum <- sum(w^2)
  psd <- numeric(nf)
  for (i in seq_len(nseg)) {
    s0 <- (i - 1) * step + 1
    e0 <- s0 + nperseg - 1
    if (e0 > n) break
    ft <- stats::fft(x[s0:e0] * w)
    psd <- psd + Mod(ft[seq_len(nf)])^2
  }
  psd <- psd / (nseg * sr * wsum)
  psd[2:(nf - 1)] <- 2 * psd[2:(nf - 1)]
  list(power = psd, freq = seq(0, sr / 2, length.out = nf))
}

wr <- welch_ref(sig_b, sr_b)
bandpower_ref <- vapply(bands, function(b) {
  idx <- which(wr$freq >= b[1] & wr$freq <= b[2])
  sum(wr$power[idx])
}, numeric(1))
names(bandpower_ref) <- names(bands)
write_golden(
  bandpower_ref, "bandpower_welch_handfft",
  source = paste0("independent from-scratch Welch PSD (base fft, hanning, 50% ",
                  "overlap, one-sided), band sums for theta/alpha/beta; ",
                  "seed-505 rnorm(512), sr=256"),
  tol = 1e-9
)

# ---------------------------------------------------------------------------
# 5. SPM t-map vs base t.test per time point (pooled/equal-variance).
# ---------------------------------------------------------------------------
set.seed(303)
nt <- 60L
n1 <- 8L
n2 <- 9L
g1 <- matrix(stats::rnorm(nt * n1), nt, n1)
g2 <- matrix(stats::rnorm(nt * n2), nt, n2) + 0.3
dat_two <- cbind(g1, g2)

# Two-sample pooled t at each time point via base t.test(var.equal=TRUE).
t_two_ref <- vapply(seq_len(nt), function(i) {
  unname(stats::t.test(dat_two[i, seq_len(n1)],
                       dat_two[i, n1 + seq_len(n2)],
                       var.equal = TRUE)$statistic)
}, numeric(1))
write_golden(
  t_two_ref, "spm_ttest_twosample_base",
  source = paste0("base stats::t.test(var.equal=TRUE) per time point, ",
                  "two-sample; seed-303 groups n1=8 vs n2=9 (+0.3), 60 points"),
  tol = 1e-8
)

# One-sample t against zero at each time point via base t.test.
t_one_ref <- vapply(seq_len(nt), function(i) {
  unname(stats::t.test(dat_two[i, seq_len(n1)])$statistic)
}, numeric(1))
write_golden(
  t_one_ref, "spm_ttest_onesample_base",
  source = paste0("base stats::t.test() per time point, one-sample vs 0; ",
                  "seed-303 group1 n=8, 60 points"),
  tol = 1e-8
)

message("PhysioAnalysis golden fixtures written to: ", golden_dir)
