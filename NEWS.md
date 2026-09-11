# PhysioAnalysis 0.5.29

## New features

* `leveneTest()` — Levene's test for homogeneity of variance across `k >= 2` groups, as a one-way ANOVA on
  the absolute deviations from each group's centre; `center = "median"` (default) is the robust
  Brown-Forsythe variant. Returns the F-statistic, p-value, degrees of freedom, `k`, `n` and `center`.
  Reproduces `scipy.stats.levene` and `car::leveneTest` to machine precision. (Its name collides with
  `car::leveneTest` — call `PhysioAnalysis::leveneTest` when `car` is attached.)

# PhysioAnalysis 0.5.28

## New features

* `intraclassCorrelation()` — the six McGraw & Wong (1996) intraclass correlation coefficients (ICC1/2/3
  single and ICC1k/2k/3k average) from a subjects-by-raters matrix, via a two-way ANOVA. The standard
  test-retest / inter-rater reliability index. Reproduces `psych::ICC(lmer = FALSE)` and
  `pingouin.intraclass_corr` to machine precision (note psych's default `lmer = TRUE` uses REML variance
  components, which differ slightly).

# PhysioAnalysis 0.5.27

## New features

* `cronbachAlpha()` — Cronbach's alpha, the standard index of a scale's internal-consistency reliability,
  computed from the item and total-score variances. Takes a respondents-by-items matrix and returns the raw
  alpha, the mean inter-item correlation, and the item / respondent counts. Reproduces base R, the `psych`
  package (`psych::alpha` raw_alpha) bit-for-bit and `pingouin.cronbach_alpha` to machine precision. Adds a
  psychometrics capability to the analysis toolkit.

# PhysioAnalysis 0.5.26

## New features

* `multipleRegression()` — multivariable OLS regression of a response (first column) on two or more
  predictors: the coefficient table (estimate, SE, t, p) plus R-squared, adjusted R-squared and the overall
  F-test, via the same QR decomposition as base R. The many-predictor generalisation of
  `linearRegression()`. Reproduces `stats::lm` / `summary.lm` bit-for-bit and `statsmodels` OLS to machine
  precision.

# PhysioAnalysis 0.5.25

## New features

* `nemenyiTest()` — the Nemenyi post-hoc test, all pairwise comparisons of the treatments after a
  Friedman repeated-measures test (the repeated-measures counterpart of `tukeyHSD()` / `dunnTest()`). From
  the within-block mean ranks and the studentized-range distribution with single-step family-wise error
  control; takes a blocks-by-treatments matrix and returns a data frame of comparisons (`statistic`,
  `p_adj`). Reproduces base R's Nemenyi formula, the `PMCMRplus` package and
  `scikit_posthocs.posthoc_nemenyi_friedman` to machine precision.

# PhysioAnalysis 0.5.24

## New features

* `dunnTest()` — Dunn's test, the non-parametric post-hoc for all pairwise comparisons after a
  Kruskal-Wallis test (the rank-based counterpart of `tukeyHSD()`). Uses the shared tie-corrected rank
  variance and a chosen multiple-comparison adjustment (bonferroni default; also none/holm/BH), and
  returns a data frame of comparisons (`z`, `p_value`, `p_adj`). Reproduces base R's Dunn formula, the
  `PMCMRplus` package and `scikit_posthocs.posthoc_dunn` to machine precision.

# PhysioAnalysis 0.5.23

## New features

* `tukeyHSD()` — Tukey's Honest Significant Difference post-hoc test: all pairwise mean comparisons after
  a one-way ANOVA, with studentized-range p-values, confidence intervals and single-step control of the
  family-wise error rate. Accepts a named list of group vectors and returns a data frame of comparisons
  (`diff`, `lwr`, `upr`, `p_adj`). Reproduces base R `stats::TukeyHSD` and `scipy.stats.tukey_hsd` to
  machine precision.

# PhysioAnalysis 0.5.22

## New features

* `partialCorrelation()` — partial correlation of `x` and `y` controlling for one or more covariates
  `z` (Pearson, or `method = "spearman"` for a rank/robust partial correlation), by residualising both
  on `[1, z]` via QR and correlating the residuals. The standard test of whether an association survives
  adjustment for a confound. Returns the partial correlation, its two-sided p-value, `df` (`n - 2 - k`),
  `n` and `k`. Reproduces base R's residual method and the `ppcor` package to machine precision.

# PhysioAnalysis 0.5.21

## New features

* `kendallTau()` — Kendall's tau-b rank correlation (concordance-based, tie-corrected), the third
  rank-correlation method alongside Pearson's r and Spearman's rho (`correlationTest`) and the robust
  choice for ordinal or heavily-tied data and small samples. Accepts `(x, y)` vectors or a two-column
  `x`, and returns the tau-b statistic, the two-sided normal-approximation p-value, the z-statistic and
  `n`. Reproduces base R `stats::cor.test(method = "kendall")` and `scipy.stats.kendalltau` bit-for-bit.

# PhysioAnalysis 0.5.20

## New features

* `linearRegression()` — simple univariate OLS regression: slope, intercept, R-squared and the two-sided
  t-test of the slope. Accepts `(x, y)` vectors or a two-column `x`, and drops incomplete pairs. The slope,
  R-squared and slope p-value reproduce `scipy.stats.linregress` bit-for-bit and base R `stats::lm` to
  machine precision.

# PhysioAnalysis 0.5.19

## New features

* `kruskalTest()` — the Kruskal-Wallis rank-sum test, the non-parametric between-groups omnibus and the
  distribution-free counterpart of `oneWayAnova()` (no normality/equal-variance assumption). Accepts a
  list of `k >= 2` group vectors (or the groups as arguments) and returns the tie-corrected H-statistic,
  p-value, degrees of freedom (`k - 1`), `k` and `n`. The H is accumulated in the same order as base R,
  so it reproduces `stats::kruskal.test` bit-for-bit and `scipy.stats.kruskal` to machine precision.

# PhysioAnalysis 0.5.18

## New features

* `friedmanTest()` — the Friedman rank-sum test, the non-parametric repeated-measures omnibus and the
  within-subject counterpart of `oneWayAnova()`. Accepts a blocks-by-treatments matrix (or a list of
  aligned treatment vectors), drops incomplete blocks, and returns the tie-corrected chi-square
  statistic, p-value, degrees of freedom (`k - 1`), the number of treatments (`k`) and complete blocks
  (`n`). Reproduces base R `stats::friedman.test` bit-for-bit and `scipy.stats.friedmanchisquare` to
  machine precision.

# PhysioAnalysis 0.5.17

## New features

* `oneWayAnova()` — the classic equal-variance one-way ANOVA (omnibus F-test), the many-group
  generalisation of `twoSampleTTest()`. Accepts a list of `k >= 2` numeric group vectors (or the
  groups as separate arguments) and returns the F-statistic, p-value, between/within degrees of
  freedom (`df1`, `df2`), the number of groups (`k`) and the total sample size (`n`). The omnibus F
  reproduces base R `stats::oneway.test(var.equal = TRUE)` bit-for-bit and `scipy.stats.f_oneway`
  to machine precision.

# PhysioAnalysis 0.5.16

## New features

* `twoSampleTTest()` — the independent two-sample t-test, Welch's (unequal-variance) by default or
  Student's pooled with `var_equal = TRUE`. Accepts two vectors or a two-element list; returns the
  statistic, two-sided p-value and df. Welch reproduces `stats::t.test` and
  `scipy.stats.ttest_ind(equal_var = FALSE)` bit-for-bit.

# PhysioAnalysis 0.5.15

## New features

* `correlationTest()` — the correlation test (Pearson, default, or Spearman) with the two-sided t-based
  p-value and df. Accepts two vectors or a matrix; returns the correlation coefficient, p-value and df.
  Pearson reproduces `stats::cor.test` and `scipy.stats.pearsonr` bit-for-bit.

# PhysioAnalysis 0.5.14

## New features

* `pairedTTest()` — the paired-samples t-test (one-sample t of the within-pair differences `x - y`), the
  standard before/after or condition-A/B test on matched subjects. Accepts two vectors or a two-column
  matrix; returns the statistic, two-sided p-value and df. Reproduces `stats::t.test(x, y, paired = TRUE)`
  and `scipy.stats.ttest_rel` bit-for-bit.

# PhysioAnalysis 0.5.13

## New features

* `ksNormalityTest()` — the Kolmogorov-Smirnov one-sample normality test, `D = sup|ECDF - fitted normal
  CDF|` with the two-sided asymptotic p-value. The empirical-CDF (distribution-shape) counterpart of the
  moment-based `jarqueBeraTest()`. The `D` statistic reproduces `stats::ks.test` and
  `scipy.stats.kstest` bit-for-bit. On a long, mildly non-Gaussian signal the two tests can disagree
  (KS reads the overall shape; Jarque-Bera flags the skew/kurtosis with large-N power). Caveat: with
  parameters estimated from the sample the naive p-value is anti-conservative (Lilliefors); the `D`
  statistic is exact.

# PhysioAnalysis 0.5.12

## New features

* `trimmedMean()` — the trimmed mean (arithmetic mean after discarding the `floor(n*trim)` most extreme
  values from each end), a robust location estimator completing the robust toolkit alongside
  `medianAbsDev()` and `interquartileRange()` (robust dispersion) and `signalMoments()` (classical
  mean/SD). Reproduces `scipy.stats.trim_mean` and base R `mean(x, trim)` bit-for-bit. On a skewed
  distribution it lies between the mean and the median, so reporting the mean, a trimmed mean and the
  median together shows how much the tail inflates the mean (certified on real Fantasia RR intervals).

# PhysioAnalysis 0.5.11

## New features

* `interquartileRange()` — the interquartile range (Q3 − Q1, type-7 quantiles), a robust, quantile-based
  dispersion measure completing the robust-dispersion toolkit alongside `medianAbsDev()` (deviation-based
  scale) and `signalMoments()` (classical SD). Reproduces `scipy.stats.iqr` and base R `IQR` bit-for-bit.
  Its normal-consistent scale `IQR/(2*qnorm(0.75))` is a robust SD estimate; comparing it with the SD is
  a distribution-shape / non-normality check (certified on real Fantasia RR intervals, where the skewed
  RR distribution makes the robust scale fall below SDNN).

# PhysioAnalysis 0.5.10

## New features

* `medianAbsDev()` — the median absolute deviation, a robust dispersion measure and the robust
  companion of `signalMoments()`'s classical SD. `MAD = constant * median(|x - median(x)|)`; the
  default `constant = 1` returns the raw MAD (reproducing `scipy.stats.median_abs_deviation` and base R
  `mad(constant = 1)` bit-for-bit), and `constant = 1 / qnorm(0.75)` (~1.4826) gives the
  normal-consistent robust estimate of the SD. Certified against both scipy and base R on real eegmmidb
  POz EEG; comparing its normal-consistent scale with the SD is a lightweight artifact check.

# PhysioAnalysis 0.5.9

## New features

* `jarqueBeraTest()` — the Jarque-Bera goodness-of-fit test for normality (Jarque & Bera 1980):
  `JB = N/6 (S^2 + K^2/4)` from the skewness S and excess kurtosis K (via `signalMoments()`), compared
  to a chi-squared(2) under the normality null. The formal, moment-based counterpart of the descriptive
  near-Gaussian check; returns the statistic, p-value and df. Same definition as
  `scipy.stats.jarque_bera`, reproduced bit-for-bit for both the statistic and the p-value on real
  eegmmidb POz EEG.

# PhysioAnalysis 0.5.8

## New features

* `signalMoments()` — the first four moments of a signal's amplitude distribution: mean, (population)
  standard deviation, skewness (asymmetry) and excess kurtosis (tailedness vs Gaussian). Standard
  shape descriptors (e.g. EEG artifact flagging via high kurtosis). Definitions match
  `scipy.stats.skew` (biased) and `scipy.stats.kurtosis` (Fisher / excess, biased) with the population
  SD (ddof = 0), reproduced bit-for-bit (max |diff| ~4e-16) on real eegmmidb POz EEG.

# PhysioAnalysis 0.5.7

## New features

* `arSpectrum()` — the autoregressive (parametric) spectral density: `S(f) = var / |1 - sum phi_k
  exp(-i 2pi f k)|^2` from a Yule-Walker AR(`order`) fit (`arYuleWalker()`), evaluated on `n_freq`
  normalized frequencies. A smooth, low-variance, sharp-peaked parametric alternative to the periodogram
  / Welch PSD. Matches `stats::spec.ar` (including its N/(N-p-1) small-sample variance convention)
  bit-for-bit on real eegmmidb POz EEG; its alpha peak agrees with the nonparametric periodogram within
  ~0.5 Hz.

# PhysioAnalysis 0.5.6

## New features

* `kpssTest()` — the KPSS test statistic for level stationarity (Kwiatkowski et al. 1992):
  `eta = N^-2 sum S_t^2 / lrv` from the partial sums of the demeaned signal, with a Bartlett
  (Newey-West) long-run variance (truncation `lag`). The complement of `adfTest()` — KPSS's null is
  stationarity, so a *small* statistic fails to reject it. The statistic reproduces
  `statsmodels.tsa.stattools.kpss(regression = "c")` bit-for-bit on real eegmmidb POz EEG. Only the
  statistic is returned; critical values / p-value come from the reference tables and are not computed.

# PhysioAnalysis 0.5.5

## New features

* `adfTest()` — the Augmented Dickey-Fuller unit-root test statistic (Dickey & Fuller 1979): the OLS
  t-statistic on the lagged level in `d y_t = a + b y_{t-1} + sum g_i d y_{t-i} + e`; a large negative
  value rejects the unit-root (non-stationarity) null. Adds the stationarity check the time-series
  methods (ACF/AR/PSD) implicitly assume. The statistic reproduces
  `statsmodels.tsa.stattools.adfuller(autolag = None)` bit-for-bit on real eegmmidb POz EEG. Only the
  statistic is returned (an exact regression quantity); critical values / p-value require the Dickey-
  Fuller / MacKinnon reference tables and are not computed.

# PhysioAnalysis 0.5.4

## New features

* `ljungBoxTest()` — the Ljung-Box portmanteau test for autocorrelation (Ljung & Box 1978):
  `Q = N(N+2) sum_{k=1}^{h} rho_k^2 / (N-k)` from the ACF (`autocorrelation()`), with a chi-squared(`lag`)
  reference under the white-noise null; returns the statistic, p-value and df. The standard formal test
  for autocorrelation (raw signal) and model adequacy (residuals). Same definition as
  `statsmodels.stats.diagnostic.acorr_ljungbox` and `stats::Box.test(type = "Ljung-Box")`, which it
  reproduces bit-for-bit for both the statistic and the p-value on real eegmmidb POz EEG.

# PhysioAnalysis 0.5.3

## New features

* `arYuleWalker()` — fits an autoregressive AR(`order`) model by the Yule-Walker method (Durbin-Levinson
  on the biased demeaned ACF, building on `autocorrelation()`); returns the AR coefficients and the
  innovation (one-step prediction) variance, the same definition as
  `statsmodels.regression.linear_model.yule_walker(method = "mle")`. Completes the ACF/PACF/AR-fit lane
  (choose the order from the `partialAutocorrelation` cut-off). Validated against
  `statsmodels.regression.yule_walker(method = "mle")` bit-for-bit for both the coefficients (max |diff|
  ~5e-15) and the innovation variance (~5e-13) on real eegmmidb POz EEG.

# PhysioAnalysis 0.5.2

## New features

* `partialAutocorrelation()` — the single-series partial autocorrelation function (PACF) at lags
  0..`lag_max`, via the Durbin-Levinson recursion on the biased demeaned ACF (building on
  `autocorrelation()`); the same definition as `stats::pacf` and `statsmodels.tsa.pacf(method = "ldb")`.
  Where the ACF decays, the PACF cuts off after the AR order — the standard AR model-order tool.
  Validated against `statsmodels.tsa.pacf(method = "ldb")` bit-for-bit (max |diff| ~6e-15) on real
  eegmmidb POz EEG.

# PhysioAnalysis 0.5.1

## New features

* `autocorrelation()` — the single-series sample autocorrelation function (ACF) at lags 0..`lag_max`,
  using the biased, demeaned estimator (the same definition as `stats::acf` and
  `statsmodels.tsa.acf(adjusted = FALSE)`). The fundamental tool for periodicity, rhythmicity and the
  decorrelation time of a signal; complements the between-two-series `crossCorrelation()`. Validated
  against `statsmodels.tsa.acf` bit-for-bit (max |diff| ~2e-16) on real eegmmidb POz EEG.

# PhysioAnalysis 0.5.0

Surfaced the reliability/agreement, functional-PCA and circular-statistics tools
into the statistics layer. These validated methods live in PhysioCore (their
single source of truth) and are now re-exported by PhysioAnalysis so they appear
in its reference alongside the other statistical methods (no new dependency —
PhysioAnalysis already imports PhysioCore):

* Reliability / agreement: `icc()`, `sem()`, `mdc()`, `blandAltman()`,
  `cohensD()`, `etaSquared()` (scalar) and `waveformCMC()`, `waveformICC()`,
  `waveformReliability()` (pointwise waveform reliability).
* Functional PCA: `fPCA()`, `reconstructFPCA()`, `registerCurves()`.
* Circular statistics: `circularSummary()`, `rayleighTest()`,
  `watsonWilliamsTest()`, `circularLinearCorrelation()`.

(The fPCA/circular/waveform implementations were relocated from PhysioMoCap down
into PhysioCore, and re-exported back into PhysioMoCap for back-compatibility.)

# PhysioAnalysis 0.3.5

## Performance

- `epochData()` extracts fixed-length epochs through a compiled kernel
  (`cpp_epoch_fixed`, OpenMP across epochs) that copies each epoch's samples
  directly, instead of allocating an NA-initialized 4D array and gathering
  per channel. The 2D input is no longer reshaped/flattened first. Results are
  identical to the pure-R gather (retained as the internal `.epochFixedR()`
  reference). On a 64-channel, 5-minute, 256 Hz recording with 150 events the
  standard epoching drops from ~0.12 s to ~0.035 s. Variable-length epoching
  (when `tmax` names an end-event type) is unchanged.

# PhysioAnalysis 0.3.4

## Performance

- `spectrogram()` now computes the short-time Fourier transform in a compiled
  radix-2 kernel (OpenMP-parallel across analysis windows) for power-of-two
  window sizes, with a batched `mvfft` path for other sizes. The pure-R
  per-window loop is retained internally as `.spectrogramLoop()` and covered by
  numerical-parity tests (agreement ~1e-9 or better). On a 64-channel, 5-minute,
  256 Hz recording the all-channel spectrogram drops from ~0.55 s to ~0.26 s.
- `fftSignals()` batches the per-channel FFT through a single `mvfft` call
  instead of a per-column `apply()`, removing the array-permute overhead;
  results are identical. The same recording drops from ~0.39 s to ~0.26 s.

## Bug fixes

- `spectrogram()` now raises a clear error when the requested `overlap` implies
  a hop of less than one sample (e.g. `overlap = 1`), instead of dividing by
  zero. The compiled STFT kernel additionally guards `step >= 1` and returns an
  empty result for signals shorter than one window.
- `fftSignals()` handles 3D single-channel arrays (`time x 1 x samples`)
  correctly (the mvfft slice is kept a matrix).
- The internal `.spectrogramLoop()` reference no longer mis-doubles the DC/
  Nyquist bin for windows of two or fewer frequency bins.

# PhysioAnalysis 0.3.3

## Validation

- Added a numeric EEG<->MNE parity test (`test-eeg-mne-parity.R`, fixture
  `tests/testthat/fixtures/eeg-mne-reference.rds`, generated by
  `data-raw/eeg_mne_reference.{py,R}`), replacing the manuscript's
  documentation-level "MNE-equivalent" correspondence (VAL-04). Welch band power
  (theta/alpha/beta/gamma) matches MNE `psd_array_welch` at correlation 1.0 and
  < 0.5% per channel, reproducing the manuscript's alpha-power cross-validation.
  Documented findings: the delta band diverges because MNE detrends each Welch
  segment while `.welchPSD` does not; and `plv()`/`wPLI()` (continuous Hilbert
  estimators) agree with `mne_connectivity`'s epoch-spectral estimator on
  coupling structure (~0.95 correlation) but are not identical, so the
  manuscript's PLV=0.9996 / wPLI=0.9994 figures are not reproduced.

# PhysioAnalysis 0.3.2

## Bug Fixes

- Corrected the SPM cluster-extent p-value (`.rftClusterPValue`, used by
  `spmTTest()`/`spmPairedTTest()`/`spmAnova()`). The expected cluster extent is
  `E[k] = rho0/rho1` (Worsley 1996 / spm12); the previous form used
  `(E[m]/E[n])^2` for the scale, adding a spurious `1/resel` term that diverged
  from `spm1d`/`rft1d`. The corrected formula now reproduces `rft1d.p_cluster`
  to machine precision (validated on a grid in `test-spm1d-parity.R`, VAL-07).
  Cluster p-values returned by the SPM functions change accordingly. Note that
  `spm1d`'s *reported* cluster P additionally evaluates the RFT cluster
  probability at an idiosyncratic internal cluster height rather than the
  cluster-defining threshold, so end-to-end cluster P still differs from
  `spm1d`; PhysioAnalysis uses the standard threshold-based convention.

# PhysioAnalysis 0.3.1

## Validation

- Added a numerical parity test for the SPM/random-field-theory functions
  against an actual `spm1d` 0.4.53 + `rft1d` 0.2.5 run (`test-spm1d-parity.R`,
  fixture `tests/testthat/fixtures/spm1d-reference.rds`, generated by
  `data-raw/spm1d_reference.{py,R}`), replacing the previous closed-form Worsley
  surrogate as the reference. The RFT critical t/F thresholds match `rft1d` to
  < 0.7%, the t/F statistic fields to ~1e-13, and the residual FWHM exactly.
  Two divergences are documented and tracked as follow-up (VAL-07): the
  cluster-extent p-values (Friston-1994 GRF approximation vs `spm1d`'s
  field-specific distribution) and `spm1d`'s `anova1` F-threshold convention.

# PhysioAnalysis 0.3.0

- Added the offline inference parity gate for SPM, regression, MANOVA,
  permutation fields, and parametric and rank-based effect sizes. Committed
  references are provenance- and SHA-256-checked before public API execution.
- Added optional Perrin spherical-spline interpolation to `plotTopomap()` and
  `plotTopomapSeries()`. The existing Shepard inverse-distance-weighted method
  remains the default and retains its numerical behavior.
- Added explicit spline stiffness, Legendre term-count, and diagonal
  regularization controls. Spherical splines are documented as spatial
  interpolation, not as a surface-Laplacian or source-localization method.

# PhysioAnalysis 0.2.0

Initial release of PhysioAnalysis as a standalone package in the
x-biosignal ecosystem, split out from the original PhysioExperiment
monolith. PhysioAnalysis provides analysis and visualization for
physiological signal data held in `PhysioExperiment` objects (built on
the Bioconductor `SummarizedExperiment` container, via PhysioCore).

## New Features

- Spectral and time-frequency analysis: `fftSignals()` for the fast
  Fourier transform, `spectrogram()` (STFT with one-sided PSD
  normalization), `waveletTransform()`, and `bandPower()` for
  band-limited power extraction.
- Analytic-signal utilities via the Hilbert transform:
  `hilbertTransform()`, `instantaneousAmplitude()`, and
  `instantaneousPhase()`.
- Epoching and averaging: `epochData()` cuts event-locked epochs,
  `averageEpochs()` and `grandAverage()` build condition and group
  averages, `epochTimes()` returns epoch time vectors, and
  `epochSliding()` produces overlapping sliding-window epochs.
- Functional connectivity: magnitude-squared `coherence()` (Welch's
  method) and `crossSpectrum()`, phase-based synchrony via `plv()`,
  `pli()`, and `wPLI()`, plus `correlationMatrix()` and a unified
  `connectivityMatrix()` interface.
- Graph-theoretic network analysis with an Rcpp/RcppArmadillo backend:
  network construction (`adjacencyMatrix()`, `thresholdNetwork()`,
  `binarizeNetwork()`) and topology measures including `nodeDegree()`,
  `clusteringCoefficient()`, `pathLength()`, `betweennessCentrality()`,
  `eigenvectorCentrality()`, `globalEfficiency()`, `localEfficiency()`,
  `smallWorldness()`, and `modularity()`.
- Spectral graph methods: `graphLaplacian()`, `spectralDecomposition()`,
  and `spectralClustering()`.
- Dynamic connectivity over time via `slidingWindowConnectivity()` and
  `temporalStability()`.

## Statistical Analysis

- Epoch-level testing: pointwise `tTestEpochs()` and `anovaEpochs()`
  across time points and channels, with `effectSize()` (Cohen's d),
  `bootstrapCI()`, `correctPValues()`, and `findSignificantWindows()`.
- Cluster-based permutation testing (`clusterPermutationTest()`) for
  family-wise error control over EEG/MEG-style multichannel data
  (Maris & Oostenveld methodology).
- Statistical Parametric Mapping for continuous 1-D waveforms
  (`spmTTest()`, `spmPairedTTest()`, `spmAnova()`, `plotSPM()`) using
  Random Field Theory correction, suited to biomechanical gait curves.

## Visualization

- Signal and evoked-response plots: `plotSignal()`, `plotMultiChannel()`,
  `plotERP()`, and `plotPSD()`.
- Scalp topographic maps: `plotTopomap()` and `plotTopomapSeries()`,
  using Shepard inverse-distance-weighted interpolation.
- Time-frequency and network displays: `plotSpectrogram()`,
  `plotNetwork()`, `plotAdjacencyMatrix()`, `plotNetworkMetrics()`,
  `plotDynamicConnectivity()`, and `plotNetworkStability()`.

## Robustness

- Graph-metric entry points now validate that they receive a square
  (n x n) adjacency matrix and fail with a clear error, preventing the
  C++ backend from indexing out of bounds when raw time-by-channel data
  is passed by mistake.
