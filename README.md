# PhysioAnalysis <img src="man/figures/logo.png" align="right" height="139" alt="PhysioAnalysis logo" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/x-biosignal/PhysioAnalysis/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/x-biosignal/PhysioAnalysis/actions/workflows/R-CMD-check.yaml)
[![CRAN status](https://www.r-pkg.org/badges/version/PhysioAnalysis)](https://CRAN.R-project.org/package=PhysioAnalysis)
[![r-universe](https://x-biosignal.r-universe.dev/badges/PhysioAnalysis)](https://x-biosignal.r-universe.dev/PhysioAnalysis)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

**Analysis and Visualization for PhysioExperiment Objects**

PhysioAnalysis provides 43 exported functions for spectral analysis, time-frequency decomposition, epoching, functional connectivity, statistical testing, and publication-quality visualization of multi-modal physiological signals. Built on PhysioCore, it operates directly on `PhysioExperiment` objects and covers the full analysis workflow from raw epochs to statistical inference and topographic visualization.

## Installation

You can install PhysioAnalysis from [r-universe](https://x-biosignal.r-universe.dev):

```r
install.packages("PhysioAnalysis",
  repos = c("https://x-biosignal.r-universe.dev", "https://cloud.r-project.org"))
```

Or install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("x-biosignal/PhysioAnalysis")
```

## Quick Start

```r
library(PhysioAnalysis)

# Create sample EEG data with events
signal_matrix <- matrix(rnorm(2500 * 32), nrow = 2500, ncol = 32)
pe <- PhysioExperiment(
  assays = list(raw = signal_matrix),
  samplingRate = 250
)
pe <- addEvents(pe, name = "stimulus", onset = c(1.0, 3.0, 5.0, 7.0))

# Epoch around stimulus events (-0.2 to 0.8 s)
pe_epoched <- epochData(pe, event = "stimulus", pre = 0.2, post = 0.8)

# Compute the average ERP
erp <- averageEpochs(pe_epoched)

# Run a cluster permutation test (condition A vs B)
result <- clusterPermutationTest(pe_epoched,
  group = "condition", n_permutations = 1000)

# Visualize results
plotERP(erp, channels = c("Fz", "Cz", "Pz"))
plotTopomap(erp, time = 0.3)
```

## Features

### FFT and Spectral Analysis

Frequency-domain analysis and spectral decomposition:

- `fftSignals()` -- compute FFT magnitude and phase spectra for all channels
- `bandPower()` -- extract power in standard frequency bands (delta, theta, alpha, beta, gamma) or custom ranges
- `hilbertTransform()` -- analytic signal via Hilbert transform for envelope extraction
- `instantaneousAmplitude()`, `instantaneousPhase()` -- extract amplitude envelope and instantaneous phase from analytic signals

### Time-Frequency Analysis

Multi-resolution time-frequency decomposition:

- `spectrogram()` -- short-time Fourier transform (STFT) with configurable window and overlap
- `waveletTransform()` -- continuous wavelet transform (Morlet) for time-frequency representations
- `plotSpectrogram()` -- visualize time-frequency power maps with customizable color scales

### Epoching

Trial segmentation and averaging:

- `epochData()` -- segment continuous data into time-locked epochs around events
- `averageEpochs()` -- compute trial-averaged waveforms (e.g., ERPs, ERFs)
- `grandAverage()` -- compute grand averages across subjects or sessions
- `epochTimes()` -- retrieve the time vector for epoched data

### Connectivity Analysis

Functional and effective connectivity metrics:

- **Spectral coherence:** `coherence()`, `crossSpectrum()` -- frequency-domain measures of linear coupling
- **Phase synchrony:** `plv()` (Phase Locking Value), `pli()` (Phase Lag Index), `wPLI()` (weighted Phase Lag Index) -- volume conduction-robust phase metrics
- **Correlation:** `correlationMatrix()` -- pairwise amplitude correlations across channels
- **General interface:** `connectivityMatrix()` -- compute any connectivity metric as a channel-by-channel matrix

### Statistical Testing

Parametric, non-parametric, and mass-univariate statistical methods:

- **Classical tests:** `tTestEpochs()`, `anovaEpochs()` -- point-by-point or window-based parametric tests across conditions
- **Cluster permutation:** `clusterPermutationTest()` -- non-parametric cluster-based permutation testing for family-wise error control over space and time
- **Effect sizes:** `effectSize()` -- Cohen's d, Hedges' g, and eta-squared
- **Confidence intervals:** `bootstrapCI()` -- bootstrap confidence intervals for any statistic
- **Multiple comparisons:** `correctPValues()` -- FDR (Benjamini-Hochberg) and Bonferroni correction
- **Temporal analysis:** `findSignificantWindows()` -- identify contiguous time windows with significant effects
- **SPM methods:** `spmTTest()`, `spmPairedTTest()`, `spmAnova()` -- Statistical Parametric Mapping for continuous signal analysis

### Visualization

Publication-quality plots for physiological signal data:

- **Time series:** `plotSignal()` -- single-channel waveform display with event markers
- **Multi-channel:** `plotMultiChannel()` -- stacked multi-channel display with vertical offset and scaling
- **ERP plots:** `plotERP()` -- event-related potential waveforms with confidence bands and condition overlays
- **Power spectra:** `plotPSD()` -- power spectral density plots with frequency band shading
- **Topographic maps:** `plotTopomap()` -- scalp topography with channel
  markers and either Shepard inverse-distance weighting (default) or Perrin
  spherical-spline interpolation
- **Topomap series:** `plotTopomapSeries()` -- temporal evolution of scalp topography at multiple time points
- **Spectrogram:** `plotSpectrogram()` -- time-frequency power maps with configurable color scales

## Dependencies

- **R** (>= 4.2)
- **PhysioCore** -- core data structures
- **signal** -- DSP primitives
- **methods**, **SummarizedExperiment**, **S4Vectors**, **stats**, **graphics**, **grDevices**
- **Suggests:** ggplot2, testthat, knitr, rmarkdown

## PhysioExperiment Ecosystem

PhysioAnalysis is part of the PhysioExperiment ecosystem, a suite of R packages for multi-modal physiological signal analysis:

| Package | Description |
|---------|-------------|
| [PhysioCore](https://github.com/x-biosignal/PhysioCore) | Core data structures and accessors |
| [PhysioIO](https://github.com/x-biosignal/PhysioIO) | File I/O (EDF, HDF5, BIDS, CSV, MAT) |
| [PhysioPreprocess](https://github.com/x-biosignal/PhysioPreprocess) | Signal preprocessing and artifact removal |
| **PhysioAnalysis** | Spectral analysis, epoching, statistics, visualization |
| [PhysioMoCap](https://github.com/x-biosignal/PhysioMoCap) | Motion capture data processing |
| [PhysioOpenSim](https://github.com/x-biosignal/PhysioOpenSim) | OpenSim biomechanical modeling integration |

Visit the [r-universe page](https://x-biosignal.r-universe.dev) to browse all available packages.

## License

MIT License. See [LICENSE](LICENSE) for details.

## Author

Yusuke Matsui

## Governance & support

Part of the [Physio ecosystem](https://x-biosignal.r-universe.dev). Community and
policy documents live in the umbrella repository:

- [Code of Conduct](https://github.com/x-biosignal/PhysioExperiment/blob/main/CODE_OF_CONDUCT.md)
- [Contributing](https://github.com/x-biosignal/PhysioExperiment/blob/main/CONTRIBUTING.md)
- [Governance](https://github.com/x-biosignal/PhysioExperiment/blob/main/GOVERNANCE.md)
- [Support](https://github.com/x-biosignal/PhysioExperiment/blob/main/SUPPORT.md)
- [Security policy](https://github.com/x-biosignal/PhysioExperiment/blob/main/SECURITY.md)
- [Deprecation & lifecycle policy](https://github.com/x-biosignal/PhysioExperiment/blob/main/DEPRECATION.md)
