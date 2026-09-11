# Statistics surfaced from PhysioCore into the PhysioAnalysis statistics layer
# (single source of truth). The implementations live in PhysioCore; these
# re-exports make them appear in the PhysioAnalysis reference so users find the
# reliability / agreement, functional-PCA and circular-statistics tools where
# the other statistical methods live. PhysioAnalysis already imports PhysioCore,
# so these add no new dependency.

# --- Reliability / agreement (scalar) ---
#' @importFrom PhysioCore icc
#' @export
PhysioCore::icc

#' @importFrom PhysioCore sem
#' @export
PhysioCore::sem

#' @importFrom PhysioCore mdc
#' @export
PhysioCore::mdc

#' @importFrom PhysioCore blandAltman
#' @export
PhysioCore::blandAltman

#' @importFrom PhysioCore cohensD
#' @export
PhysioCore::cohensD

#' @importFrom PhysioCore etaSquared
#' @export
PhysioCore::etaSquared

# --- Reliability / agreement (waveform) ---
#' @importFrom PhysioCore waveformCMC
#' @export
PhysioCore::waveformCMC

#' @importFrom PhysioCore waveformICC
#' @export
PhysioCore::waveformICC

#' @importFrom PhysioCore waveformReliability
#' @export
PhysioCore::waveformReliability

# --- Functional PCA ---
#' @importFrom PhysioCore fPCA
#' @export
PhysioCore::fPCA

#' @importFrom PhysioCore reconstructFPCA
#' @export
PhysioCore::reconstructFPCA

#' @importFrom PhysioCore registerCurves
#' @export
PhysioCore::registerCurves

# --- Circular statistics ---
#' @importFrom PhysioCore circularSummary
#' @export
PhysioCore::circularSummary

#' @importFrom PhysioCore rayleighTest
#' @export
PhysioCore::rayleighTest

#' @importFrom PhysioCore watsonWilliamsTest
#' @export
PhysioCore::watsonWilliamsTest

#' @importFrom PhysioCore circularLinearCorrelation
#' @export
PhysioCore::circularLinearCorrelation
