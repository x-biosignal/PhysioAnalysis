# Statistics surfaced from PhysioExperiment into the PhysioAnalysis statistics layer
# (single source of truth). The implementations live in PhysioExperiment; these
# re-exports make them appear in the PhysioAnalysis reference so users find the
# reliability / agreement, functional-PCA and circular-statistics tools where
# the other statistical methods live. PhysioAnalysis already imports PhysioExperiment,
# so these add no new dependency.

# --- Reliability / agreement (scalar) ---
#' @importFrom PhysioExperiment icc
#' @export
PhysioExperiment::icc

#' @importFrom PhysioExperiment sem
#' @export
PhysioExperiment::sem

#' @importFrom PhysioExperiment mdc
#' @export
PhysioExperiment::mdc

#' @importFrom PhysioExperiment blandAltman
#' @export
PhysioExperiment::blandAltman

#' @importFrom PhysioExperiment cohensD
#' @export
PhysioExperiment::cohensD

#' @importFrom PhysioExperiment etaSquared
#' @export
PhysioExperiment::etaSquared

# --- Reliability / agreement (waveform) ---
#' @importFrom PhysioExperiment waveformCMC
#' @export
PhysioExperiment::waveformCMC

#' @importFrom PhysioExperiment waveformICC
#' @export
PhysioExperiment::waveformICC

#' @importFrom PhysioExperiment waveformReliability
#' @export
PhysioExperiment::waveformReliability

# --- Functional PCA ---
#' @importFrom PhysioExperiment fPCA
#' @export
PhysioExperiment::fPCA

#' @importFrom PhysioExperiment reconstructFPCA
#' @export
PhysioExperiment::reconstructFPCA

#' @importFrom PhysioExperiment registerCurves
#' @export
PhysioExperiment::registerCurves

# --- Circular statistics ---
#' @importFrom PhysioExperiment circularSummary
#' @export
PhysioExperiment::circularSummary

#' @importFrom PhysioExperiment rayleighTest
#' @export
PhysioExperiment::rayleighTest

#' @importFrom PhysioExperiment watsonWilliamsTest
#' @export
PhysioExperiment::watsonWilliamsTest

#' @importFrom PhysioExperiment circularLinearCorrelation
#' @export
PhysioExperiment::circularLinearCorrelation
