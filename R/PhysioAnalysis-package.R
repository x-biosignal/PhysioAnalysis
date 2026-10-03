#' PhysioAnalysis: Analysis and Visualization for PhysioExperiment Objects
#'
#' Spectral and time-frequency analysis, epoching, functional connectivity,
#' network metrics, statistical testing, and publication-quality plotting for
#' multi-modal physiological signals. Every function operates on a
#' `PhysioExperiment` object (from \pkg{PhysioExperiment}) and, where it produces
#' a result, returns it either as a plain value or back in the object.
#'
#' @section Spectral and time-frequency analysis:
#' [fftSignals()], [bandPower()], [hilbertTransform()],
#' [instantaneousAmplitude()], [instantaneousPhase()], [spectrogram()],
#' [waveletTransform()], [specparam()] (aperiodic/oscillatory parameterization).
#'
#' @section Epoching and averaging:
#' [epochData()], [averageEpochs()], [grandAverage()], [epochTimes()],
#' [epochSliding()].
#'
#' @section Connectivity and networks:
#' Pairwise coupling with [coherence()], [crossSpectrum()], [plv()], [pli()],
#' [wPLI()], [correlationMatrix()] and the general [connectivityMatrix()]; graph
#' metrics with [adjacencyMatrix()], [thresholdNetwork()], [nodeDegree()],
#' [clusteringCoefficient()], [pathLength()], [betweennessCentrality()],
#' [globalEfficiency()], [smallWorldness()], [modularity()].
#'
#' @section Statistical testing:
#' Epoch-wise tests [tTestEpochs()], [anovaEpochs()]; mass-univariate
#' [clusterPermutationTest()], [tfce()], [findSignificantWindows()]; effect sizes
#' [effectSize()], [cohensD()], [rankBiserial()], [cliffsDelta()]; inference
#' helpers [bootstrapCI()], [correctPValues()]. Statistical Parametric Mapping:
#' [spmTTest()], [spmPairedTTest()], [spmAnova()], [spmRegression()],
#' [spmMANOVA()], [spmSnPM()]. Functional-data models [functionalRegression()],
#' [scalarOnFunctionRegression()], [functionalMixedModel()]. Classical tests
#' ([twoSampleTTest()], [oneWayAnova()], [kruskalTest()], [linearRegression()],
#' ...) and time-series diagnostics ([autocorrelation()], [arYuleWalker()],
#' [ljungBoxTest()], [adfTest()], [kpssTest()]) are also exported.
#'
#' @section Visualization:
#' [plotSignal()], [plotMultiChannel()], [plotERP()], [plotPSD()],
#' [plotSpectrogram()], [plotTopomap()], [plotTopomapSeries()], [plotNetwork()],
#' [plotAdjacencyMatrix()], [plotSPM()].
#'
#' @section Where to go next:
#' The data model, accessors and reliability/effect-size statistics come from
#' \pkg{PhysioExperiment} (re-exported here). Filtering, referencing and ICA live
#' upstream in \pkg{PhysioPreprocess}; coupling between modalities recorded at
#' different rates lives in \pkg{PhysioCrossModal}. See the package vignettes,
#' e.g. `vignette("spectral-analysis", package = "PhysioAnalysis")`.
#'
#' @keywords internal
"_PACKAGE"
