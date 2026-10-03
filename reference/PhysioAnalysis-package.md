# PhysioAnalysis: Analysis and Visualization for PhysioExperiment Objects

Spectral and time-frequency analysis, epoching, functional connectivity,
network metrics, statistical testing, and publication-quality plotting
for multi-modal physiological signals. Every function operates on a
`PhysioExperiment` object (from PhysioExperiment) and, where it produces
a result, returns it either as a plain value or back in the object.

## Spectral and time-frequency analysis

[`fftSignals()`](https://x-biosignal.github.io/PhysioAnalysis/reference/fftSignals.md),
[`bandPower()`](https://x-biosignal.github.io/PhysioAnalysis/reference/bandPower.md),
[`hilbertTransform()`](https://x-biosignal.github.io/PhysioAnalysis/reference/hilbertTransform.md),
[`instantaneousAmplitude()`](https://x-biosignal.github.io/PhysioAnalysis/reference/instantaneousAmplitude.md),
[`instantaneousPhase()`](https://x-biosignal.github.io/PhysioAnalysis/reference/instantaneousPhase.md),
[`spectrogram()`](https://x-biosignal.github.io/PhysioAnalysis/reference/spectrogram.md),
[`waveletTransform()`](https://x-biosignal.github.io/PhysioAnalysis/reference/waveletTransform.md),
[`specparam()`](https://x-biosignal.github.io/PhysioAnalysis/reference/specparam.md)
(aperiodic/oscillatory parameterization).

## Epoching and averaging

[`epochData()`](https://x-biosignal.github.io/PhysioAnalysis/reference/epochData.md),
[`averageEpochs()`](https://x-biosignal.github.io/PhysioAnalysis/reference/averageEpochs.md),
[`grandAverage()`](https://x-biosignal.github.io/PhysioAnalysis/reference/grandAverage.md),
[`epochTimes()`](https://x-biosignal.github.io/PhysioAnalysis/reference/epochTimes.md),
[`epochSliding()`](https://x-biosignal.github.io/PhysioAnalysis/reference/epochSliding.md).

## Connectivity and networks

Pairwise coupling with
[`coherence()`](https://x-biosignal.github.io/PhysioAnalysis/reference/coherence.md),
[`crossSpectrum()`](https://x-biosignal.github.io/PhysioAnalysis/reference/crossSpectrum.md),
[`plv()`](https://x-biosignal.github.io/PhysioAnalysis/reference/plv.md),
[`pli()`](https://x-biosignal.github.io/PhysioAnalysis/reference/pli.md),
[`wPLI()`](https://x-biosignal.github.io/PhysioAnalysis/reference/wPLI.md),
[`correlationMatrix()`](https://x-biosignal.github.io/PhysioAnalysis/reference/correlationMatrix.md)
and the general
[`connectivityMatrix()`](https://x-biosignal.github.io/PhysioAnalysis/reference/connectivityMatrix.md);
graph metrics with
[`adjacencyMatrix()`](https://x-biosignal.github.io/PhysioAnalysis/reference/adjacencyMatrix.md),
[`thresholdNetwork()`](https://x-biosignal.github.io/PhysioAnalysis/reference/thresholdNetwork.md),
[`nodeDegree()`](https://x-biosignal.github.io/PhysioAnalysis/reference/nodeDegree.md),
[`clusteringCoefficient()`](https://x-biosignal.github.io/PhysioAnalysis/reference/clusteringCoefficient.md),
[`pathLength()`](https://x-biosignal.github.io/PhysioAnalysis/reference/pathLength.md),
[`betweennessCentrality()`](https://x-biosignal.github.io/PhysioAnalysis/reference/betweennessCentrality.md),
[`globalEfficiency()`](https://x-biosignal.github.io/PhysioAnalysis/reference/globalEfficiency.md),
[`smallWorldness()`](https://x-biosignal.github.io/PhysioAnalysis/reference/smallWorldness.md),
[`modularity()`](https://x-biosignal.github.io/PhysioAnalysis/reference/modularity.md).

## Statistical testing

Epoch-wise tests
[`tTestEpochs()`](https://x-biosignal.github.io/PhysioAnalysis/reference/tTestEpochs.md),
[`anovaEpochs()`](https://x-biosignal.github.io/PhysioAnalysis/reference/anovaEpochs.md);
mass-univariate
[`clusterPermutationTest()`](https://x-biosignal.github.io/PhysioAnalysis/reference/clusterPermutationTest.md),
[`tfce()`](https://x-biosignal.github.io/PhysioAnalysis/reference/tfce.md),
[`findSignificantWindows()`](https://x-biosignal.github.io/PhysioAnalysis/reference/findSignificantWindows.md);
effect sizes
[`effectSize()`](https://x-biosignal.github.io/PhysioAnalysis/reference/effectSize.md),
[`cohensD()`](https://x-biosignal.r-universe.dev/PhysioExperiment/reference/cohensD.html),
[`rankBiserial()`](https://x-biosignal.github.io/PhysioAnalysis/reference/rankBiserial.md),
[`cliffsDelta()`](https://x-biosignal.github.io/PhysioAnalysis/reference/cliffsDelta.md);
inference helpers
[`bootstrapCI()`](https://x-biosignal.github.io/PhysioAnalysis/reference/bootstrapCI.md),
[`correctPValues()`](https://x-biosignal.github.io/PhysioAnalysis/reference/correctPValues.md).
Statistical Parametric Mapping:
[`spmTTest()`](https://x-biosignal.github.io/PhysioAnalysis/reference/spmTTest.md),
[`spmPairedTTest()`](https://x-biosignal.github.io/PhysioAnalysis/reference/spmPairedTTest.md),
[`spmAnova()`](https://x-biosignal.github.io/PhysioAnalysis/reference/spmAnova.md),
[`spmRegression()`](https://x-biosignal.github.io/PhysioAnalysis/reference/spmRegression.md),
[`spmMANOVA()`](https://x-biosignal.github.io/PhysioAnalysis/reference/spmMANOVA.md),
[`spmSnPM()`](https://x-biosignal.github.io/PhysioAnalysis/reference/spmSnPM.md).
Functional-data models
[`functionalRegression()`](https://x-biosignal.github.io/PhysioAnalysis/reference/functionalRegression.md),
[`scalarOnFunctionRegression()`](https://x-biosignal.github.io/PhysioAnalysis/reference/scalarOnFunctionRegression.md),
[`functionalMixedModel()`](https://x-biosignal.github.io/PhysioAnalysis/reference/functionalMixedModel.md).
Classical tests
([`twoSampleTTest()`](https://x-biosignal.github.io/PhysioAnalysis/reference/twoSampleTTest.md),
[`oneWayAnova()`](https://x-biosignal.github.io/PhysioAnalysis/reference/oneWayAnova.md),
[`kruskalTest()`](https://x-biosignal.github.io/PhysioAnalysis/reference/kruskalTest.md),
[`linearRegression()`](https://x-biosignal.github.io/PhysioAnalysis/reference/linearRegression.md),
...) and time-series diagnostics
([`autocorrelation()`](https://x-biosignal.github.io/PhysioAnalysis/reference/autocorrelation.md),
[`arYuleWalker()`](https://x-biosignal.github.io/PhysioAnalysis/reference/arYuleWalker.md),
[`ljungBoxTest()`](https://x-biosignal.github.io/PhysioAnalysis/reference/ljungBoxTest.md),
[`adfTest()`](https://x-biosignal.github.io/PhysioAnalysis/reference/adfTest.md),
[`kpssTest()`](https://x-biosignal.github.io/PhysioAnalysis/reference/kpssTest.md))
are also exported.

## Visualization

[`plotSignal()`](https://x-biosignal.github.io/PhysioAnalysis/reference/plotSignal.md),
[`plotMultiChannel()`](https://x-biosignal.github.io/PhysioAnalysis/reference/plotMultiChannel.md),
[`plotERP()`](https://x-biosignal.github.io/PhysioAnalysis/reference/plotERP.md),
[`plotPSD()`](https://x-biosignal.github.io/PhysioAnalysis/reference/plotPSD.md),
[`plotSpectrogram()`](https://x-biosignal.github.io/PhysioAnalysis/reference/plotSpectrogram.md),
[`plotTopomap()`](https://x-biosignal.github.io/PhysioAnalysis/reference/plotTopomap.md),
[`plotTopomapSeries()`](https://x-biosignal.github.io/PhysioAnalysis/reference/plotTopomapSeries.md),
[`plotNetwork()`](https://x-biosignal.github.io/PhysioAnalysis/reference/plotNetwork.md),
[`plotAdjacencyMatrix()`](https://x-biosignal.github.io/PhysioAnalysis/reference/plotAdjacencyMatrix.md),
[`plotSPM()`](https://x-biosignal.github.io/PhysioAnalysis/reference/plotSPM.md).

## Where to go next

The data model, accessors and reliability/effect-size statistics come
from PhysioExperiment (re-exported here). Filtering, referencing and ICA
live upstream in PhysioPreprocess; coupling between modalities recorded
at different rates lives in PhysioCrossModal. See the package vignettes,
e.g.
[`vignette("spectral-analysis", package = "PhysioAnalysis")`](https://x-biosignal.github.io/PhysioAnalysis/articles/spectral-analysis.md).

## See also

Useful links:

- <https://github.com/x-biosignal/PhysioAnalysis>

- <https://x-biosignal.r-universe.dev/PhysioAnalysis>

- Report bugs at <https://github.com/x-biosignal/PhysioAnalysis/issues>

## Author

**Maintainer**: Yusuke Matsui <mail.to.matsui@gmail.com>
