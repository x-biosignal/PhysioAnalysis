library(testthat)
library(PhysioAnalysis)

make_example <- function() {
  n <- 100
  n_ch <- 2
  n_samples <- 1
  assays <- S4Vectors::SimpleList(raw = array(rnorm(n * n_ch * n_samples), dim = c(n, n_ch, n_samples)))
  rowData <- S4Vectors::DataFrame(time_idx = seq_len(n))
  colData <- S4Vectors::DataFrame(
    label = paste0("Ch", seq_len(n_ch)),
    sensor_type = rep("EEG", n_ch)
  )
  PhysioExperiment(assays, rowData, colData, samplingRate = 500)
}

test_that("fftSignals adds an assay", {
  x <- make_example()
  y <- fftSignals(x)
  expect_true("fft" %in% SummarizedExperiment::assayNames(y))
})
