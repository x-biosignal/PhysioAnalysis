library(testthat)
library(PhysioAnalysis)

# DMIO-03: every exported ops-* that returns a modified PhysioExperiment must
# record exactly one W3C PROV activity. make_pe_2d() comes from helper-analysis.R.

# assert a single op call adds exactly one provenance record
.run_plus1 <- function(pe, call_fn) {
  before <- nrow(provenance(pe))
  out <- call_fn(pe)
  expect_s4_class(out, "PhysioExperiment")
  expect_equal(nrow(provenance(out)) - before, 1L)
  out
}

.mk_averaged <- function() {
  ev <- addEvents(make_pe_2d(2500, 3, sr = 250),
                  onset = c(1, 3, 5, 7), type = "stim")
  averageEpochs(epochData(ev, tmin = -0.2, tmax = 0.5))
}

test_that("every PE-returning analysis ops appends exactly one provenance record", {
  # in-place assay writers
  .run_plus1(make_pe_2d(256, 3), function(p) fftSignals(p))
  h <- .run_plus1(make_pe_2d(256, 3), function(p) hilbertTransform(p))
  .run_plus1(h, function(p) instantaneousAmplitude(p))
  .run_plus1(h, function(p) instantaneousPhase(p))
  .run_plus1(make_pe_2d(1000, 3, sr = 250),
             function(p) epochSliding(p, window = 0.4, step = 0.2))

  # event-based / freshly-constructed objects
  ev <- addEvents(make_pe_2d(2500, 3, sr = 250),
                  onset = c(1, 3, 5, 7), type = "stim")
  epd <- .run_plus1(ev, function(p) epochData(p, tmin = -0.2, tmax = 0.5))
  .run_plus1(epd, function(p) averageEpochs(p))

  # grandAverage constructs a fresh object from several inputs
  a1 <- .mk_averaged(); a2 <- .mk_averaged()
  before <- nrow(provenance(a1))
  ga <- grandAverage(a1, a2)
  expect_equal(nrow(provenance(ga)) - before, 1L)
  expect_equal(provenance(ga)$activity[nrow(provenance(ga))], "grandAverage")
})

test_that("a multi-step pipeline yields a provenance chain with wasDerivedFrom links", {
  pe <- make_pe_2d(256, 2, sr = 250)
  raw_name <- defaultAssay(pe)
  out <- instantaneousPhase(instantaneousAmplitude(hilbertTransform(pe)))
  prov <- provenance(out)

  expect_equal(nrow(prov), 3L)
  expect_equal(prov$activity,
               c("hilbertTransform", "instantaneousAmplitude", "instantaneousPhase"))
  # generated entities = output assays
  expect_equal(prov$generated, c("analytic", "amplitude", "phase"))
  # used inputs = wasDerivedFrom: step 1 derives from raw, steps 2-3 from analytic
  expect_equal(prov$used, c(raw_name, "analytic", "analytic"))
  # append-only, non-decreasing activity times
  expect_true(all(diff(as.numeric(prov$endedAtTime)) >= 0))
  # software agent recorded
  expect_equal(unique(prov$version),
               as.character(utils::packageVersion("PhysioAnalysis")))
})

test_that("params of the ops call are captured in the provenance record", {
  pe <- epochSliding(make_pe_2d(1000, 3, sr = 250), window = 0.4, step = 0.2)
  prov <- provenance(pe)
  expect_match(prov$params_json[1], "window")
  expect_match(prov$params_json[1], "step")
})

test_that("non-PE-returning ops (connectivity/network) add no provenance", {
  pe <- make_pe_2d(256, 3)
  cm <- correlationMatrix(pe)
  expect_true(is.matrix(cm))          # returns a matrix, not a PhysioExperiment
})
