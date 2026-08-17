library(testthat)
library(PhysioAnalysis)

# 4D epoched PhysioExperiment (time x channel x epoch x 1); `effect` is a
# time x channel mean added to every epoch.
make_4d_pe <- function(n_time, n_ch, n_ep, effect = NULL, sd = 1, seed = 1) {
  set.seed(seed)
  arr <- array(stats::rnorm(n_time * n_ch * n_ep, sd = sd),
               c(n_time, n_ch, n_ep, 1))
  if (!is.null(effect)) for (ep in seq_len(n_ep)) arr[, , ep, 1] <- arr[, , ep, 1] + effect
  PhysioExperiment(
    assays = list(raw = arr),
    colData = S4Vectors::DataFrame(label = paste0("Ch", seq_len(n_ch)),
                                   type = rep("EEG", n_ch)),
    samplingRate = 250)
}

test_that("TFCE matches the closed-form integral on a 1-D step function (<1e-6)", {
  E <- 0.5; H <- 2; a <- 1.7; w <- 25L
  m <- matrix(0, 60, 1); m[10:(10 + w - 1), 1] <- a
  tf <- tfce(m, E = E, H = H, tail = 1)          # dh = NULL -> exact integral
  closed <- w^E * a^(H + 1) / (H + 1)
  expect_lt(abs(tf[15, 1] - closed), 1e-6)
})

test_that("TFCE enhances both a broad-low and a narrow-high effect", {
  m <- matrix(0, 200, 1); m[20:80, 1] <- 2; m[120:125, 1] <- 6
  tf <- tfce(m, tail = 1)
  expect_gt(tf[50, 1], 0)                          # broad enhanced (via extent)
  expect_gt(tf[122, 1], 0)                         # narrow enhanced (via height)
  expect_equal(tf[150, 1], 0)                      # background stays zero
  # a single cluster-forming threshold captures only one of the two effects
  thr <- 3
  expect_false(m[50, 1] > thr)                     # broad below threshold
  expect_true(m[122, 1] > thr)                     # narrow above threshold
})

test_that("the dh grid converges to the exact integral", {
  m <- matrix(0, 100, 1); m[20:80, 1] <- 2; m[90:93, 1] <- 6
  expect_lt(max(abs(tfce(m, tail = 1) - tfce(m, dh = 0.005, tail = 1))), 0.01)
})

test_that("tail selects the sign of the enhancement", {
  m <- matrix(0, 50, 1); m[10:20, 1] <- 3; m[30:35, 1] <- -4
  expect_gt(tfce(m, tail = 1)[15, 1], 0)
  expect_equal(tfce(m, tail = 1)[32, 1], 0)
  expect_lt(tfce(m, tail = -1)[32, 1], 0)
  s <- tfce(m, tail = 0)
  expect_gt(s[15, 1], 0); expect_lt(s[32, 1], 0)
})

test_that("clusterPermutationTest(method='tfce') flags a true effect and returns a TFCE map", {
  n_time <- 40; n_ch <- 3
  effect <- matrix(0, n_time, n_ch)
  effect[8:28, ] <- 0.6                            # broad, low
  effect[33:34, ] <- 1.6                           # narrow, high
  pe <- make_4d_pe(n_time, n_ch, 20, effect = effect, seed = 1)
  res <- clusterPermutationTest(pe, method = "tfce", tail = 1,
                                n_permutations = 200, seed = 3)
  expect_equal(res$method, "tfce")
  expect_equal(dim(res$tfce), c(n_time, n_ch))
  expect_equal(dim(res$p_values), c(n_time, n_ch))
  # both the broad and the narrow regions contain significant points
  expect_true(any(res$significant[10:26, ]))
  expect_true(any(res$significant[33:34, ]))
})

test_that("TFCE permutation FWER is controlled under the global null", {
  hits <- vapply(1:12, function(sd) {
    pe <- make_4d_pe(30, 2, 18, effect = NULL, seed = sd)
    res <- clusterPermutationTest(pe, method = "tfce", tail = 1,
                                  n_permutations = 100, seed = sd)
    any(res$significant)
  }, logical(1))
  expect_lt(mean(hits), 0.2)                       # near-nominal (MC tolerance)
})

test_that("grid dh handles a map maximum smaller than dh without error", {
  m <- matrix(c(0, 0.5, 0.5, 0.5), 4, 1)
  expect_equal(tfce(m, dh = 1.0, tail = 1), matrix(0, 4, 1))  # negligible -> 0
  expect_error(tfce(m, dh = 1.0, tail = 0), NA)               # two-sided: no crash
})

test_that("asymmetric channel adjacency is symmetrized (same extent as symmetric)", {
  m <- matrix(0, 5, 3); m[2:4, 1:2] <- 3
  adj_sym <- matrix(0, 3, 3)
  adj_sym[1, 2] <- adj_sym[2, 1] <- 1; adj_sym[2, 3] <- adj_sym[3, 2] <- 1
  adj_asym <- matrix(0, 3, 3); adj_asym[1, 2] <- 1; adj_asym[2, 3] <- 1
  expect_equal(tfce(m, adjacency = adj_asym, tail = 1),
               tfce(m, adjacency = adj_sym, tail = 1))
})

test_that("tfce validates and handles empty / vector input", {
  expect_equal(tfce(matrix(0, 10, 1), tail = 1), matrix(0, 10, 1))
  v <- numeric(20); v[5:10] <- 2
  expect_true(is.matrix(tfce(v, tail = 1)))        # vector coerced to a column
})
