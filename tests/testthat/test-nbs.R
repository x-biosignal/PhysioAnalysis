library(testthat)
library(PhysioAnalysis)

# A group of N symmetric n x n connectivity matrices; a planted subnetwork
# (all edges among planted nodes) is elevated by `effect`.
sim_group <- function(N, n, planted = NULL, effect = 0, seed = 1) {
  set.seed(seed)
  lapply(seq_len(N), function(s) {
    M <- matrix(0, n, n); ut <- upper.tri(M)
    M[ut] <- stats::rnorm(sum(ut))
    if (!is.null(planted)) for (e in seq_len(nrow(planted)))
      M[planted[e, 1], planted[e, 2]] <- M[planted[e, 1], planted[e, 2]] + effect
    M + t(M)
  })
}

n <- 20; planted_nodes <- 1:6
planted <- t(utils::combn(planted_nodes, 2))

test_that("planted differential subnetwork is recovered as a significant component", {
  g1 <- sim_group(20, n, planted, effect = 1.0, seed = 1)
  g2 <- sim_group(20, n, NULL, effect = 0, seed = 2)
  res <- networkBasedStatistic(g1, g2, thresh = 3, n_perm = 500, tail = "right",
                               seed = 7)
  expect_s3_class(res$components, "data.frame")
  expect_gt(nrow(res$components), 0)
  expect_lt(min(res$components$p_value), 0.05)
  # the significant component's nodes lie within the planted set
  sig_k <- which.min(res$components$p_value)
  sig_nodes <- unique(as.vector(res$component_edges[[sig_k]]))
  expect_true(all(sig_nodes %in% planted_nodes))
  # adjacency mask marks significant edges only among planted nodes
  expect_true(sum(res$adjacency) > 0)
  expect_true(all(which(res$adjacency == 1L, arr.ind = TRUE) <= max(planted_nodes)))
})

test_that("no significant component under the global null (both groups null)", {
  g1 <- sim_group(20, n, NULL, 0, seed = 11)
  g2 <- sim_group(20, n, NULL, 0, seed = 22)
  res <- networkBasedStatistic(g1, g2, thresh = 3, n_perm = 500, tail = "both",
                               seed = 3)
  if (nrow(res$components)) expect_gt(min(res$components$p_value), 0.05)
  expect_equal(sum(res$adjacency), 0)
})

test_that("component p-values decrease with effect size", {
  ps <- vapply(c(0.4, 0.7, 1.0, 1.5), function(eff) {
    g1 <- sim_group(20, n, planted, effect = eff, seed = 1)
    g2 <- sim_group(20, n, NULL, 0, seed = 2)
    res <- networkBasedStatistic(g1, g2, thresh = 3, n_perm = 300,
                                 tail = "right", seed = 5)
    if (nrow(res$components)) min(res$components$p_value) else 1
  }, numeric(1))
  expect_true(all(diff(ps) <= 1e-9))            # monotone non-increasing
  expect_lt(ps[length(ps)], 0.05)               # strongest effect significant
})

test_that("false-positive rate at alpha=0.05 is controlled under the global null", {
  hits <- vapply(1:30, function(sd) {
    g1 <- sim_group(18, n, NULL, 0, seed = sd)
    g2 <- sim_group(18, n, NULL, 0, seed = sd + 100)
    res <- networkBasedStatistic(g1, g2, thresh = 3, n_perm = 200, tail = "right",
                                 seed = sd)
    if (nrow(res$components)) min(res$components$p_value) < 0.05 else FALSE
  }, logical(1))
  expect_lt(mean(hits), 0.12)                   # near-nominal control
})

test_that("paired (within-subject) NBS recovers a within-subject contrast", {
  set.seed(1)
  base <- sim_group(20, n, NULL, 0, seed = 1)   # subject baselines
  # condition 2 = baseline + planted effect (paired)
  cond1 <- base
  cond2 <- lapply(base, function(M) {
    for (e in seq_len(nrow(planted)))
      M[planted[e, 1], planted[e, 2]] <- M[planted[e, 1], planted[e, 2]] + 0.9
    M[lower.tri(M)] <- t(M)[lower.tri(M)]; M
  })
  res <- networkBasedStatistic(cond2, cond1, thresh = 3, n_perm = 500,
                               tail = "right", paired = TRUE, seed = 4)
  expect_lt(min(res$components$p_value), 0.05)
})

test_that("mass component measure and plotNBSnetwork work", {
  g1 <- sim_group(20, n, planted, effect = 1.0, seed = 1)
  g2 <- sim_group(20, n, NULL, 0, seed = 2)
  res <- networkBasedStatistic(g1, g2, thresh = 3, n_perm = 300, tail = "right",
                               component = "mass", seed = 7)
  expect_true("mass" %in% names(res$components))
  expect_lt(min(res$components$p_value), 0.05)
  skip_if_not_installed("ggplot2")
  expect_s3_class(plotNBSnetwork(res), "ggplot")
})

test_that("networkBasedStatistic validates input", {
  expect_error(networkBasedStatistic(matrix(0, 5, 5), sim_group(3, 5)), "list")
})
