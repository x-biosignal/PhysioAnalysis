# Functional mixed-effects (multilevel curves), verified on synthetic data with
# a known between/within-subject variance structure.

mk_subject <- function(offset, n_stride, P, within_sd) {
  t <- seq(0, 1, length.out = P)
  t(sapply(seq_len(n_stride), function(k) offset + sin(2 * pi * t) + rnorm(P, 0, within_sd)))
}

test_that("functional ICC is high when between-subject variance dominates", {
  set.seed(1)
  P <- 51
  offs <- rnorm(10, 0, 1)                              # big between-subject spread
  Y <- do.call(rbind, lapply(offs, mk_subject, n_stride = 8, P = P, within_sd = 0.1))
  subj <- rep(seq_along(offs), each = 8)
  fm <- functionalMixedModel(Y, subj, n_perm = 0)
  expect_s3_class(fm, "functional_mixed")
  expect_gt(mean(fm$icc), 0.9)                         # subject identity stable
  expect_length(fm$icc, P)
  expect_equal(fm$n_subjects, 10)
})

test_that("functional ICC is low when within-subject noise dominates", {
  set.seed(2)
  P <- 51
  offs <- rnorm(10, 0, 0.05)                           # subjects nearly identical
  Y <- do.call(rbind, lapply(offs, mk_subject, n_stride = 8, P = P, within_sd = 1.0))
  subj <- rep(seq_along(offs), each = 8)
  fm <- functionalMixedModel(Y, subj, n_perm = 0)
  expect_lt(mean(fm$icc), 0.4)
})

test_that("subject-level group test detects a real group shift", {
  set.seed(3)
  P <- 61
  ctrl_off <- rnorm(7, 0, 0.4); pat_off <- rnorm(7, 1.2, 0.4)   # shifted group
  ctrl <- do.call(rbind, lapply(ctrl_off, mk_subject, n_stride = 10, P = P, within_sd = 0.15))
  pat  <- do.call(rbind, lapply(pat_off,  mk_subject, n_stride = 10, P = P, within_sd = 0.15))
  Y <- rbind(ctrl, pat)
  subj <- rep(1:14, each = 10); grp <- rep(c("ctrl", "pat"), each = 70)
  fm <- functionalMixedModel(Y, subj, grp, n_perm = 300, seed = 1)
  expect_lt(fm$p_global, 0.05)
  expect_length(fm$group_diff, P)
  expect_output(print(fm), "group effect")
})

test_that("subject-level test does NOT flag a null group (no pseudoreplication)", {
  set.seed(4)
  P <- 61
  # both groups from the SAME subject distribution; many strides each
  offs <- rnorm(14, 0, 0.5)
  Y <- do.call(rbind, lapply(offs, mk_subject, n_stride = 20, P = P, within_sd = 0.15))
  subj <- rep(1:14, each = 20)
  grp <- rep(c("a", "b"), each = 140)                  # arbitrary split, no real effect
  fm <- functionalMixedModel(Y, subj, grp, n_perm = 300, seed = 2)
  expect_gt(fm$p_global, 0.05)                         # correctly non-significant
})

test_that("group must be constant within subject", {
  set.seed(5)
  Y <- rbind(mk_subject(0, 4, 20, 0.2), mk_subject(1, 4, 20, 0.2))
  subj <- rep(c(1, 2), each = 4); grp <- rep(c("a", "b"), 4)   # varies within subject
  expect_error(functionalMixedModel(Y, subj, grp, n_perm = 0), "constant within")
})
