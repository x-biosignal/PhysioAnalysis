library(testthat)
library(PhysioAnalysis)

ncpci <- function(d, df, nfac, conf = 0.95) PhysioAnalysis:::.ncpCI(d, df, nfac, conf)

test_that("non-central-t CI satisfies its defining equation (MBESS::ci.smd algorithm)", {
  # MBESS/effectsize are not installed here; the CI limits are the ncp values
  # placing the observed t at the alpha/2 and 1-alpha/2 non-central-t quantiles.
  for (n in c(10, 20, 50)) for (d in c(0.2, 0.5, 0.8)) {
    df <- 2 * n - 2; nfac <- n * n / (2 * n)            # two-sample, n1 = n2 = n
    ci <- ncpci(d, df, nfac, 0.95)
    tval <- d * sqrt(nfac)
    expect_equal(stats::pt(tval, df, ncp = ci["upper"] * sqrt(nfac)), 0.025,
                 tolerance = 1e-6)
    expect_equal(stats::pt(tval, df, ncp = ci["lower"] * sqrt(nfac)), 0.975,
                 tolerance = 1e-6)
    expect_lt(ci["lower"], d); expect_gt(ci["upper"], d)
  }
})

test_that("the non-central-t 95% interval attains nominal coverage; normal under-covers", {
  set.seed(123)
  n1 <- 8; n2 <- 8; true_d <- 0.8; nrep <- 2000
  cn <- 0L; co <- 0L
  for (r in seq_len(nrep)) {
    x1 <- rnorm(n1, true_d, 1); x2 <- rnorm(n2, 0, 1)
    sp <- sqrt(((n1 - 1) * var(x1) + (n2 - 1) * var(x2)) / (n1 + n2 - 2))
    d <- (mean(x1) - mean(x2)) / sp
    ci <- ncpci(d, n1 + n2 - 2, n1 * n2 / (n1 + n2), 0.95)
    if (ci["lower"] <= true_d && true_d <= ci["upper"]) cn <- cn + 1L
    se <- sqrt((n1 + n2) / (n1 * n2) + d^2 / (2 * (n1 + n2)))   # old normal approx
    if (d - 1.96 * se <= true_d && true_d <= d + 1.96 * se) co <- co + 1L
  }
  nct <- cn / nrep; norm <- co / nrep
  expect_gte(nct, 0.93); expect_lte(nct, 0.97)         # attains nominal coverage
  expect_lt(norm, nct)                                 # normal approx under-covers
  # deterministic reason: the exact interval is wider for small n
  ci <- ncpci(0.8, 14, 8 * 8 / 16, 0.95)
  expect_gt(ci["upper"] - ci["lower"], 2 * 1.96 * sqrt(16 / 64 + 0.8^2 / 32))
})

test_that("Hedges bias correction scales the estimate and interval by J", {
  set.seed(1)
  epochs <- array(rnorm(3 * 1 * 24), dim = c(3, 1, 24, 1))
  pe <- PhysioExperiment(assays = list(epoched = epochs), samplingRate = 100)
  d_res <- effectSize(pe, condition1 = 1:12, condition2 = 13:24, correction = "none")
  g_res <- effectSize(pe, condition1 = 1:12, condition2 = 13:24, correction = "hedges")
  df <- 12 + 12 - 2; jj <- 1 - 3 / (4 * df - 1)
  expect_lt(jj, 1)
  expect_equal(g_res$d, d_res$d * jj, tolerance = 1e-8)
  expect_equal(g_res$ci_lower, d_res$ci_lower * jj, tolerance = 1e-8)
  expect_equal(g_res$correction, "hedges")
})

test_that("effectSize uses the non-central-t interval (not the old normal approx)", {
  set.seed(2)
  epochs <- array(rnorm(2 * 1 * 30, mean = 0.6), dim = c(2, 1, 30, 1))
  pe <- PhysioExperiment(assays = list(epoched = epochs), samplingRate = 100)
  res <- effectSize(pe, condition1 = 1:15, condition2 = 16:30)
  d <- res$d[1, 1]
  ci <- ncpci(d, 15 + 15 - 2, 15 * 15 / 30, 0.95)
  expect_equal(res$ci_lower[1, 1], unname(ci["lower"]), tolerance = 1e-6)
  expect_equal(res$ci_upper[1, 1], unname(ci["upper"]), tolerance = 1e-6)
})

test_that("rankBiserial satisfies r = 1 - 2*U/(n1*n2) and equals Cliff's delta", {
  set.seed(7)
  a <- rnorm(15, 1); b <- rnorm(18, 0)
  rb <- rankBiserial(a, b)
  u2 <- rb$n1 * rb$n2 - rb$u
  expect_equal(rb$r, 1 - 2 * u2 / (rb$n1 * rb$n2), tolerance = 1e-8)
  cd <- cliffsDelta(a, b)
  expect_equal(cd$delta, rb$r, tolerance = 1e-12)      # identical for independent samples
  expect_gte(cd$delta, -1); expect_lte(cd$delta, 1)
  expect_equal(sign(cd$delta), sign(median(a) - median(b)))
  expect_gte(cd$ci_lower, -1); expect_lte(cd$ci_upper, 1)
})

test_that("rank effect sizes span the full range and validate input", {
  expect_equal(cliffsDelta(c(4, 5, 6), c(1, 2, 3))$delta, 1)     # complete dominance
  expect_equal(cliffsDelta(c(1, 2, 3), c(4, 5, 6))$delta, -1)
  expect_equal(cliffsDelta(c(1, 2, 3), c(1, 2, 3))$delta, 0)     # identical
  expect_error(cliffsDelta(numeric(0), 1:3), "at least one")
  # ties are counted at half
  rb <- rankBiserial(c(1, 2, 3), c(2, 2, 2))
  expect_true(is.finite(rb$r) && rb$r >= -1 && rb$r <= 1)
})

# --- regression tests for adversarial-review findings (WS8-02) -----------------

test_that("Cliff's delta CI attains near-nominal coverage (consistent variance)", {
  # guards the SS-vs-var bug that collapsed coverage to ~47%
  set.seed(7)
  n1 <- 12; n2 <- 12; shift <- 0.8; nrep <- 2000
  true_delta <- 2 * pnorm(shift / sqrt(2)) - 1
  cov <- 0L
  for (r in seq_len(nrep)) {
    cd <- cliffsDelta(rnorm(n1, shift), rnorm(n2))
    if (cd$ci_lower <= true_delta && true_delta <= cd$ci_upper) cov <- cov + 1L
  }
  expect_gt(cov / nrep, 0.90)                    # near 0.95, not the broken ~0.47
  expect_lt(cov / nrep, 0.99)
})

test_that("rank effect sizes validate conf_level and .ncpCI degrades gracefully", {
  expect_error(cliffsDelta(1:5, 6:10, conf_level = 2), "conf_level")
  expect_error(rankBiserial(1:5, 6:10, conf_level = 0), "conf_level")
  # an extreme standardized difference must not abort the interval solver
  ci <- PhysioAnalysis:::.ncpCI(500, 8, 5, 0.95)
  expect_true(all(is.na(ci) | is.finite(ci)))
})

test_that("Glass's delta CI uses the control-sample degrees of freedom", {
  set.seed(1)
  epochs <- array(rnorm(2 * 2 * 30, 0.6), dim = c(2, 2, 30, 1))
  pe <- PhysioExperiment(assays = list(epoched = epochs), samplingRate = 100)
  g <- effectSize(pe, condition1 = 1:15, condition2 = 16:30, pooled = FALSE)
  d <- g$d[1, 1]
  ci_control <- PhysioAnalysis:::.ncpCI(d, 15 - 1, 15 * 15 / 30, 0.95)  # df = n2-1
  expect_equal(g$ci_lower[1, 1], unname(ci_control["lower"]), tolerance = 1e-6)
})
