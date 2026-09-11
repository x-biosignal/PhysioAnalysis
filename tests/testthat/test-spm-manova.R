test_that("single-component two-group T2 reduces to the squared t-field", {
  set.seed(7)
  n_time <- 40; grp <- rep(c("A", "B"), c(10, 12))
  mat <- matrix(rnorm(n_time * 22), nrow = n_time)
  mat[15:22, grp == "B"] <- mat[15:22, grp == "B"] + 1.3
  mv <- spmMANOVA(mat, groups = grp, vector_components = 1)
  tt <- spmTTest(mat, group1 = which(grp == "A"), group2 = which(grp == "B"))
  expect_equal(mv$T2, tt$t^2, tolerance = 1e-9)          # statistic
  expect_equal(mv$fwhm, tt$fwhm, tolerance = 1e-9)        # smoothness
  expect_equal(mv$threshold, tt$threshold^2, tolerance = 1e-6)  # RFT threshold
  expect_equal(mv$stat_name, "T2")
})

test_that("the T2 field equals the direct two-sample Hotelling T2 at a node", {
  set.seed(3)
  n_time <- 30; n1 <- 10; n2 <- 12; grp <- rep(c("A", "B"), c(n1, n2))
  arr <- array(rnorm(n_time * 22 * 3), c(n_time, 22, 3))
  arr[10:18, grp == "B", ] <- arr[10:18, grp == "B", ] + 1.0
  mv <- spmMANOVA(arr, groups = grp)
  node <- 14
  Y <- matrix(arr[node, , ], 22, 3)
  Sp <- ((n1 - 1) * cov(Y[grp == "A", ]) + (n2 - 1) * cov(Y[grp == "B", ])) /
    (22 - 2)
  d <- colMeans(Y[grp == "A", ]) - colMeans(Y[grp == "B", ])
  t2 <- (n1 * n2 / 22) * as.numeric(t(d) %*% solve(Sp) %*% d)
  expect_equal(mv$T2[node], t2, tolerance = 1e-8)
})

test_that("the >2-group X2 field equals R's Bartlett-transformed Wilks Lambda", {
  set.seed(5)
  n_time <- 30; grp <- rep(c("A", "B", "C"), c(8, 8, 9)); N <- 25
  arr <- array(rnorm(n_time * N * 2), c(n_time, N, 2))
  mv <- spmMANOVA(arr, groups = grp)
  expect_equal(mv$stat_name, "X2")
  node <- 20
  Y <- matrix(arr[node, , ], N, 2)
  lam <- summary(manova(Y ~ factor(grp)),
                 test = "Wilks")$stats["factor(grp)", "Wilks"]
  x2 <- -(N - 1 - (2 + 3) / 2) * log(lam)
  expect_equal(mv$X2[node], x2, tolerance = 1e-8)
  expect_equal(unname(mv$df["df"]), 2L * (3L - 1L))       # I*(g-1)
})

test_that("spmMANOVA accepts a list of component matrices and validates input", {
  set.seed(1)
  comps <- list(matrix(rnorm(20 * 12), 20), matrix(rnorm(20 * 12), 20))
  grp <- rep(c("A", "B"), each = 6)
  mv <- spmMANOVA(comps, groups = grp)
  expect_s3_class(mv, "spm_result")
  expect_equal(mv$n_comp, 2L)
  expect_error(spmMANOVA(comps, groups = rep("A", 12)), "at least 2 groups")
  # too few residual df for the number of components
  arr <- array(rnorm(20 * 4 * 3), c(20, 4, 3))
  expect_error(spmMANOVA(arr, groups = rep(c("A", "B"), each = 2)),
               "residual df")
})
