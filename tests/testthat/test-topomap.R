# Tests for topographic map visualization

test_that("plotTopomap works with electrode positions", {
  pe <- PhysioExperiment(
    assays = list(raw = matrix(rnorm(400), nrow = 100, ncol = 4)),
    colData = S4Vectors::DataFrame(label = c("Fz", "Cz", "Pz", "Oz")),
    samplingRate = 100
  )

  # Apply montage to get positions

pe <- applyMontage(pe, "10-20")

  # Should work with default time point
  p <- plotTopomap(pe)
  expect_s3_class(p, "ggplot")
})

test_that("plotTopomap works with custom values", {
  pe <- PhysioExperiment(
    assays = list(raw = matrix(rnorm(400), nrow = 100, ncol = 4)),
    colData = S4Vectors::DataFrame(label = c("Fz", "Cz", "Pz", "Oz")),
    samplingRate = 100
  )
  pe <- applyMontage(pe, "10-20")

  # Should work with custom values
  p <- plotTopomap(pe, values = c(1, 0.5, -0.5, -1))
  expect_s3_class(p, "ggplot")
})

test_that("plotTopomap works with specified time point", {
  pe <- PhysioExperiment(
    assays = list(raw = matrix(rnorm(400), nrow = 100, ncol = 4)),
    colData = S4Vectors::DataFrame(label = c("Fz", "Cz", "Pz", "Oz")),
    samplingRate = 100
  )
  pe <- applyMontage(pe, "10-20")

  p <- plotTopomap(pe, time = 0.5)
  expect_s3_class(p, "ggplot")
})

test_that("plotTopomap errors without electrode positions", {
  pe <- PhysioExperiment(
    assays = list(raw = matrix(rnorm(400), nrow = 100, ncol = 4)),
    colData = S4Vectors::DataFrame(label = c("Fz", "Cz", "Pz", "Oz")),
    samplingRate = 100
  )

  expect_error(plotTopomap(pe), "positions")
})

test_that("plotTopomap respects visual options", {
  pe <- PhysioExperiment(
    assays = list(raw = matrix(rnorm(400), nrow = 100, ncol = 4)),
    colData = S4Vectors::DataFrame(label = c("Fz", "Cz", "Pz", "Oz")),
    samplingRate = 100
  )
  pe <- applyMontage(pe, "10-20")

  # Without contours
  p <- plotTopomap(pe, contours = FALSE)
  expect_s3_class(p, "ggplot")

  # Without head shape
  p <- plotTopomap(pe, head_shape = FALSE)
  expect_s3_class(p, "ggplot")

  # Without electrodes
  p <- plotTopomap(pe, electrodes = FALSE)
  expect_s3_class(p, "ggplot")
})

test_that("plotTopomapSeries creates multiple plots", {
  pe <- PhysioExperiment(
    assays = list(raw = matrix(rnorm(400), nrow = 100, ncol = 4)),
    colData = S4Vectors::DataFrame(label = c("Fz", "Cz", "Pz", "Oz")),
    samplingRate = 100
  )
  pe <- applyMontage(pe, "10-20")

  plots <- plotTopomapSeries(pe, times = c(0.1, 0.2, 0.3))

  expect_type(plots, "list")
  expect_length(plots, 3)
  expect_s3_class(plots[[1]], "ggplot")
})

test_that(".interpolateIDW is exact at electrodes and symmetric at midpoints (Shepard)", {
  idw <- PhysioAnalysis:::.interpolateIDW
  ex <- c(-1, 1); ey <- c(0, 0)

  # equal-valued electrodes: exact at a node and at the midpoint
  val <- c(5, 5)
  expect_equal(idw(ex, ey, val, -1, 0), 5, tolerance = 1e-8)  # at electrode
  expect_equal(idw(ex, ey, val,  0, 0), 5, tolerance = 1e-8)  # midpoint

  # distinct values: Shepard exactness at the nodes
  val2 <- c(2, 8)
  expect_equal(idw(ex, ey, val2, -1, 0), 2, tolerance = 1e-8)
  expect_equal(idw(ex, ey, val2,  1, 0), 8, tolerance = 1e-8)
  # equidistant midpoint under power-2 IDW -> mean of the two node values
  expect_equal(idw(ex, ey, val2,  0, 0), 5, tolerance = 1e-8)
})


# --- IDW/Shepard golden reference (WSCB-06 / WSCB-09) ---

test_that(".interpolateIDW matches the independent Shepard reference", {
  fx <- readRDS(test_path("fixtures", "idw_reference.rds"))
  got <- PhysioAnalysis:::.interpolateIDW(fx$x, fx$y, fx$values,
                                          fx$grid_x, fx$grid_y, power = fx$power)
  # matches the closed-form Shepard estimate computed by an independent route
  expect_equal(got, fx$expected, tolerance = 1e-10)
  # and is exact at the electrode positions (Shepard interpolation property)
  at_nodes <- PhysioAnalysis:::.interpolateIDW(fx$x, fx$y, fx$values,
                                               fx$x, fx$y, power = fx$power)
  expect_equal(at_nodes, fx$values, tolerance = 1e-10)
})

test_that("spherical spline is exact, constant preserving, and invariant", {
  spline <- PhysioAnalysis:::.interpolateSphericalSpline
  theta <- seq(0, 2 * pi, length.out = 9)[-9]
  x <- 0.8 * cos(theta)
  y <- 0.8 * sin(theta)
  values <- c(1, 0.5, -0.2, -0.8, -1, -0.4, 0.3, 0.9)
  query_x <- c(x, -0.3, 0, 0.25)
  query_y <- c(y, 0.2, 0, -0.4)

  got <- spline(
    x, y, values, query_x, query_y,
    center = c(0, 0), radius = 1
  )
  expect_equal(got[seq_along(values)], values, tolerance = 1e-6)

  constant <- spline(
    x, y, rep(3.25, length(x)), query_x, query_y,
    center = c(0, 0), radius = 1
  )
  expect_equal(constant, rep(3.25, length(constant)), tolerance = 1e-8)

  permutation <- c(5, 2, 8, 1, 7, 3, 6, 4)
  permuted <- spline(
    x[permutation], y[permutation], values[permutation], query_x, query_y,
    center = c(0, 0), radius = 1
  )
  expect_equal(permuted, got, tolerance = 1e-10)

  transformed <- spline(
    4 * x + 7, 4 * y - 2, values,
    4 * query_x + 7, 4 * query_y - 2,
    center = c(7, -2), radius = 4
  )
  expect_equal(transformed, got, tolerance = 1e-10)
})

test_that("spherical spline kernel and geometry contracts are enforced", {
  kernel <- PhysioAnalysis:::.sphericalSplineKernel
  cosine <- matrix(c(1, 0, 0, -1), 2, 2)
  got <- kernel(cosine)
  expect_true(all(is.finite(got)))
  expect_equal(got, t(got), tolerance = 1e-15)

  spline <- PhysioAnalysis:::.interpolateSphericalSpline
  args <- list(
    x = c(-0.7, 0.7, 0), y = c(0, 0, 0.7), values = c(1, 2, 3),
    xi = c(0, 1, 1 + 1e-4), yi = c(0, 0, 0),
    center = c(0, 0), radius = 1
  )
  boundary <- do.call(spline, args)
  expect_true(all(is.finite(boundary[1:2])))
  expect_true(is.na(boundary[3]))

  expect_error(
    spline(c(0, 0, 0.5), c(0, 0, 0.5), 1:3, 0, 0, c(0, 0), 1),
    "unique source"
  )
  expect_error(
    spline(
      c(0, 1e-12, 0.5), c(0, 0, 0.5), 1:3,
      0, 0, c(0, 0), 1
    ),
    "numerically indistinguishable"
  )
  expect_error(
    spline(c(0, 0.5), c(0, 0.5), 1:2, 0, 0, c(0, 0), 1),
    "three unique"
  )
  expect_error(
    do.call(spline, c(args, list(regularization = -1))),
    "non-negative"
  )
  expect_error(
    do.call(spline, c(args, list(stiffness = 1))),
    "at least 2"
  )
  expect_error(
    do.call(spline, c(args, list(n_terms = 0))),
    "positive integer"
  )

  regularized <- do.call(
    spline, c(args, list(regularization = 1e-5))
  )
  expect_true(all(is.finite(regularized[1:2])))

  expect_error(
    spline(1:3, 1:2, 1:3, 0, 0, c(0, 0), 3),
    "equal lengths"
  )
  expect_error(
    spline(1:3, 1:3, c(1, NA, 3), 0, 0, c(0, 0), 3),
    "finite numeric"
  )
  expect_error(
    spline(1:3, 1:3, 1:3, Inf, 0, c(0, 0), 3),
    "Query coordinates"
  )
  expect_error(
    spline(1:3, 1:3, 1:3, 0, 0, c(NA, 0), 3),
    "center"
  )
  expect_error(
    spline(1:3, 1:3, 1:3, 0, 0, c(0, 0), 0),
    "radius"
  )
  expect_error(
    spline(c(0, 0.5, 2), c(0, 0.5, 0), 1:3, 0, 0, c(0, 0), 1),
    "outside"
  )
})

test_that("spherical spline matches the pinned MNE reference", {
  fx <- readRDS(test_path("fixtures", "spherical_spline_mne_reference.rds"))
  got <- PhysioAnalysis:::.interpolateSphericalSpline(
    fx$x, fx$y, fx$values, fx$query_x, fx$query_y,
    center = fx$center,
    radius = fx$radius,
    stiffness = fx$stiffness,
    n_terms = fx$n_terms,
    regularization = fx$regularization
  )
  scale <- max(diff(range(fx$expected)), sqrt(.Machine$double.eps))
  nrmse <- sqrt(mean((got - fx$expected)^2)) / scale
  expect_lt(nrmse, 0.05)
  expect_equal(got, fx$expected, tolerance = 1e-8)
})

test_that("plotTopomap preserves IDW default and accepts spline controls", {
  set.seed(310)
  pe <- PhysioExperiment(
    assays = list(raw = matrix(rnorm(800), nrow = 100, ncol = 8)),
    colData = S4Vectors::DataFrame(
      label = c("Fp1", "Fp2", "F3", "F4", "C3", "C4", "P3", "P4")
    ),
    samplingRate = 100
  )
  pe <- applyMontage(pe, "10-20")
  values <- c(1, -1, 0.7, -0.6, 0.4, -0.3, 0.2, -0.1)

  default <- plotTopomap(pe, values = values, resolution = 21)
  explicit <- plotTopomap(
    pe, values = values, resolution = 21, interpolation = "idw"
  )
  expect_equal(default$data, explicit$data)

  positions <- getElectrodePositions(pe)
  expected <- PhysioAnalysis:::.interpolateIDW(
    positions$x, positions$y, values, default$data$x, default$data$y
  )
  expected[default$data$distance >
             max(diff(range(positions$x)), diff(range(positions$y))) * 0.6] <- NA
  expect_equal(default$data$value, expected, tolerance = 1e-10)

  spline <- plotTopomap(
    pe, values = values, resolution = 21, interpolation = "spline",
    spline_regularization = 1e-8
  )
  expect_s3_class(spline, "ggplot")
  expect_error(plotTopomap(pe, values = values, interpolation = "unknown"))

  series <- plotTopomapSeries(
    pe, times = c(0.01, 0.02), resolution = 15,
    interpolation = "spline", spline_regularization = 1e-8
  )
  expect_length(series, 2)
  expect_true(all(vapply(series, inherits, logical(1), "ggplot")))
})
