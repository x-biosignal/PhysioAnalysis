#' Topographic Map Visualization
#'
#' Functions for plotting scalp topography maps showing the spatial
#' distribution of signal values across electrode positions.

#' Plot topographic map (scalp topography)
#'
#' Creates a 2D topographic map showing the spatial distribution of values
#' across electrode positions on the scalp.
#'
#' @param x A PhysioExperiment object with electrode positions.
#' @param values Optional numeric vector of values to plot. If NULL, uses
#'   values from the specified time point.
#' @param time Time point in seconds to extract values (if values is NULL).
#' @param channel_values Named vector of channel values (alternative to values).
#' @param assay_name Optional assay name. If NULL, uses the default assay.
#' @param resolution Grid resolution for interpolation. Default is 100.
#' @param contours Logical. If TRUE, adds contour lines. Default is TRUE
#' @param head_shape Logical. If TRUE, draws head outline. Default is TRUE.
#' @param electrodes Logical. If TRUE, shows electrode positions. Default is TRUE.
#' @param palette Color palette name or vector of colors.
#' @param limits Numeric vector of length 2 for color scale limits.
#' @param title Plot title. If NULL, auto-generated.
#' @param interpolation Interpolation method. `"idw"` preserves the Shepard
#'   inverse-distance-weighted default; `"spline"` uses Perrin spherical
#'   splines on an upper-hemisphere lift.
#' @param spline_stiffness Positive integer of at least 2 controlling the
#'   spherical-spline kernel. Default is 4.
#' @param spline_terms Positive integer number of Legendre terms. Default is 50.
#' @param spline_regularization Non-negative diagonal regularization applied to
#'   the electrode kernel. The default 0 gives an interpolating spline; positive
#'   values improve conditioning but need not reproduce electrode values.
#' @return A ggplot object.
#' @details With `interpolation = "idw"`, scalp values are interpolated by
#'   inverse distance weighting (Shepard's method) with power 2. With
#'   `interpolation = "spline"`, planar montage coordinates are lifted to a
#'   shared upper unit hemisphere and evaluated with the Perrin spherical-spline
#'   kernel. This is spatial interpolation, not a surface Laplacian,
#'   current-source-density estimate, reference transformation, or source
#'   localization.
#' @references Shepard, D. (1968). "A two-dimensional interpolation function
#'   for irregularly-spaced data." \emph{Proceedings of the 1968 23rd ACM
#'   National Conference}, 517-524. \doi{10.1145/800186.810616}
#'
#'   Perrin F, Pernier J, Bertrand O, Echallier J. (1989). Spherical splines
#'   for scalp potential and current density mapping.
#'   \emph{Electroencephalography and Clinical Neurophysiology}, 72(2),
#'   184-187. \doi{10.1016/0013-4694(89)90180-6}
#' @seealso [plotTopomapSeries()] for topographic maps across time,
#'   [plotMultiChannel()] for multi-channel signal visualization,
#'   [plotERP()] for event-related potential plots.
#' @export
#' @examples
#' # Create example with 10-20 electrode positions
#' pe <- PhysioExperiment(
#'   assays = list(raw = matrix(rnorm(400), nrow = 100, ncol = 4)),
#'   colData = S4Vectors::DataFrame(label = c("Fz", "Cz", "Pz", "Oz")),
#'   samplingRate = 100
#' )
#'
#' # Apply 10-20 montage to get electrode positions
#' pe <- applyMontage(pe, "10-20")
#'
#' # Plot topographic map at time = 0.5s
#' plotTopomap(pe, time = 0.5)
#'
#' # Plot with custom values
#' plotTopomap(pe, values = c(1, 0.5, -0.5, -1))
plotTopomap <- function(x, values = NULL, time = NULL, channel_values = NULL,
                        assay_name = NULL, resolution = 100L,
                        contours = TRUE, head_shape = TRUE, electrodes = TRUE,
                        palette = "RdBu", limits = NULL, title = NULL,
                        interpolation = c("idw", "spline"),
                        spline_stiffness = 4L, spline_terms = 50L,
                        spline_regularization = 0) {
  stopifnot(inherits(x, "PhysioExperiment"))
  interpolation <- match.arg(interpolation)

  # Get electrode positions
  positions <- getElectrodePositions(x)

  if (is.null(positions)) {
    stop("Electrode positions not set. Use setElectrodePositions() or applyMontage() first.",
         call. = FALSE)
  }

  # Remove electrodes without positions
  valid_idx <- !is.na(positions$x) & !is.na(positions$y)
  if (sum(valid_idx) < 3) {
    stop("At least 3 electrodes with valid positions are required", call. = FALSE)
  }

  positions <- positions[valid_idx, ]
  n_electrodes <- nrow(positions)

  # Get values to plot
  if (!is.null(values)) {
    if (length(values) != n_electrodes) {
      stop("Length of values must match number of electrodes with positions", call. = FALSE)
    }
    plot_values <- values
  } else if (!is.null(channel_values)) {
    ch_names <- positions$channel
    if (!all(ch_names %in% names(channel_values))) {
      stop("channel_values must contain values for all electrodes with positions", call. = FALSE)
    }
    plot_values <- channel_values[ch_names]
  } else {
    # Extract values from data at specified time
    if (is.null(assay_name)) {
      assay_name <- defaultAssay(x)
    }
    data <- SummarizedExperiment::assay(x, assay_name)
    dims <- dim(data)

    # Get time index
    if (is.null(time)) {
      time_idx <- 1L
    } else {
      sr <- samplingRate(x)
      if (is.na(sr) || sr <= 0) {
        stop("Valid sampling rate required", call. = FALSE)
      }
      time_idx <- max(1L, min(dims[1], as.integer(round(time * sr)) + 1L))
    }

    # Extract values for the valid channels
    ch_indices <- which(valid_idx)
    if (length(dims) == 2) {
      plot_values <- data[time_idx, ch_indices]
    } else if (length(dims) >= 3) {
      plot_values <- data[time_idx, ch_indices, 1]
    }
  }

  if (!is.numeric(plot_values) || length(plot_values) != n_electrodes ||
      any(!is.finite(plot_values))) {
    stop("Topomap values must be finite numeric values for every electrode",
         call. = FALSE)
  }

  # Convert 3D to 2D projection (simple azimuthal equidistant)
  # Using x, y coordinates directly (assuming already in 2D head space)
  pos_x <- positions$x
  pos_y <- positions$y

  # Create interpolation grid
  x_range <- range(pos_x, na.rm = TRUE)
  y_range <- range(pos_y, na.rm = TRUE)
  margin <- 0.2 * max(diff(x_range), diff(y_range))

  grid_x <- seq(x_range[1] - margin, x_range[2] + margin, length.out = resolution)
  grid_y <- seq(y_range[1] - margin, y_range[2] + margin, length.out = resolution)
  grid <- expand.grid(x = grid_x, y = grid_y)

  # Create circular mask for head
  head_center <- c(mean(x_range), mean(y_range))
  head_radius <- max(diff(x_range), diff(y_range)) / 2 + margin * 0.5

  if (identical(interpolation, "idw")) {
    grid$value <- .interpolateIDW(
      pos_x, pos_y, plot_values, grid$x, grid$y
    )
  } else {
    source_radius <- max(sqrt(
      (pos_x - head_center[1])^2 + (pos_y - head_center[2])^2
    ))
    projection_radius <- max(head_radius, source_radius)
    grid$value <- .interpolateSphericalSpline(
      pos_x, pos_y, plot_values, grid$x, grid$y,
      center = head_center,
      radius = projection_radius,
      stiffness = spline_stiffness,
      n_terms = spline_terms,
      regularization = spline_regularization
    )
  }

  grid$distance <- sqrt((grid$x - head_center[1])^2 + (grid$y - head_center[2])^2)
  grid$value[grid$distance > head_radius] <- NA

  # Build plot
  p <- ggplot2::ggplot(grid, ggplot2::aes(x = x, y = y, fill = value))

  # Add interpolated surface
  p <- p + ggplot2::geom_raster(interpolate = TRUE)

  # Add contour lines
  if (contours) {
    p <- p + ggplot2::geom_contour(ggplot2::aes(z = value), color = "black",
                                    alpha = 0.5, bins = 10)
  }

  # Add head shape
  if (head_shape) {
    # Create head outline
    theta <- seq(0, 2 * pi, length.out = 100)
    head_outline <- data.frame(
      x = head_center[1] + head_radius * cos(theta),
      y = head_center[2] + head_radius * sin(theta)
    )

    # Nose
    nose_tip <- head_center[2] + head_radius * 1.1
    nose <- data.frame(
      x = c(head_center[1] - 0.1 * head_radius, head_center[1],
            head_center[1] + 0.1 * head_radius),
      y = c(head_center[2] + head_radius, nose_tip,
            head_center[2] + head_radius)
    )

    # Ears
    ear_left <- data.frame(
      x = c(head_center[1] - head_radius, head_center[1] - head_radius * 1.1,
            head_center[1] - head_radius),
      y = c(head_center[2] + 0.1 * head_radius, head_center[2],
            head_center[2] - 0.1 * head_radius)
    )
    ear_right <- data.frame(
      x = c(head_center[1] + head_radius, head_center[1] + head_radius * 1.1,
            head_center[1] + head_radius),
      y = c(head_center[2] + 0.1 * head_radius, head_center[2],
            head_center[2] - 0.1 * head_radius)
    )

    p <- p +
      ggplot2::geom_path(data = head_outline, ggplot2::aes(x = x, y = y),
                         inherit.aes = FALSE, linewidth = 1) +
      ggplot2::geom_path(data = nose, ggplot2::aes(x = x, y = y),
                         inherit.aes = FALSE, linewidth = 1) +
      ggplot2::geom_path(data = ear_left, ggplot2::aes(x = x, y = y),
                         inherit.aes = FALSE, linewidth = 1) +
      ggplot2::geom_path(data = ear_right, ggplot2::aes(x = x, y = y),
                         inherit.aes = FALSE, linewidth = 1)
  }

  # Add electrode positions
  if (electrodes) {
    electrode_df <- data.frame(
      x = pos_x,
      y = pos_y,
      label = positions$channel
    )

    p <- p +
      ggplot2::geom_point(data = electrode_df, ggplot2::aes(x = x, y = y),
                          inherit.aes = FALSE, size = 2, shape = 21, fill = "white") +
      ggplot2::geom_text(data = electrode_df, ggplot2::aes(x = x, y = y, label = label),
                         inherit.aes = FALSE, size = 2.5, vjust = -1)
  }

  # Set color scale
  if (!is.null(limits)) {
    p <- p + ggplot2::scale_fill_distiller(
      palette = palette, direction = -1,
      limits = limits, na.value = "transparent"
    )
  } else {
    # Symmetric limits for diverging scale
    max_abs <- max(abs(plot_values), na.rm = TRUE)
    p <- p + ggplot2::scale_fill_distiller(
      palette = palette, direction = -1,
      limits = c(-max_abs, max_abs), na.value = "transparent"
    )
  }

  # Final styling
  plot_title <- if (!is.null(title)) {
    title
  } else if (!is.null(time)) {
    sprintf("Topographic Map (t = %.3f s)", time)
  } else {
    "Topographic Map"
  }

  p <- p +
    ggplot2::coord_fixed() +
    ggplot2::labs(title = plot_title, fill = "Value") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    )

  p
}

#' Inverse distance weighting interpolation
#' @noRd
.interpolateIDW <- function(x, y, values, xi, yi, power = 2) {
  n <- length(xi)
  result <- numeric(n)

  for (i in seq_len(n)) {
    distances <- sqrt((x - xi[i])^2 + (y - yi[i])^2)

    # Handle exact matches
    if (any(distances < 1e-10)) {
      result[i] <- values[which.min(distances)]
    } else {
      weights <- 1 / (distances^power)
      result[i] <- sum(weights * values) / sum(weights)
    }
  }

  result
}

#' Perrin spherical-spline kernel
#' @noRd
.sphericalSplineKernel <- function(cosine, stiffness = 4L, n_terms = 50L) {
  if (!is.numeric(cosine) || any(!is.finite(cosine))) {
    stop("cosine values must be finite numeric values", call. = FALSE)
  }
  if (length(stiffness) != 1L || !is.finite(stiffness) ||
      stiffness != floor(stiffness) || stiffness < 2 ||
      stiffness > .Machine$integer.max) {
    stop("stiffness must be one integer of at least 2", call. = FALSE)
  }
  if (length(n_terms) != 1L || !is.finite(n_terms) ||
      n_terms != floor(n_terms) || n_terms < 1 ||
      n_terms > .Machine$integer.max) {
    stop("n_terms must be one positive integer", call. = FALSE)
  }

  stiffness <- as.integer(stiffness)
  n_terms <- as.integer(n_terms)
  cosine <- pmin(pmax(cosine, -1), 1)
  p_previous <- cosine * 0 + 1
  p_current <- cosine
  result <- 3 / (2^stiffness) * p_current

  if (n_terms >= 2L) {
    for (degree in 2:n_terms) {
      p_next <- ((2 * degree - 1) * cosine * p_current -
                   (degree - 1) * p_previous) / degree
      coefficient <- (2 * degree + 1) /
        ((degree * (degree + 1))^stiffness)
      result <- result + coefficient * p_next
      p_previous <- p_current
      p_current <- p_next
    }
  }

  result / (4 * pi)
}

#' Lift planar coordinates to a shared upper unit hemisphere
#' @noRd
.liftTopomapCoordinates <- function(x, y, center, radius, tolerance,
                                    query = FALSE) {
  if (!is.numeric(x) || !is.numeric(y) || length(x) != length(y) ||
      any(!is.finite(x)) || any(!is.finite(y))) {
    stop("Coordinate vectors must have equal lengths and finite numeric values",
         call. = FALSE)
  }
  if (!is.numeric(center) || length(center) != 2L ||
      any(!is.finite(center))) {
    stop("center must contain two finite numeric coordinates", call. = FALSE)
  }
  if (!is.numeric(radius) || length(radius) != 1L ||
      !is.finite(radius) || radius <= 0) {
    stop("radius must be one positive finite number", call. = FALSE)
  }

  u <- (x - center[1]) / radius
  v <- (y - center[2]) / radius
  radial_squared <- u^2 + v^2
  outside <- radial_squared > 1 + tolerance
  if (!query && any(outside)) {
    stop("Source coordinates lie outside the spherical projection disk",
         call. = FALSE)
  }

  inside <- !outside
  radial_squared[inside] <- pmin(1, radial_squared[inside])
  points <- matrix(NA_real_, nrow = length(x), ncol = 3L)
  points[inside, ] <- cbind(
    u[inside],
    v[inside],
    sqrt(pmax(0, 1 - radial_squared[inside]))
  )
  list(points = points, inside = inside)
}

#' Fit the augmented Perrin spherical-spline system
#' @noRd
.fitSphericalSpline <- function(source, values, stiffness, n_terms,
                                regularization, tolerance) {
  source_kernel <- .sphericalSplineKernel(
    tcrossprod(source), stiffness = stiffness, n_terms = n_terms
  )
  n_source <- nrow(source)
  system <- rbind(
    cbind(source_kernel + diag(regularization, n_source), rep(1, n_source)),
    c(rep(1, n_source), 0)
  )
  if (any(!is.finite(system))) {
    stop("Spherical-spline system is not finite", call. = FALSE)
  }
  rhs <- c(values, 0)
  decomposition <- svd(system)
  keep <- decomposition$d > tolerance * max(decomposition$d)
  if (!any(keep)) {
    stop("Spherical-spline system has no numerically stable solution",
         call. = FALSE)
  }

  projected <- crossprod(decomposition$u[, keep, drop = FALSE], rhs)
  projected <- projected / decomposition$d[keep]
  coefficients <- decomposition$v[, keep, drop = FALSE] %*% projected
  residual <- sqrt(sum((system %*% coefficients - rhs)^2)) /
    max(sqrt(sum(rhs^2)), .Machine$double.eps)
  residual_limit <- max(100 * tolerance, 1e-10)
  if (!is.finite(residual) || residual > residual_limit) {
    stop(
      sprintf(
        paste0(
          "Spherical-spline system residual %.3g exceeds %.3g; ",
          "use positive spline_regularization"
        ),
        residual, residual_limit
      ),
      call. = FALSE
    )
  }

  list(
    weights = as.numeric(coefficients[seq_len(n_source)]),
    constant = as.numeric(coefficients[n_source + 1L]),
    source_kernel = source_kernel,
    residual = residual
  )
}

#' Perrin spherical-spline interpolation
#' @noRd
.interpolateSphericalSpline <- function(
    x, y, values, xi, yi, center, radius, stiffness = 4L, n_terms = 50L,
    regularization = 0, tolerance = sqrt(.Machine$double.eps)) {
  if (!is.numeric(x) || !is.numeric(y) || length(x) != length(y) ||
      any(!is.finite(x)) || any(!is.finite(y))) {
    stop("Source coordinates must have equal lengths and be finite numeric values",
         call. = FALSE)
  }
  if (!is.numeric(values) || length(values) != length(x) ||
      any(!is.finite(values))) {
    stop("values must be finite numeric values matching source coordinates",
         call. = FALSE)
  }
  if (!is.numeric(xi) || !is.numeric(yi) || length(xi) != length(yi) ||
      any(!is.finite(xi)) || any(!is.finite(yi))) {
    stop("Query coordinates must have equal lengths and be finite",
         call. = FALSE)
  }
  if (length(regularization) != 1L || !is.numeric(regularization) ||
      !is.finite(regularization) || regularization < 0) {
    stop("regularization must be one non-negative finite number",
         call. = FALSE)
  }
  if (length(tolerance) != 1L || !is.numeric(tolerance) ||
      !is.finite(tolerance) || tolerance <= 0) {
    stop("tolerance must be one positive finite number", call. = FALSE)
  }
  if (length(x) < 3L || any(duplicated(data.frame(x = x, y = y)))) {
    stop("At least three unique source positions are required", call. = FALSE)
  }

  source_lift <- .liftTopomapCoordinates(
    x, y, center, radius, tolerance, query = FALSE
  )
  source_distances <- as.matrix(stats::dist(source_lift$points))
  diag(source_distances) <- Inf
  if (min(source_distances) <= tolerance) {
    stop("Source positions are duplicated or numerically indistinguishable",
         call. = FALSE)
  }
  query_lift <- .liftTopomapCoordinates(
    xi, yi, center, radius, tolerance, query = TRUE
  )

  ordering <- order(x, y, method = "radix")
  source <- source_lift$points[ordering, , drop = FALSE]
  values <- values[ordering]
  fit <- .fitSphericalSpline(
    source, values, stiffness, n_terms, regularization, tolerance
  )

  result <- rep(NA_real_, length(xi))
  if (any(query_lift$inside)) {
    query <- query_lift$points[query_lift$inside, , drop = FALSE]
    query_kernel <- .sphericalSplineKernel(
      query %*% t(source), stiffness = stiffness, n_terms = n_terms
    )
    result[query_lift$inside] <- as.numeric(
      query_kernel %*% fit$weights + fit$constant
    )
  }
  result
}

#' Plot topographic map animation
#'
#' Creates a series of topographic maps across time.
#'
#' @param x A PhysioExperiment object with electrode positions.
#' @param times Numeric vector of time points to plot.
#' @param ... Additional arguments passed to [plotTopomap()], including
#'   `interpolation` and the spherical-spline controls.
#' @return A list of ggplot objects.
#' @details The interpolation method and controls are forwarded unchanged to
#'   [plotTopomap()].
#' @references Shepard, D. (1968). "A two-dimensional interpolation function
#'   for irregularly-spaced data." \emph{Proceedings of the 1968 23rd ACM
#'   National Conference}, 517-524. \doi{10.1145/800186.810616}
#'
#'   Perrin F, Pernier J, Bertrand O, Echallier J. (1989). Spherical splines
#'   for scalp potential and current density mapping.
#'   \emph{Electroencephalography and Clinical Neurophysiology}, 72(2),
#'   184-187. \doi{10.1016/0013-4694(89)90180-6}
#' @seealso [plotTopomap()] for a single topographic map, [plotERP()] for
#'   event-related potential plots, [plotMultiChannel()] for multi-channel
#'   signal visualization.
#' @export
#' @examples
#' \donttest{
#' pe <- PhysioExperiment(
#'   assays = list(raw = matrix(rnorm(400), nrow = 100, ncol = 4)),
#'   colData = S4Vectors::DataFrame(label = c("Fz", "Cz", "Pz", "Oz")),
#'   samplingRate = 100
#' )
#' pe <- applyMontage(pe, "10-20")
#'
#' # Create topomaps at multiple time points
#' plots <- plotTopomapSeries(pe, times = c(0.1, 0.2, 0.3, 0.4))
#' }
plotTopomapSeries <- function(x, times, ...) {
  stopifnot(inherits(x, "PhysioExperiment"))

  plots <- lapply(times, function(t) {
    plotTopomap(x, time = t, title = sprintf("t = %.3f s", t), ...)
  })

  plots
}
