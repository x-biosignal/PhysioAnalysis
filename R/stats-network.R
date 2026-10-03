# Network-Based Statistic (NBS; Zalesky, Fornito & Bullmore 2010): mass-
# univariate edge tests, suprathreshold connected-component identification, and
# an FWER-corrected permutation null on the largest component's size or mass.

utils::globalVariables(c("x", "y", "xend", "yend", "component", "node"))

# Coerce a list of matrices or a 3D array to an n x n x N array.
.nbs_to_array <- function(mats) {
  if (is.array(mats) && length(dim(mats)) == 3L) return(mats)
  if (is.list(mats)) {
    n <- nrow(mats[[1]])
    A <- array(0, c(n, n, length(mats)))
    for (s in seq_along(mats)) A[, , s] <- as.matrix(mats[[s]])
    return(A)
  }
  stop("Provide a list of matrices or an n x n x N array.", call. = FALSE)
}

# Pooled two-sample t per edge; grp is 1/2. X: subjects x edges.
.nbs_t_twosample <- function(X, grp) {
  g1 <- X[grp == 1L, , drop = FALSE]; g2 <- X[grp == 2L, , drop = FALSE]
  n1 <- nrow(g1); n2 <- nrow(g2)
  m1 <- colMeans(g1); m2 <- colMeans(g2)
  v1 <- .col_var(g1); v2 <- .col_var(g2)
  sp <- sqrt(((n1 - 1) * v1 + (n2 - 1) * v2) / (n1 + n2 - 2))
  (m1 - m2) / (sp * sqrt(1 / n1 + 1 / n2) + 1e-300)
}

# One-sample t per edge (paired differences D: subjects x edges).
.nbs_t_onesample <- function(D) {
  n <- nrow(D)
  m <- colMeans(D)
  v <- .col_var(D)
  m / (sqrt(v / n) + 1e-300)
}

.col_var <- function(M) {
  mu <- colMeans(M)
  colSums((M - rep(mu, each = nrow(M)))^2) / (nrow(M) - 1)
}

# Connected components over an edge list (union-find); returns edge-index groups.
.uf_components <- function(edges, n_nodes) {
  parent <- seq_len(n_nodes)
  find <- function(v) { while (parent[v] != v) v <- parent[v]; v }
  for (e in seq_len(nrow(edges))) {
    ra <- find(edges[e, 1]); rb <- find(edges[e, 2])
    if (ra != rb) parent[rb] <- ra
  }
  roots <- vapply(seq_len(nrow(edges)), function(e) find(edges[e, 1]), integer(1))
  split(seq_len(nrow(edges)), roots)
}

# Suprathreshold selection given per-edge stats and the tail.
.nbs_supra <- function(stats, thresh, tail) {
  switch(tail, right = stats > thresh, left = stats < -thresh,
         both = abs(stats) > thresh)
}

# Largest-component measure (size = edge count; mass = summed exceedance).
.nbs_max_measure <- function(stats, eidx, thresh, tail, component) {
  sel <- .nbs_supra(stats, thresh, tail)
  if (!any(sel)) return(0)
  comps <- .uf_components(eidx[sel, , drop = FALSE], max(eidx))
  sel_ix <- which(sel)
  vals <- vapply(comps, function(g) {
    if (component == "size") length(g) else sum(abs(stats[sel_ix[g]]) - thresh)
  }, numeric(1))
  max(vals)
}

#' Network-Based Statistic (NBS)
#'
#' Identifies connected subnetworks of edges that differ between two groups (or
#' two within-subject conditions) of connectivity matrices, with family-wise
#' error control from a permutation null on the largest suprathreshold component
#' (Zalesky, Fornito & Bullmore 2010). Each edge is first tested with a
#' mass-univariate t-test; edges whose statistic exceeds \code{thresh} form a
#' suprathreshold graph whose connected components are the candidate
#' subnetworks. Their size (edge count) or mass (summed statistic exceedance,
#' Smith 2009) is compared against the permutation distribution of the largest
#' component to obtain FWER-corrected component p-values.
#'
#' @param mats_group1 A list of n x n connectivity matrices or an n x n x N
#'   array (group 1, or condition 1 when \code{paired}).
#' @param mats_group2 The corresponding matrices for group 2 (or condition 2 for
#'   a paired design).
#' @param thresh Primary edge-statistic (t) threshold (default: 3).
#' @param n_perm Number of permutations (default: 1000).
#' @param tail \code{"both"}, \code{"right"} (group1 > group2), or \code{"left"}.
#' @param paired Logical; \code{TRUE} for a within-subject (paired) design, in
#'   which case the two inputs are paired condition matrices (default:
#'   \code{FALSE}).
#' @param component Component measure: \code{"size"} (edge count) or
#'   \code{"mass"} (summed statistic exceedance).
#' @param directed Logical; use all off-diagonal (directed) edges rather than the
#'   upper triangle (default: \code{FALSE}).
#' @param alpha Significance level for the returned adjacency mask (default:
#'   0.05).
#' @param seed Optional RNG seed for reproducible permutations.
#' @return A list with \code{components} (a data.frame of component
#'   \code{size}, \code{mass}, and \code{p_value}), \code{component_edges} (a
#'   list of node-pair matrices, one per component), the \code{adjacency} mask of
#'   edges in significant components, the \code{edge_stats} matrix, the
#'   \code{suprathreshold} mask, the permutation \code{null_distribution}, and
#'   the settings used.
#' @references
#' Zalesky, A., Fornito, A., & Bullmore, E. T. (2010). Network-based statistic:
#' identifying differences in brain networks. NeuroImage, 53(4), 1197-1207.
#' @seealso [plotNBSnetwork()]
#' @export
#' @examples
#' \dontrun{
#' g1 <- replicate(20, matrix(rnorm(400), 20), simplify = FALSE)
#' g2 <- replicate(20, matrix(rnorm(400), 20), simplify = FALSE)
#' nbs <- networkBasedStatistic(g1, g2, thresh = 3, n_perm = 500)
#' nbs$components
#' }
networkBasedStatistic <- function(mats_group1, mats_group2, thresh = 3,
                                  n_perm = 1000L, tail = c("both", "right", "left"),
                                  paired = FALSE, component = c("size", "mass"),
                                  directed = FALSE, alpha = 0.05, seed = NULL) {
  tail <- match.arg(tail); component <- match.arg(component)
  if (!is.null(seed)) set.seed(seed)
  A1 <- .nbs_to_array(mats_group1)
  A2 <- .nbs_to_array(mats_group2)
  n <- dim(A1)[1]
  stopifnot(dim(A1)[1] == dim(A1)[2], dim(A2)[1] == n)

  mask <- if (directed) {
    m <- matrix(TRUE, n, n); diag(m) <- FALSE; m
  } else {
    upper.tri(matrix(0, n, n))
  }
  eidx <- which(mask, arr.ind = TRUE)                 # E x 2 node pairs

  flat <- function(A) t(vapply(seq_len(dim(A)[3]), function(s) A[, , s][mask],
                               numeric(nrow(eidx))))
  if (paired) {
    stopifnot(dim(A1)[3] == dim(A2)[3])
    D <- flat(A1) - flat(A2)                          # N x E differences
    obs_stats <- .nbs_t_onesample(D)
    perm_stats <- function() .nbs_t_onesample(D * sample(c(-1, 1), nrow(D), TRUE))
  } else {
    X <- rbind(flat(A1), flat(A2))
    grp <- rep(1:2, c(dim(A1)[3], dim(A2)[3]))
    obs_stats <- .nbs_t_twosample(X, grp)
    perm_stats <- function() .nbs_t_twosample(X, sample(grp))
  }

  # observed components
  sel <- .nbs_supra(obs_stats, thresh, tail)
  comps <- if (any(sel)) .uf_components(eidx[sel, , drop = FALSE], n) else list()
  sel_ix <- which(sel)
  sizes <- vapply(comps, length, integer(1))
  masses <- vapply(comps, function(g) sum(abs(obs_stats[sel_ix[g]]) - thresh),
                   numeric(1))
  measure <- if (component == "size") sizes else masses

  # permutation null on the largest component measure
  n_perm <- as.integer(n_perm)
  null_max <- vapply(seq_len(n_perm), function(i)
    .nbs_max_measure(perm_stats(), eidx, thresh, tail, component), numeric(1))

  pvals <- if (length(comps)) {
    vapply(measure, function(m) (1 + sum(null_max >= m)) / (n_perm + 1), numeric(1))
  } else numeric(0)

  component_edges <- lapply(comps, function(g) eidx[sel_ix[g], , drop = FALSE])
  edge_stats <- matrix(0, n, n)
  edge_stats[mask] <- obs_stats
  supra_mask <- matrix(FALSE, n, n); supra_mask[mask] <- sel
  adjacency <- matrix(0L, n, n)
  if (length(pvals)) for (k in which(pvals < alpha)) {
    ee <- component_edges[[k]]
    for (r in seq_len(nrow(ee))) { adjacency[ee[r, 1], ee[r, 2]] <- 1L }
  }
  if (!directed) adjacency <- adjacency + t(adjacency)

  comp_df <- if (length(comps)) {
    data.frame(component = seq_along(comps), size = sizes, mass = masses,
               p_value = pvals)
  } else {
    data.frame(component = integer(0), size = integer(0), mass = numeric(0),
               p_value = numeric(0))
  }

  list(components = comp_df, component_edges = component_edges,
       adjacency = adjacency, edge_stats = edge_stats,
       suprathreshold = supra_mask, null_distribution = null_max,
       threshold = thresh, tail = tail, component = component,
       n_perm = n_perm, alpha = alpha, paired = paired)
}

#' Plot a significant NBS subnetwork
#'
#' Draws the edges belonging to significant NBS components on a circular node
#' layout.
#'
#' @param result A result from [networkBasedStatistic()].
#' @param alpha Significance level for the components to draw (default: uses the
#'   value stored in \code{result}).
#' @param node_labels Optional character node labels.
#' @param title Plot title.
#' @return A \code{ggplot} object.
#' @seealso [networkBasedStatistic()]
#' @export
#' @examples
#' \dontrun{
#' plotNBSnetwork(networkBasedStatistic(g1, g2))
#' }
plotNBSnetwork <- function(result, alpha = NULL, node_labels = NULL,
                           title = "NBS subnetwork") {
  stopifnot(is.list(result), !is.null(result$edge_stats))
  if (is.null(alpha)) alpha <- result$alpha
  n <- nrow(result$edge_stats)
  ang <- 2 * pi * (seq_len(n) - 1) / n
  nodes <- data.frame(node = seq_len(n), x = cos(ang), y = sin(ang))
  if (!is.null(node_labels)) nodes$label <- node_labels else
    nodes$label <- as.character(seq_len(n))

  sig <- which(result$components$p_value < alpha)
  edge_rows <- do.call(rbind, lapply(sig, function(k) {
    ee <- result$component_edges[[k]]
    data.frame(x = nodes$x[ee[, 1]], y = nodes$y[ee[, 1]],
               xend = nodes$x[ee[, 2]], yend = nodes$y[ee[, 2]],
               component = factor(k))
  }))

  p <- ggplot2::ggplot()
  if (!is.null(edge_rows) && nrow(edge_rows)) {
    p <- p + ggplot2::geom_segment(
      data = edge_rows,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend, colour = component),
      linewidth = 0.8)
  }
  p +
    ggplot2::geom_point(data = nodes, ggplot2::aes(x = x, y = y), size = 3) +
    ggplot2::geom_text(data = nodes,
                       ggplot2::aes(x = x * 1.1, y = y * 1.1, label = node),
                       size = 3) +
    ggplot2::coord_equal() +
    ggplot2::labs(title = title, colour = "component") +
    ggplot2::theme_void()
}
