# =============================================================================
# 05_weights.R
# Step 5. Alternative spatial weights matrices (Section 4.6).
#
# The local dates are held at the baseline estimates; only W (and therefore the
# spatial lag, the instruments and the effect decomposition) changes.
# Feeds: Table 10 (contiguity and nearest-neighbour columns) and Figure 6.
# =============================================================================

# -----------------------------------------------------------------------------
# contiguity_matrix()
#   What it does : Queen contiguity from the ONS LAD boundaries: two LADs are
#                  neighbours when their boundaries share at least one vertex.
#                  Island LADs without a land neighbour are linked to their
#                  nearest LAD. Rows are standardised.
#   Called by    : weights_analysis()
#   Feeds        : Table 10; Figure 6; Section 4.6 text (865 pairs, 4.96)
# -----------------------------------------------------------------------------
contiguity_matrix <- function(panel, data_dir) {
  x <- readLines(file.path(data_dir, "ons_lad24_boundaries_bgc.geojsonl"), warn = FALSE)
  codes <- sub('^.*"LAD24CD": *"([^"]+)".*$', "\\1", x)
  stopifnot(setequal(codes, panel$units$areacode))
  # All vertices "[lon,lat]" of each LAD, as text
  coords <- regmatches(x, gregexpr("\\[-?[0-9.]+,-?[0-9.]+\\]", x))
  owner <- rep(match(codes, panel$units$areacode), lengths(coords))
  vertex <- unlist(coords)
  keep <- !duplicated(paste(vertex, owner))
  vertex <- vertex[keep]
  owner <- owner[keep]
  # Vertices that appear in more than one LAD link those LADs
  shared <- vertex %in% vertex[duplicated(vertex)]
  owners_by_vertex <- split(owner[shared], vertex[shared])
  A <- matrix(0, panel$N, panel$N)
  for (g in owners_by_vertex) {
    if (length(g) > 1) A[as.matrix(expand.grid(g, g))] <- 1
  }
  diag(A) <- 0
  isolated <- which(rowSums(A) == 0)
  for (i in isolated) {
    j <- which.min(panel$distance[i, ])
    A[i, j] <- 1
    A[j, i] <- 1
  }
  report("Section 4.6: contiguity has %d neighbour pairs, %.2f neighbours per LAD",
      sum(A) / 2, mean(rowSums(A)))
  W <- A / rowSums(A)
  dimnames(W) <- dimnames(panel$W)
  W
}

# -----------------------------------------------------------------------------
# knn_matrix()
#   What it does : Five nearest neighbours by representative-point distance,
#                  weight 1/5 each.
#   Called by    : weights_analysis()
#   Feeds        : Table 10; Figure 6
# -----------------------------------------------------------------------------
knn_matrix <- function(panel, k = 5) {
  W <- matrix(0, panel$N, panel$N, dimnames = dimnames(panel$W))
  for (i in 1:panel$N) W[i, order(panel$distance[i, ])[1:k]] <- 1 / k
  W
}

# -----------------------------------------------------------------------------
# weights_analysis()
#   What it does : For the contiguity and nearest-neighbour matrices: CCEX-IV
#                  on the baseline local dates, 199 bootstrap draws (same
#                  multipliers as the baseline), regional medians over all
#                  LAD-months with SEs, and median direct and indirect effects
#                  with every LAD given its region's median spatial
#                  coefficient for that matrix (recomputed in every draw;
#                  bands are +/- 1.96 SE).
#   Called by    : run_all.R (Step 5)
#   Feeds        : Table 10; Figure 6; Section 4.6 text
# -----------------------------------------------------------------------------
weights_analysis <- function(panel, local_breaks, xi, data_dir) {
  gi <- summary_groups(panel)
  matrices <- list(contiguity = contiguity_matrix(panel, data_dir), knn = knn_matrix(panel))
  out <- list()
  for (w in names(matrices)) {
    pw <- panel
    pw$W <- matrices[[w]]
    fit <- fit_panel_iv(pw, local_breaks)
    point <- region_medians_all(fit$coef, gi)
    op <- fixed_operator(pw, pw$W, point[, 1])
    resp_point <- response_paths(fit$coef, pw, op)
    noise <- spatial_noise(pw, pw$W, fit)
    base <- pw$Y - noise
    draws <- array(NA, c(length(xi), 5, 3))
    resp_draws <- array(NA, c(length(xi), panel$T, 5, 2, 3))
    for (b in 1:length(xi)) {
      theta <- boot_coefficients(base + noise * xi[[b]], pw, pw$W, fit)
      draws[b, , ] <- region_medians_all(theta, gi)
      op_b <- fixed_operator(pw, pw$W, draws[b, , 1])             # recalibrated
      resp_draws[b, , , , ] <- response_paths(theta, pw, op_b)
    }
    se <- apply(draws, 2:3, sd)
    band <- normal_band(resp_point, resp_draws)
    report("Section 4.6 (%s): employment direct time averages by region %s; indirect %s", w,
        paste(sprintf("%.4f", apply(resp_point[, , 2, 1], 2, mean)), collapse = ", "),
        paste(sprintf("%.4f", apply(resp_point[, , 2, 2], 2, mean)), collapse = ", "))
    out[[w]] <- list(point = point, se = se, p = 2 * pnorm(-abs(point / se)),
        resp_point = resp_point, resp_lower = band$lower, resp_upper = band$upper)
    cat("  Alternative W done:", w, "\n")
  }
  out
}
