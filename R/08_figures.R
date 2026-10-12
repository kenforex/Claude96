# =============================================================================
# 08_figures.R
# Draws the seven figures of the paper into Figures/figureNN_<label>.pdf; the
# file number NN is the figure number in the paper:
#   Figure 1 regional_indices             Figure 5 fixed_effects_mean
#   Figure 2 break_timing_four_breaks     Figure 6 fixed_effects_spatial_weights
#   Figure 3 event_price_paths            Figure 7 break_timing (Appendix A.2)
#   Figure 4 coefficient_paths
# PDF time stamps are fixed so that reruns give identical files.
# =============================================================================

# ---- Drawing helpers (used only inside this file) ---------------------------
# open_figure()      : opens Figures/<name>.pdf with the common settings.
#                      Called by every figure function. Feeds all figures.
# close_figure()     : closes the device and fixes the PDF time stamps, so
#                      reruns give identical files. Called by every figure
#                      function. Feeds all figures.
# add_common_dates() : dotted vertical lines at the common dates. Called by
#                      figure01/04 and response_figure(). Feeds Figures 1, 4-6.
# calendar_axis()    : date axis with a tick every three years. Called by the
#                      same functions. Feeds Figures 1, 4-6.
# region_cols        : one colour per region (Figures 1-3, 7).
region_cols <- c("#375E97", "#008B8B", "#8A5B9E", "#D35D36", "#7B8B31")
open_figure <- function(name, width, height) {
  pdf(file.path("Figures", paste0(name, ".pdf")), width = width, height = height,
      family = "Helvetica", useDingbats = FALSE)
  par(family = "sans", col.axis = "#333333", col.lab = "#222222", fg = "#555555", las = 1)
}
close_figure <- function(name) {
  dev.off()
  path <- file.path("Figures", paste0(name, ".pdf"))
  bytes <- readBin(path, "raw", file.info(path)$size)
  for (key in c("/CreationDate (D:", "/ModDate (D:")) {
    for (i in grepRaw(key, bytes, fixed = TRUE, all = TRUE)) {
      bytes[i + nchar(key) + 0:13] <- charToRaw("20260101000000")
    }
  }
  writeBin(bytes, path)
}
add_common_dates <- function(dates) {
  abline(v = as.numeric(dates), col = "#B5B5B5", lty = 3, lwd = 0.7)
}
calendar_axis <- function(dates) {
  years <- as.integer(format(range(dates), "%Y"))
  axis.Date(1, at = as.Date(paste0(seq(years[1], years[2], 3), "-01-01")), format = "%Y")
}

# -----------------------------------------------------------------------------
# year_counts()
#   What it does : Number of local breaks by region and calendar year of the
#                  first month of the new regime. 5 x years matrix.
#   Called by    : figure02_four_breaks(); figure07_break_timing(); run_all.R
#   Feeds        : Figures 2, 7; Section 4.2 text (annual counts)
# -----------------------------------------------------------------------------
year_counts <- function(panel, local_breaks, years) {
  gi <- region_index(panel)
  out <- matrix(0, 5, length(years))
  for (g in 1:5) {
    y <- c()
    for (i in gi[[g]]) y <- c(y, format(panel$dates[local_breaks[[i]] + 1], "%Y"))
    out[g, ] <- as.numeric(table(factor(y, levels = years)))
  }
  out
}

# -----------------------------------------------------------------------------
# heat_panel()
#   What it does : Panel (b) of Figures 2 and 7: breaks per 100 LADs in each
#                  region and year, as a shaded grid with the numbers printed.
#   Called by    : figure02_four_breaks(); figure07_break_timing()
#   Feeds        : Figures 2, 7
# -----------------------------------------------------------------------------
heat_panel <- function(counts, n_lads, years, title) {
  par(mar = c(3.5, 8, 2.5, 1))
  z <- t(counts / n_lads * 100)
  image(years, 1:5, z, col = hcl.colors(20, "Blues 3", rev = TRUE), axes = FALSE,
      xlab = "", ylab = "", main = title)
  axis(1, at = years, labels = years, cex.axis = 0.66)
  axis(2, at = 1:5, labels = short_names, cex.axis = 0.77)
  for (g in 1:5) for (j in 1:length(years)) {
    text(years[j], g, sprintf("%.0f", z[j, g]), cex = 0.65,
        col = if (z[j, g] > max(z) * 0.55) "white" else "#222222")
  }
  box()
}

# -----------------------------------------------------------------------------
# figure01_regional()   Figure 1: regional paths of prices, sales, employment
#   Called by: run_all.R
# -----------------------------------------------------------------------------
figure01_regional <- function(panel, common_dates) {
  series <- list(panel$Y, panel$X[, , 1], panel$X[, , 2])
  titles <- c("(a) Real house prices", "(b) Sales volume", "(c) Employment measure")
  open_figure("figure01_regional_indices", 8.3, 8.5)
  par(mfrow = c(3, 1), mar = c(2.3, 4.3, 2.2, 1), oma = c(1.5, 0, 0, 0), cex = 0.87)
  for (k in 1:3) {
    m <- regional_log_means(series[[k]], panel)
    v <- 100 * exp(sweep(m, 2, m[1, ], "-"))
    matplot(panel$dates, v, type = "n", xaxt = "n", ylim = range(v), xlab = "",
        ylab = "Index: first month = 100", main = titles[k], cex.main = 1)
    add_common_dates(common_dates)
    matlines(panel$dates, v, col = region_cols, lty = 1, lwd = 1.6)
    calendar_axis(panel$dates)
    if (k == 1) legend("topleft", legend = short_names, col = region_cols, lty = 1, ncol = 3,
        bty = "n", cex = 0.72)
  }
  close_figure("figure01_regional_indices")
}

# -----------------------------------------------------------------------------
# figure02_four_breaks()   Figure 2: timing of local breaks when every LAD has
#                          the common number of breaks (four); the black line
#                          is the mean count in the common-date simulation
#   Called by: run_all.R
# -----------------------------------------------------------------------------
figure02_four_breaks <- function(panel, local4, sim_local4, common_dates) {
  years <- as.integer(format(panel$dates[1], "%Y")):as.integer(format(panel$dates[panel$T], "%Y"))
  counts <- year_counts(panel, local4, years)
  # Mean all-region count per year over the simulation draws
  sim_counts <- matrix(0, length(sim_local4), length(years))
  for (b in 1:length(sim_local4)) sim_counts[b, ] <- colSums(year_counts(panel, sim_local4[[b]], years))
  sim_mean <- colMeans(sim_counts)
  open_figure("figure02_break_timing_four_breaks", 8.3, 7.4)
  layout(matrix(1:2, 2, 1), heights = c(1, 1))
  par(mar = c(3.3, 4.8, 2.2, 1), cex = 0.88)
  top <- max(colSums(counts), sim_mean) * 1.45
  bx <- barplot(counts, names.arg = years, col = region_cols, border = NA,
      ylab = "Number of local breaks", cex.names = 0.67, ylim = c(0, top),
      main = "(a) Local breaks by first year of the new regime: four breaks per LAD")
  # Dotted vertical lines at the months of the common dates (within each year's bar)
  for (d in as.character(common_dates)) {
    j <- match(as.integer(substr(d, 1, 4)), years)
    month <- as.integer(substr(d, 6, 7))
    x <- bx[j] - 0.5 + (month - 0.5) / 12
    segments(x, 0, x, max(colSums(counts), sim_mean) * 1.08, col = "#555555", lty = 3, lwd = 1)
  }
  lines(bx, sim_mean, type = "b", pch = 19, cex = 0.55, lwd = 1.2, col = "#222222")
  legend("topright", legend = c(short_names, "All regions, common-date simulation",
      "Common dates"), fill = c(region_cols, NA, NA), border = NA,
      col = c(rep(NA, 5), "#222222", "#555555"), lty = c(rep(NA, 5), 1, 3),
      pch = c(rep(NA, 5), 19, NA), bty = "n", ncol = 3, cex = 0.66)
  heat_panel(counts, lengths(region_index(panel)), years,
      "(b) Breaks per 100 LADs in each region: four breaks per LAD")
  close_figure("figure02_break_timing_four_breaks")
}

# -----------------------------------------------------------------------------
# figure03_events()   Figure 3: regional real-price paths around four events
#   Called by: run_all.R
# -----------------------------------------------------------------------------
figure03_events <- function(panel) {
  m <- regional_log_means(panel$Y, panel)
  month_no <- function(x) as.integer(format(x, "%Y")) * 12 + as.integer(format(x, "%m"))
  events <- as.Date(c("2007-09-01", "2013-04-01", "2016-06-01", "2020-03-01"))
  titles <- c("Northern Rock support", "Help to Buy equity loan", "EU referendum",
      "First lockdown")
  open_figure("figure03_event_price_paths", 8.3, 6.8)
  layout(matrix(c(1, 2, 3, 4, 5, 5), 3, 2, byrow = TRUE), heights = c(1, 1, 0.15))
  par(mar = c(3.8, 4.5, 2.5, 1), oma = c(0, 0, 1, 0), mgp = c(2.4, 0.65, 0), cex = 0.9)
  for (e in 1:4) {
    ix <- which(abs(month_no(panel$dates) - month_no(events[e])) <= 18)
    xx <- month_no(panel$dates[ix]) - month_no(events[e])
    base <- match(-1, xx)
    v <- 100 * expm1(sweep(m[ix, ], 2, m[ix[base], ], "-"))
    matplot(xx, v, type = "n", xlab = "Months from event", ylab = "Real price change (%)",
        main = titles[e], cex.main = 0.85)
    abline(h = 0, v = 0, col = "#AAAAAA", lty = 3)
    matlines(xx, v, col = region_cols, lty = 1, lwd = 1.7)
  }
  par(mar = rep(0, 4))
  plot.new()
  legend("center", legend = short_names, col = region_cols, lty = 1, ncol = 5, bty = "n",
      cex = 0.72)
  close_figure("figure03_event_price_paths")
}

# -----------------------------------------------------------------------------
# figure04_coefficients()   Figure 4: regional median coefficient paths, local
#                           and common dates, with 95% bands (median +/- 1.96
#                           bootstrap SE)
#   Called by: run_all.R
# -----------------------------------------------------------------------------
figure04_coefficients <- function(panel, bl, bc, common_dates) {
  cols <- c("#375E97", "#777777")
  titles <- c("Spatial coefficient", "Sales slope", "Employment slope")
  open_figure("figure04_coefficient_paths", 8.3, 9.4)
  par(mfrow = c(5, 3), mar = c(2.5, 4.5, 1.8, 0.7), oma = c(1.2, 0, 0, 0),
      mgp = c(2.4, 0.65, 0), cex = 0.78)
  for (g in 1:5) for (k in 1:3) {
    v <- cbind(bl$md_point[, g, k], bc$md_point[, g, k])
    lo <- cbind(bl$md_lower[, g, k], bc$md_lower[, g, k])
    hi <- cbind(bl$md_upper[, g, k], bc$md_upper[, g, k])
    matplot(panel$dates, v, type = "n", xaxt = "n", xlab = "",
        ylab = if (k == 1) short_names[g] else "", ylim = range(v, lo, hi),
        main = if (g == 1) titles[k] else "", cex.main = 1)
    abline(h = 0, col = "#DDDDDD")
    add_common_dates(common_dates)
    for (j in c(2, 1)) polygon(c(panel$dates, rev(panel$dates)), c(lo[, j], rev(hi[, j])),
        col = adjustcolor(cols[j], alpha.f = 0.17), border = NA)
    matlines(panel$dates, v, col = cols, lty = c(1, 3), lwd = c(1.25, 1.15))
    calendar_axis(panel$dates)
  }
  close_figure("figure04_coefficient_paths")
}

# -----------------------------------------------------------------------------
# response_figure()
#   What it does : Five regions x four columns (sales direct, sales indirect,
#                  employment direct, employment indirect) of median effect
#                  paths with pale 95% bands (+/- 1.96 bootstrap SE); one or
#                  several specifications.
#   Called by    : figure05_fixed_mean(); figure06_spatial_weights()
#   Feeds        : Figures 5, 6
# -----------------------------------------------------------------------------
response_figure <- function(name, panel, specs, cols, ltys, lwds, common_dates, single) {
  titles <- c("Sales Direct\nEffect", "Sales Indirect\nEffect", "Employment Direct\nEffect",
      "Employment Indirect\nEffect")
  kk <- c(1, 1, 2, 2)                 # 1 = sales, 2 = employment
  jj <- c(1, 2, 1, 2)
  open_figure(name, 8.3, 9.4)
  par(mfrow = c(5, 4), mar = c(2.6, 5.4, 2.2, 0.6), oma = c(0.7, 0, 0, 0),
      mgp = c(3.6, 0.6, 0), cex = 0.79)
  for (g in 1:5) for (col in 1:4) {
    k <- kk[col]
    j <- jj[col]
    v <- lo <- hi <- c()
    for (s in specs) {
      v <- cbind(v, s$resp_point[, g, k, j])
      lo <- cbind(lo, s$resp_lower[, g, k, j])
      hi <- cbind(hi, s$resp_upper[, g, k, j])
    }
    ylab <- ""
    if (col == 1) ylab <- short_names[g]
    main <- ""
    if (g == 1) main <- titles[col]
    if (single) {
      plot(panel$dates, v[, 1], type = "n", xaxt = "n", xlab = "", ylab = ylab,
          ylim = range(v, lo, hi), main = main, cex.main = 0.95)
    } else {
      matplot(panel$dates, v, type = "n", xaxt = "n", xlab = "", ylab = ylab,
          ylim = range(v, lo, hi), main = main, cex.main = 0.95)
    }
    abline(h = 0, col = "#DDDDDD")
    add_common_dates(common_dates)
    if (single) {
      # Figure 5: direct (blue, solid) or indirect (orange, dashed)
      polygon(c(panel$dates, rev(panel$dates)), c(lo[, 1], rev(hi[, 1])),
          col = adjustcolor(cols[j], alpha.f = 0.17), border = NA)
      lines(panel$dates, v[, 1], col = cols[j], lty = j, lwd = 1.3)
    } else {
      # Figure 6: one line per spatial weights matrix
      for (m in ncol(v):1) polygon(c(panel$dates, rev(panel$dates)), c(lo[, m], rev(hi[, m])),
          col = adjustcolor(cols[m], alpha.f = 0.15), border = NA)
      matlines(panel$dates, v, col = cols, lty = ltys, lwd = lwds)
    }
    calendar_axis(panel$dates)
  }
  close_figure(name)
}

# -----------------------------------------------------------------------------
# figure05_fixed_mean()   Figure 5: median direct and indirect effects, with
#                         every LAD given its region's median spatial
#                         coefficient (Table 5, Panel B)
#   Called by: run_all.R
# -----------------------------------------------------------------------------
figure05_fixed_mean <- function(panel, bl, common_dates) {
  response_figure("figure05_fixed_effects_mean", panel, list(bl), c("#375E97", "#D35D36"), NULL, NULL,
      common_dates, TRUE)
}

# -----------------------------------------------------------------------------
# figure06_spatial_weights()   Figure 6: median effects under the three spatial
#                              weights matrices, each with its own regional
#                              median spatial coefficients
#   Called by: run_all.R
# -----------------------------------------------------------------------------
figure06_spatial_weights <- function(panel, bl, wa, common_dates) {
  response_figure("figure06_fixed_effects_spatial_weights", panel, list(bl, wa$contiguity, wa$knn),
      c("#375E97", "#D35D36", "#777777"), c(1, 2, 3), c(1.25, 1.25, 1.15), common_dates, FALSE)
}

# -----------------------------------------------------------------------------
# figure07_break_timing()   Figure 7 (Appendix A.2): timing of the local breaks
#                           selected by the LAD-level sequential tests
#   Called by: run_all.R
# -----------------------------------------------------------------------------
figure07_break_timing <- function(panel, local_breaks) {
  years <- as.integer(format(panel$dates[1], "%Y")):as.integer(format(panel$dates[panel$T], "%Y"))
  counts <- year_counts(panel, local_breaks, years)
  open_figure("figure07_break_timing", 8.3, 7.4)
  layout(matrix(1:2, 2, 1), heights = c(1, 1))
  par(mar = c(3.3, 4.8, 2.2, 1), cex = 0.88)
  barplot(counts, names.arg = years, col = region_cols, border = NA,
      ylab = "Number of local breaks", cex.names = 0.67,
      main = "(a) Local breaks by first year of the new regime",
      ylim = c(0, max(colSums(counts)) * 1.18))
  legend("topright", legend = short_names, fill = region_cols, bty = "n", ncol = 3, cex = 0.7)
  heat_panel(counts, lengths(region_index(panel)), years,
      "(b) Breaks per 100 LADs in each region")
  close_figure("figure07_break_timing")
}
