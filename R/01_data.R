# =============================================================================
# 01_data.R
# Step 1. Read the data, seasonally adjust prices and sales, build W.
#
# The four files in data/ are the only inputs of the whole package.
# The output of this step is the list `panel`, used by every later step:
#   panel$Y        T x N matrix of log real house prices (seasonally adjusted)
#   panel$X        T x N x 2 array: log sales (seasonally adjusted), log employment rate
#   panel$W        N x N inverse-squared-distance weights, rows sum to one
#   panel$F        T x 2 cross-sectional averages of sales and employment
#   panel$units    one row per LAD (code, name, coordinates, region)
#   panel$groups   names of the five regions; panel$dates, panel$distance, panel$N, panel$T
# Feeds: every table (1-13) and every figure (1-7).
# =============================================================================

# -----------------------------------------------------------------------------
# stl_adjust()
#   What it does : Robust periodic STL for every column of a T x N matrix of
#                  monthly series. The seasonal component minus its mean is
#                  subtracted, so each LAD keeps its sample mean (Section 4.1).
#   Called by    : read_panel()
#   Feeds        : all tables and figures (through the adjusted data)
# -----------------------------------------------------------------------------
stl_adjust <- function(z, start_year, start_month) {
  nt <- nrow(z)
  adjusted <- z
  # These are the stats::stl() defaults for a monthly series, written out.
  # With s.window = "periodic" the seasonal window is 10 * T + 1, so the
  # default s.jump is ceiling((10 * T + 1) / 10).
  s_jump <- ceiling((10 * nt + 1) / 10)
  for (i in 1:ncol(z)) {
    series <- ts(z[, i], start = c(start_year, start_month), frequency = 12)
    fit <- stl(series, s.window = "periodic", s.degree = 0, t.window = 19, t.degree = 1,
        l.window = 13, l.degree = 1, s.jump = s_jump, t.jump = 2, l.jump = 2,
        robust = TRUE, inner = 1, outer = 15)
    seasonal <- as.numeric(fit$time.series[, "seasonal"])
    adjusted[, i] <- z[, i] - (seasonal - mean(seasonal))
  }
  adjusted
}

# -----------------------------------------------------------------------------
# month_dummy_ss()
#   What it does : Sum of squares explained by eleven month dummies (after a
#                  linear trend) in a pooled regression of LAD-demeaned series.
#                  Computed before and after seasonal adjustment.
#   Called by    : read_panel()
#   Feeds        : Section 4.1, footnote on STL (98.2% and 92.7%; console only)
# -----------------------------------------------------------------------------
month_dummy_ss <- function(z, dates) {
  nt <- nrow(z)
  n <- ncol(z)
  month <- as.integer(format(dates, "%m"))
  trend <- 1:nt - mean(1:nt)
  dummies <- model.matrix(~ factor(month))[, -1]
  dummies <- sweep(dummies, 2, colMeans(dummies), "-")
  design <- cbind(trend, dummies)
  design <- design[rep(1:nt, n), ]
  y <- as.vector(sweep(z, 2, colMeans(z), "-"))
  ssr_trend <- sum(lm.fit(design[, 1, drop = FALSE], y)$residuals^2)
  ssr_full <- sum(lm.fit(design, y)$residuals^2)
  max(0, ssr_trend - ssr_full)
}

# -----------------------------------------------------------------------------
# read_panel()
#   What it does : Reads data/, orders LADs by fid, matches the ONS geography by
#                  LAD code, forms the five regions, adjusts prices and sales,
#                  and builds Y, X, W and the cross-sectional averages F.
#   Called by    : run_all.R (Step 1)
#   Feeds        : every table and figure
# -----------------------------------------------------------------------------
read_panel <- function(data_dir) {
  d <- read.csv(file.path(data_dir, "data_set_update.csv"), fileEncoding = "UTF-8-BOM",
      stringsAsFactors = FALSE, check.names = FALSE)
  geo <- read.csv(file.path(data_dir, "ons_lad24_centroids.csv"), stringsAsFactors = FALSE)
  reg <- read.csv(file.path(data_dir, "ons_lad24_regions.csv"), stringsAsFactors = FALSE)

  # Month index (12 * year + month - 1) from strings such as "2004m1"
  year <- as.integer(sub("m.*$", "", d$mdate))
  mon <- as.integer(sub("^.*m", "", d$mdate))
  d$month_id <- year * 12L + mon - 1L

  # LADs in the order of fid; months in calendar order
  units <- d[!duplicated(d$areacode), c("area", "areacode", "fid")]
  units <- units[order(units$fid), ]
  rownames(units) <- NULL
  ids <- sort(unique(d$month_id))
  n <- nrow(units)
  nt <- length(ids)
  stopifnot(nrow(d) == n * nt, all(diff(ids) == 1))
  dates <- as.Date(sprintf("%d-%02d-01", ids %/% 12, ids %% 12 + 1))
  d <- d[order(match(d$areacode, units$areacode), d$month_id), ]

  # Geography, matched by ONS code
  k <- match(units$areacode, geo$LAD24CD)
  stopifnot(!anyNA(k))
  units$easting <- geo$BNG_E[k]
  units$northing <- geo$BNG_N[k]
  units$region <- reg$RGN24NM[match(units$areacode, reg$LAD24CD)]
  units$region[substr(units$areacode, 1, 1) == "S"] <- "Scotland"
  units$region[substr(units$areacode, 1, 1) == "W"] <- "Wales"
  # Five regions (Section 4.1 and its footnote; Table 2)
  groups <- c("The North and Scotland", "Midlands and Wales", "East of England", "London",
      "The South")
  units$group <- units$region
  units$group[units$region %in% c("North East", "North West", "Yorkshire and The Humber",
      "Scotland")] <- groups[1]
  units$group[units$region %in% c("East Midlands", "West Midlands", "Wales")] <- groups[2]
  units$group[units$region %in% c("South East", "South West")] <- groups[5]
  stopifnot(all(units$group %in% groups))

  # Unadjusted series as T x N matrices (column i = LAD i)
  price_nsa <- matrix(d$realp_nsa, nt, n)
  sales_nsa <- matrix(d$lnsales_nsa, nt, n)
  employ <- matrix(d$lnemploy, nt, n)

  # Seasonal adjustment of prices and sales (employment is not adjusted)
  y0 <- as.integer(format(dates[1], "%Y"))
  m0 <- as.integer(format(dates[1], "%m"))
  price <- stl_adjust(price_nsa, y0, m0)
  sales <- stl_adjust(sales_nsa, y0, m0)
  # The data file also contains adjusted columns (realp, lnsales). They are
  # used only as a check: the file keeps seven significant digits, to which
  # robust STL reacts slightly (largest gap about 7e-4).
  gap <- max(abs(as.vector(price) - d$realp), abs(as.vector(sales) - d$lnsales))
  stopifnot(gap < 1e-3)
  ss_price <- c(month_dummy_ss(price_nsa, dates), month_dummy_ss(price, dates))
  ss_sales <- c(month_dummy_ss(sales_nsa, dates), month_dummy_ss(sales, dates))
  report("Section 4.1 (STL footnote): month-dummy contribution falls by %.1f%% (prices), %.1f%% (sales)",
      100 * (1 - ss_price[2] / ss_price[1]), 100 * (1 - ss_sales[2] / ss_sales[1]))

  X <- array(0, c(nt, n, 2))
  X[, , 1] <- sales
  X[, , 2] <- employ
  colnames(price) <- units$areacode

  # Spatial weights: inverse squared distance between representative points
  # (British National Grid metres), zero diagonal, rows sum to one (Section 3)
  distance <- as.matrix(dist(cbind(units$easting, units$northing)))
  diag(distance) <- Inf
  W <- 1 / distance^2
  W <- W / rowSums(W)
  dimnames(W) <- list(units$areacode, units$areacode)

  list(Y = price, X = X, W = W, distance = distance, dates = dates, units = units,
      groups = groups, F = cbind(rowMeans(X[, , 1]), rowMeans(X[, , 2])), N = n, T = nt)
}

# -----------------------------------------------------------------------------
# summary_groups()
#   What it does : LAD indices of each region for the coefficient and effect
#                  summaries. The City of London (E09000001) is left out: its
#                  employment rate rests on a very small survey sample
#                  (Section 3). It stays in the data, W, the cross-sectional
#                  averages, the date search and the estimation.
#   Called by    : bootstrap_local(), bootstrap_common(), response_paths() [04];
#                  weights_analysis() [05]
#   Feeds        : Tables 5-10; Figures 4-6
# -----------------------------------------------------------------------------
summary_groups <- function(panel) {
  out <- list()
  for (g in 1:5) {
    out[[g]] <- which(panel$units$group == panel$groups[g] &
        panel$units$areacode != "E09000001")
  }
  out
}

# -----------------------------------------------------------------------------
# region_index()
#   What it does : LAD indices of each region, all 349 LADs (for the data,
#                  break and coverage summaries).
#   Called by    : table and figure functions [07, 08]
#   Feeds        : Tables 2, 3, 12; Figures 1-3, 7
# -----------------------------------------------------------------------------
region_index <- function(panel) {
  out <- list()
  for (g in 1:5) out[[g]] <- which(panel$units$group == panel$groups[g])
  out
}

# -----------------------------------------------------------------------------
# report()
#   What it does : Prints a number that the paper quotes only in the text (not
#                  in a table or figure) to the console, marked "[text]".
#                  Nothing is written to disk.
#   Called by    : all steps
#   Feeds        : nothing (console only)
# -----------------------------------------------------------------------------
report <- function(...) cat("  [text]", sprintf(...), "\n")
