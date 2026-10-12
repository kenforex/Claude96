# =============================================================================
# 06_benchmarks.R
# Step 6. Benchmarks for the timing of local breaks and checks of the common
# dates (Sections 3 and 4.2; Appendices A.1, A.2 and B).
#   * random_placement()        : each LAD keeps its break count, dates random
#   * common_date_simulation()  : local search on pseudo-data from the fitted
#                                 common-date dating regression (199 draws),
#                                 with four breaks per LAD (baseline) and with
#                                 the sequential tests (check); also the 95%
#                                 intervals of the common dates
#   * penalised_comparison()    : BIC and LWZ comparison of the chronologies
#   * min_regime_checks()       : common and local searches with 12- and
#                                 6-month minimum regimes
# Feeds: Tables 4, 11, 13 (Panel B); Figure 2 (simulation line).
# =============================================================================

# -----------------------------------------------------------------------------
# coverage_share()
#   What it does : Share of LADs with at least one local break within +/- h
#                  months of each common date (all dates stored as the last
#                  month of the old regime).
#   Called by    : random_placement(); common_date_simulation();
#                  min_regime_checks(); run_all.R
#   Feeds        : Table 4; Appendix A.1, A.2 text
# -----------------------------------------------------------------------------
coverage_share <- function(break_list, common_breaks, h) {
  out <- numeric(length(common_breaks))
  for (j in 1:length(common_breaks)) {
    hit <- 0
    for (b in break_list) if (any(abs(b - common_breaks[j]) <= h)) hit <- hit + 1
    out[j] <- hit / length(break_list)
  }
  out
}

# -----------------------------------------------------------------------------
# draw_random_partition()
#   What it does : One partition into m + 1 regimes of at least min_length
#                  months, drawn uniformly over all admissible partitions
#                  (stars and bars on the months left over).
#   Called by    : random_placement()
#   Feeds        : Table 4 (random column); Appendix A.2 text (random
#                  placement of the sequential counts)
# -----------------------------------------------------------------------------
draw_random_partition <- function(m, months, min_length) {
  if (m == 0) return(integer(0))
  slack <- months - min_length * (m + 1)
  bars <- sort(sample.int(slack + m, m))
  extra <- diff(c(0, bars, slack + m + 1)) - 1
  as.integer(cumsum(min_length + extra)[1:m])
}

# -----------------------------------------------------------------------------
# random_placement()
#   What it does : 199 draws in which every LAD keeps its number of local
#                  breaks and the dates are placed at random; coverage of the
#                  common dates within 6 months. Matrix draws x dates.
#   Called by    : run_all.R (Step 6, with four breaks per LAD and with the
#                  sequential counts)
#   Feeds        : Table 4 (random column); Appendix A.2 text
# -----------------------------------------------------------------------------
random_placement <- function(counts, months, min_length, common_breaks, draws, seed) {
  set.seed(seed)
  out <- matrix(NA, draws, length(common_breaks))
  for (b in 1:draws) {
    placed <- list()
    for (i in 1:length(counts)) placed[[i]] <- draw_random_partition(counts[i], months, min_length)
    out[b, ] <- coverage_share(placed, common_breaks, 6)
  }
  out
}

# -----------------------------------------------------------------------------
# regime_fit()
#   What it does : OLS fitted values and residuals of each LAD's dating
#                  regression (Wy, s, e, sbar, ebar) within the regimes given
#                  by `breaks` (the same breaks for every LAD).
#   Called by    : common_date_simulation()
#   Feeds        : Table 4; Figure 2
# -----------------------------------------------------------------------------
regime_fit <- function(panel, WY, breaks) {
  bounds <- c(0, breaks, panel$T)
  fitted <- matrix(NA, panel$T, panel$N)
  residual <- matrix(NA, panel$T, panel$N)
  for (i in 1:panel$N) {
    design <- cbind(WY[, i], panel$X[, i, ], panel$F)
    for (j in 1:(length(bounds) - 1)) {
      rows <- (bounds[j] + 1):bounds[j + 1]
      q <- qr(design[rows, , drop = FALSE], tol = 1e-10)
      fitted[rows, i] <- qr.fitted(q, panel$Y[rows, i])
      residual[rows, i] <- qr.resid(q, panel$Y[rows, i])
    }
  }
  list(fitted = fitted, residual = residual)
}

# -----------------------------------------------------------------------------
# simulate_one_draw()
#   What it does : One draw of the common-date simulation: runs the local
#                  search (dynamic programming + sequential tests) for every
#                  LAD on the pseudo-data y_star, and the pooled search with
#                  the observed number of common breaks. Also keeps each LAD's
#                  four-break partition (Figure 2).
#   Called by    : common_date_simulation()
#   Feeds        : Table 4; Figure 2
# -----------------------------------------------------------------------------
simulate_one_draw <- function(y_star, panel, WY, n_common, settings, cv) {
  ml <- settings$min_length
  mb <- settings$max_breaks
  pooled <- matrix(0, panel$T, panel$T)
  local <- list()
  local4 <- list()
  for (i in 1:panel$N) {
    design <- cbind(WY[, i], panel$X[, i, ], panel$F)
    cost <- segment_costs(y_star[, i], design, ml)
    pooled <- pooled + cost
    dp <- run_dp(cost, mb, ml)
    m <- sequential_local_breaks(y_star[, i], design[, 1:3], design[, 4:5], dp, ml, mb, cv)
    local[[i]] <- dp$breaks[[m + 1]]
    local4[[i]] <- dp$breaks[[n_common + 1]]
  }
  common <- run_dp(pooled, mb, ml)$breaks[[n_common + 1]]
  list(local = local, local4 = local4, common = common)
}

# -----------------------------------------------------------------------------
# common_date_simulation()
#   What it does : Section 3 and Appendix A.2. Pseudo-data = fitted common-date dating
#                  regression + residuals x shared block multipliers (seed + 3),
#                  with all regressors (including the observed spatial lag)
#                  held fixed. 199 draws, run in parallel; each draw is
#                  independent of the number of cores because the multipliers
#                  are drawn in advance.
#   Called by    : run_all.R (Step 6)
#   Feeds        : Table 4 (simulation column with four breaks per LAD,
#                  tail probabilities, date intervals); Figure 2 (simulation
#                  line); Appendix A.2 text (sequential counts)
# -----------------------------------------------------------------------------
common_date_simulation <- function(panel, WY, common_breaks, settings, cv) {
  fit <- regime_fit(panel, WY, common_breaks)
  xi <- mammen_multipliers(settings$reps, panel$T, settings$block, settings$seed + 3)
  n_common <- length(common_breaks)
  draws <- parallel::mclapply(1:settings$reps, function(b)
    simulate_one_draw(fit$fitted + fit$residual * xi[[b]], panel, WY, n_common, settings, cv),
    mc.cores = settings$cores)
  coverage <- matrix(NA, settings$reps, n_common)       # within 6 months
  coverage4 <- matrix(NA, settings$reps, n_common)
  common_dates <- matrix(NA, settings$reps, n_common)
  total_breaks <- numeric(settings$reps)
  at_cap <- numeric(settings$reps)
  for (b in 1:settings$reps) {
    coverage[b, ] <- coverage_share(draws[[b]]$local, common_breaks, 6)
    coverage4[b, ] <- coverage_share(draws[[b]]$local4, common_breaks, 6)
    common_dates[b, ] <- draws[[b]]$common
    total_breaks[b] <- sum(lengths(draws[[b]]$local))
    at_cap[b] <- mean(lengths(draws[[b]]$local) == settings$max_breaks)
  }
  exact <- numeric(n_common)
  for (j in 1:n_common) exact[j] <- mean(common_dates[, j] == common_breaks[j])
  report("Section 4.2 and Appendix A.2: with the sequential tests the simulation finds %.0f local breaks on average, %.1f%% of LADs at the cap",
      mean(total_breaks), 100 * mean(at_cap))
  report("Appendix A.2: common dates recovered exactly in %.0f-%.0f%% of draws",
      100 * min(exact), 100 * max(exact))
  local4 <- list()
  for (b in 1:settings$reps) local4[[b]] <- draws[[b]]$local4
  list(coverage = coverage, coverage4 = coverage4, common_dates = common_dates,
      local4 = local4)
}

# -----------------------------------------------------------------------------
# penalised_comparison()
#   What it does : IC = sum_i T log(SSR_i / T) + c_T (five coefficients per
#                  regime + break dates), with c_T = log T (BIC) or
#                  0.299 (log T)^2.1 (LWZ). Compares the selected
#                  chronologies, five breaks in both, and each chronology with
#                  its break count re-selected by the criterion.
#   Called by    : run_all.R (Step 6)
#   Feeds        : Table 13, Panel B; Appendix B text
# -----------------------------------------------------------------------------
penalised_comparison <- function(panel, search, m_local, common_m) {
  nt <- panel$T
  n <- panel$N
  m_grid <- 0:(length(search$common_dp$loss) - 1)
  local_ssr <- matrix(NA, n, length(m_grid))
  common_ssr <- matrix(NA, n, length(m_grid))
  for (i in 1:n) {
    local_ssr[i, ] <- search$local_dp[[i]]$loss
    for (m in m_grid) {
      b <- c(0, search$common_dp$breaks[[m + 1]], nt)
      pieces <- c()
      for (j in 1:(length(b) - 1)) pieces <- c(pieces, search$costs[[i]][b[j] + 1, b[j + 1]])
      common_ssr[i, m + 1] <- sum(pieces)
    }
  }
  out <- data.frame()
  penalties <- c(BIC = log(nt), LWZ = 0.299 * log(nt)^2.1)
  for (k in names(penalties)) {
    c_t <- penalties[[k]]
    # LAD-level criterion for each local break count (rows LADs, columns m)
    lad_ic <- nt * log(local_ssr / nt) + c_t * outer(rep(1, n), 5 * (m_grid + 1) + m_grid)
    ic_local <- function(m) nt * sum(log(local_ssr[cbind(1:n, m + 1)] / nt)) +
      c_t * sum(5 * (m + 1) + m)
    ic_common <- function(m) nt * sum(log(common_ssr[, m + 1] / nt)) +
      c_t * (n * 5 * (m + 1) + m)
    common_path <- sapply(m_grid, ic_common)
    lad_local <- lad_ic[cbind(1:n, m_local + 1)]
    lad_common <- nt * log(common_ssr[, common_m + 1] / nt) + c_t * 5 * (common_m + 1)
    out <- rbind(out, data.frame(criterion = k,
        selected = ic_local(m_local) - ic_common(common_m),
        five_breaks = ic_local(rep(5, n)) - ic_common(5),
        reselected = sum(apply(lad_ic, 1, min)) - min(common_path),
        reselected_local_mean = mean(apply(lad_ic, 1, which.min) - 1),
        reselected_common = which.min(common_path) - 1,
        lads_favouring_local = sum(lad_local < lad_common)))
  }
  out
}

# -----------------------------------------------------------------------------
# min_regime_checks()
#   What it does : Appendix A.1. Repeats the searches with a 12- and a 6-month
#                  minimum regime (the minimum split distance of the tests
#                  changes with it):
#                  (a) common search: summed LAD losses, global partitions and
#                      the pooled sequential tests (same critical values);
#                  (b) local search: each LAD's partitions and its break count
#                      from the sequential tests (Bai-Perron values).
#   Called by    : run_all.R (Step 6)
#   Feeds        : Table 11 (common search); Appendix A.1 text (local search:
#                  sequential counts and four breaks per LAD)
# -----------------------------------------------------------------------------
min_regime_checks <- function(panel, WY, settings, cv, common_breaks) {
  out <- list()
  for (ml in c(12, 6)) {
    cat("  Minimum regime of", ml, "months\n")
    pooled <- matrix(0, panel$T, panel$T)
    dps <- list()
    for (i in 1:panel$N) {
      design <- cbind(WY[, i], panel$X[, i, ], panel$F)
      cost <- segment_costs(panel$Y[, i], design, ml)
      pooled <- pooled + cost
      dps[[i]] <- run_dp(cost, settings$max_breaks, ml)
    }
    common_dp <- run_dp(pooled, settings$max_breaks, ml)
    tests <- select_common_breaks(panel, common_dp, WY, ml, settings$max_breaks)
    # (b) local search with the shorter minimum
    m <- unlist(parallel::mclapply(1:panel$N, function(i)
      sequential_local_breaks(panel$Y[, i], cbind(WY[, i], panel$X[, i, ]), panel$F,
          dps[[i]], ml, settings$max_breaks, cv), mc.cores = settings$cores))
    cov6 <- coverage_share(breaks_with_m(dps, m), common_breaks, 6)
    cov6_four <- coverage_share(breaks_with_m(dps, length(common_breaks)), common_breaks, 6)
    report("Appendix A.1 (%d-month minimum): sequential tests select %d local breaks, %.1f%% of LADs at the cap, %d without a break",
        ml, sum(m), 100 * mean(m == settings$max_breaks), sum(m == 0))
    report("Appendix A.1 (%d-month minimum): coverage within 6 months, sequential counts %s; four breaks per LAD %s",
        ml, paste(sprintf("%.1f", 100 * cov6), collapse = ", "),
        paste(sprintf("%.1f", 100 * cov6_four), collapse = ", "))
    out[[paste0("m", ml)]] <- list(min_length = ml, common_dp = common_dp,
        selection = tests)
  }
  out
}
