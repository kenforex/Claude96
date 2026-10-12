#!/usr/bin/env Rscript
# =============================================================================
# run_all.R -- replication script for
#   "Heterogeneous Structural Breaks and Spatial Dependence in British Local
#    Housing Markets" (Pang and Chaudhuri), Ch2_Edit_FinalVersion.tex
#
# Usage, from this folder:   Rscript run_all.R
# Needs base R only (stats, graphics, parallel). About two hours on two cores.
# Input : the four files in data/.
# Output: Tables/table01_data.tex ... Tables/table13_fit.tex (13 tables) and
#         Figures/figure01_regional_indices.pdf ... figure07_break_timing.pdf
#         (7 figures). Nothing else is written to disk. Numbers that the paper
#         quotes only in the text are printed to the console, marked "[text]".
#
# Steps (the R/ file that holds the functions is in brackets):
#   1 data, seasonal adjustment and W                         [01_data.R]
#   2 date search; common dates; local dates with four breaks [02_breaks.R]
#     per LAD; break counts of the LAD-level sequential tests
#   3 CCEX-IV estimates on local and common dates             [03_iv.R]
#   4 bootstrap; regional medians; direct and indirect effects [04_bootstrap.R]
#   5 alternative spatial weights matrices                    [05_weights.R]
#   6 coverage benchmarks, penalised comparison, 12- and      [06_benchmarks.R]
#     6-month minimum regime lengths
#   7 tables and figures                                      [07_tables.R, 08_figures.R]
# Table and figure numbers below are those of the paper.
# =============================================================================

# ---- Session settings ---------------------------------------------------------
stopifnot(file.exists("data/data_set_update.csv"))    # run from the package folder
options(stringsAsFactors = FALSE, warn = 1)
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
Sys.setenv(TZ = "UTC")
invisible(Sys.setlocale("LC_TIME", "C"))             # month names in English
dir.create("Tables", showWarnings = FALSE)
dir.create("Figures", showWarnings = FALSE)
for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) source(f)

# ---- Estimation settings (Section 3) -----------------------------------------
#   min_length  minimum regime length in months (about 9% of T = 264)
#   max_breaks  zero to nine breaks in every search
#   reps, block 199 bootstrap / simulation draws, twelve-month blocks
#   seed        base seed; each random step adds its own offset
#   cores       cores for parallel::mclapply (results do not depend on it)
settings <- list(min_length = 24, max_breaks = 9, reps = 199, block = 12, seed = 20260908,
    cores = if (.Platform$OS.type == "unix") min(2, parallel::detectCores()) else 1)
start_time <- Sys.time()

# ---- Step 1: panel ------------------------------------------------------------
cat("Step 1: data\n")
panel <- read_panel("data")

# ---- Step 2: break search and the common dates ----------------------------------
cat("Step 2: date search\n")
search <- search_breaks(panel, settings)
selection <- select_common_breaks(panel, search$common_dp, search$WY, settings$min_length,
    settings$max_breaks)
common_breaks <- selection$breaks                     # last month of each old regime
common_dates <- panel$dates[common_breaks + 1]        # first month of each new regime
m_common <- length(common_breaks)                     # four common breaks
for (j in 1:(m_common + 2)) {                         # up to one step past the stop
  t <- selection$tests[j, ]
  report("Section 4.2 footnote: test %d vs %d searches %d of %d regimes, trimming %.0f-%.0f%%",
      t$null_breaks, t$null_breaks + 1, t$regimes_tested, t$regimes, 100 * t$trim_min,
      100 * t$trim_max)
}
report("Section 4.2 footnote: tabulated Bai-Perron value for 4 vs 5 breaks %.3f",
    bp_critical_values(4))
# Baseline local dates: every LAD gets the number of common breaks (four), dated
# by its own global four-break partition. The break counts of the LAD-level
# sequential tests (seq_breaks) and of the information criterion (ic_breaks) are
# checks (Table 12, Figure 7, Appendix A).
local_breaks <- breaks_with_m(search$local_dp, m_common)
seq_breaks <- breaks_with_m(search$local_dp, search$m_seq)
ic_breaks <- breaks_with_m(search$local_dp, search$m_ic)
report("Section 4.2: sequential tests select %d local breaks (%.2f per LAD); %d LADs without a break; %d at the cap",
    sum(search$m_seq), mean(search$m_seq), sum(search$m_seq == 0),
    sum(search$m_seq == settings$max_breaks))
report("Section 4.2: the information criterion selects %d breaks (%.2f per LAD), %.1f%% at the cap",
    sum(search$m_ic), mean(search$m_ic), 100 * mean(search$m_ic == settings$max_breaks))

# ---- Step 3: CCEX-IV estimates and the SSRs of Table 13 --------------------------
cat("Step 3: CCEX-IV\n")
fit_local <- fit_panel_iv(panel, local_breaks)
fit_common <- fit_panel_iv(panel, rep(list(common_breaks), panel$N))
cce_breaks <- search$cce_dp$breaks[[search$cce_m + 1]]
fit_cce <- fit_panel_cce(panel, cce_breaks)
# Five breaks in every specification (Table 13, Panel A)
fit_local5 <- fit_panel_iv(panel, breaks_with_m(search$local_dp, 5))
fit_common5 <- fit_panel_iv(panel, rep(list(search$common_dp$breaks[[6]]), panel$N))
fit_cce5 <- fit_panel_cce(panel, search$cce_dp$breaks[[6]])
local_loss <- numeric(panel$N)                        # each LAD's dating SSR
local_loss5 <- numeric(panel$N)
for (i in 1:panel$N) {
  local_loss[i] <- search$local_dp[[i]]$loss[m_common + 1]
  local_loss5[i] <- search$local_dp[[i]]$loss[6]
}
ssr_dating <- c(sum(local_loss), search$common_dp$loss[m_common + 1],
    search$cce_dp$loss[search$cce_m + 1], sum(local_loss5), search$common_dp$loss[6],
    search$cce_dp$loss[6])
ssr_final <- c(sum(fit_local$residual^2), sum(fit_common$residual^2), sum(fit_cce$residual^2),
    sum(fit_local5$residual^2), sum(fit_common5$residual^2), sum(fit_cce5$residual^2))

# ---- Step 4: bootstrap ------------------------------------------------------------
cat("Step 4: bootstrap\n")
xi <- mammen_multipliers(settings$reps, panel$T, settings$block, settings$seed + 1)
boot_local <- bootstrap_local(panel, fit_local, xi)
boot_common <- bootstrap_common(panel, fit_common, common_breaks, boot_local, xi)
report("Section 4.4: local-date median spatial coefficients in the last month by region %s",
    paste(sprintf("%.4f", boot_local$md_point[panel$T, , 1]), collapse = ", "))
report("Section 4.5 footnote: spectral radius at the regional median spatial coefficients %.4f",
    boot_local$op$spectral_radius)
# Time averages of the monthly median effects; the 95% interval is the
# average +/- 1.96 x its bootstrap SE (Section 4.5 text)
for (k in 1:2) for (j in 1:2) {
  avg <- apply(boot_local$resp_point[, , k, j], 2, mean)
  draws <- apply(boot_local$resp_draws[, , , k, j], c(1, 3), mean)
  se <- apply(draws, 2, sd)
  report("Section 4.5: %s %s effect, time averages by region %s (95%% interval excludes zero: %s)",
      c("sales", "employment")[k], c("direct", "indirect")[j],
      paste(sprintf("%.4f", avg), collapse = ", "), paste(abs(avg) > 1.96 * se, collapse = ", "))
}
boot_local$resp_draws <- NULL                         # large; not needed any more

# ---- Step 5: alternative spatial weights -------------------------------------------
cat("Step 5: alternative spatial weights\n")
weights <- weights_analysis(panel, local_breaks, xi, "data")

# ---- Step 6: benchmarks and the minimum-regime checks --------------------------------
cat("Step 6: benchmarks\n")
observed_cov <- coverage_share(local_breaks, common_breaks, 6)     # within 6 months
random <- random_placement(rep(m_common, panel$N), panel$T, settings$min_length,
    common_breaks, settings$reps, settings$seed + 2)
sim <- common_date_simulation(panel, search$WY, common_breaks, settings, search$local_cv)
penalised <- penalised_comparison(panel, search, rep(m_common, panel$N), m_common)
checks <- min_regime_checks(panel, search$WY, settings, search$local_cv, common_breaks)
# Appendix A.1: the global common partitions with 1-5 breaks under the three
# minimum regime lengths, and the shortest regime among them
same <- TRUE
shortest <- panel$T
for (m in 1:5) {
  b <- search$common_dp$breaks[[m + 1]]
  same <- same && identical(b, checks$m12$common_dp$breaks[[m + 1]]) &&
      identical(b, checks$m6$common_dp$breaks[[m + 1]])
  shortest <- min(shortest, diff(c(0, b, panel$T)))
}
report("Appendix A.1: common partitions with 1-5 breaks identical under 24, 12 and 6 months: %s; shortest regime %d months",
    same, shortest)
# Appendix A.2: the same benchmarks with the break counts of the sequential tests
observed_cov_seq <- coverage_share(seq_breaks, common_breaks, 6)
random_seq <- random_placement(search$m_seq, panel$T, settings$min_length, common_breaks,
    settings$reps, settings$seed + 2)
p_seq <- c()                                          # one-sided simulation tail probabilities
for (j in 1:m_common) p_seq <- c(p_seq,
    (1 + sum(sim$coverage[, j] <= observed_cov_seq[j])) / (settings$reps + 1))
report("Appendix A.2 (sequential counts): observed %s; random %s; simulation %s",
    paste(sprintf("%.1f", 100 * observed_cov_seq), collapse = ", "),
    paste(sprintf("%.1f", 100 * apply(random_seq, 2, mean)), collapse = ", "),
    paste(sprintf("%.1f", 100 * apply(sim$coverage, 2, mean)), collapse = ", "))
report("Appendix A.2 (sequential counts): random 97.5th percentiles %s; simulation tail probabilities %s",
    paste(sprintf("%.1f", 100 * apply(random_seq, 2, quantile, probs = 0.975, type = 8)),
    collapse = ", "), paste(sprintf("%.3f", p_seq), collapse = ", "))
report("Appendix A.1: information-criterion counts, coverage within 6 months %s",
    paste(sprintf("%.1f", 100 * coverage_share(ic_breaks, common_breaks, 6)), collapse = ", "))

# ---- Step 7: tables and figures ---------------------------------------------------------
cat("Step 7: tables and figures\n")
table01_data(panel)
table02_regional_composition(panel)
table03_prices(panel)
table04_common(panel, common_breaks, observed_cov, random, sim$coverage4, sim$common_dates)
table05_coefficients(boot_local)
table06_fullcoef(boot_local)
table07_09_common(panel, boot_common)
table10_spatial_weights(boot_local, weights)
table11_min_regime(panel, settings, selection, search$common_dp, checks)
table12_counts(panel, search, settings$max_breaks)
table13_fit(ssr_dating, ssr_final, penalised)

figure01_regional(panel, common_dates)
figure02_four_breaks(panel, local_breaks, sim$local4, common_dates)
figure03_events(panel)
figure04_coefficients(panel, boot_local, boot_common, common_dates)
figure05_fixed_mean(panel, boot_local, common_dates)
figure06_spatial_weights(panel, boot_local, weights, common_dates)
figure07_break_timing(panel, seq_breaks)
# Section 4.2 text: annual counts behind Figure 2
years <- 2004:2025
obs_year <- colSums(year_counts(panel, local_breaks, years))
sim_year <- rep(0, length(years))
for (b in 1:settings$reps) sim_year <- sim_year + colSums(year_counts(panel, sim$local4[[b]], years))
sim_year <- sim_year / settings$reps
report("Section 4.2: observed local breaks by year %s",
    paste(paste0(years, ":", obs_year), collapse = ", "))
report("Section 4.2: mean simulated local breaks by year %s",
    paste(paste0(years, ":", round(sim_year)), collapse = ", "))
cat("Done in", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 1),
    "minutes.\n")
