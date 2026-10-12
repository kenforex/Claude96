# =============================================================================
# 07_tables.R
# Writes the 13 tables of the paper to Tables/tableNN_<label>.tex, one tabular
# per file (captions and notes are in Ch2_Edit_FinalVersion.tex). The file
# number NN is the table number in the paper.
# =============================================================================

# ---- Small formatting helpers (used only inside this file) ------------------
# stars()        : significance stars from a p-value.
#                  Called by table01_data(), stacked_cell(), coef_cell().
#                  Feeds Tables 1, 5-10.
# stacked_cell() : estimate with stars above its SE in parentheses.
#                  Called by table05/06/07_09/10 functions. Feeds Tables 5, 6, 8, 9, 10.
# coef_cell()    : estimate and SE through the \coef macro of the manuscript.
#                  Called by table07_09_common(). Feeds Table 7.
# two_line()     : two lines inside one cell.
#                  Called by table04_common(). Feeds Table 4.
# tex_row()      : one table row from a vector of cells.
#                  Called by every table function. Feeds all tables.
# write_table()  : writes Tables/<name>.tex with booktabs rules.
#                  Called by every table function. Feeds all tables.
# fmt_p()        : p-value printed as "<0.001" when small.
#                  Called by table11_min_regime(). Feeds Table 11.
stars <- function(p) {
  if (p < 0.01) return("***")
  if (p < 0.05) return("**")
  if (p < 0.10) return("*")
  ""
}
stacked_cell <- function(est, se, p) {
  sprintf("\\begin{tabular}[t]{@{}c@{}}$%.4f^{%s}$\\\\[-1pt]$(%.4f)$\\end{tabular}",
      est, stars(p), se)
}
coef_cell <- function(est, se, p) {
  s <- stars(p)
  if (s == "") return(sprintf("\\coef{%.4f}{%.4f}", est, se))
  sprintf("\\coef{%.4f^{%s}}{%.4f}", est, s, se)
}
two_line <- function(top, bottom) {
  paste0("\\begin{tabular}[t]{@{}c@{}}", top, "\\\\[-1pt]", bottom, "\\end{tabular}")
}
tex_row <- function(cells, gap = "") paste0(paste(cells, collapse = " & "), "\\\\", gap)
write_table <- function(name, open, header, body, close) {
  lines <- c(open, "\\toprule", header, "\\midrule", body, "\\bottomrule", close)
  writeLines(lines, file.path("Tables", paste0(name, ".tex")), useBytes = TRUE)
}
fmt_p <- function(x) {
  if (!is.finite(x)) return("--")
  if (x < 0.001) return("$<$0.001")
  sprintf("%.3f", x)
}
short_names <- c("North / Scotland", "Midlands / Wales", "East", "London", "South")

# -----------------------------------------------------------------------------
# cadf_statistics()
#   What it does : Cross-sectionally augmented Dickey-Fuller t-statistics
#                  (intercept, trend, own and average lagged differences), one
#                  per LAD.
#   Called by    : table01_data()
#   Feeds        : Table 1 (CIPS)
# -----------------------------------------------------------------------------
cadf_statistics <- function(Y, lags) {
  nt <- nrow(Y)
  DY <- diff(Y)
  barY <- rowMeans(Y)
  barDY <- diff(barY)
  ix <- (lags + 2):nt
  out <- numeric(ncol(Y))
  for (i in 1:ncol(Y)) {
    own_lags <- sapply(1:lags, function(l) DY[ix - 1 - l, i])
    avg_lags <- sapply(1:lags, function(l) barDY[ix - 1 - l])
    design <- cbind(1, Y[ix - 1, i], own_lags, barDY[ix - 1], barY[ix - 1], avg_lags, ix)
    fit <- lm.fit(design, DY[ix - 1, i])
    inv <- chol2inv(qr.R(fit$qr))
    se <- sqrt(sum(fit$residuals^2) / (length(ix) - fit$rank) * diag(inv))
    out[i] <- fit$coefficients[2] / se[which(fit$qr$pivot == 2)]
  }
  out
}

# -----------------------------------------------------------------------------
# table01_data()   Table 1: descriptive statistics, Pesaran CD and CIPS
#   Called by: run_all.R
# -----------------------------------------------------------------------------
table01_data <- function(panel) {
  series <- list(panel$Y, panel$X[, , 1], panel$X[, , 2])
  labels <- c("\\textit{realp} (deseasonalised)", "\\textit{lnsales} (deseasonalised)",
      "\\textit{lnemploy}")
  rows <- c()
  for (k in 1:3) {
    z <- series[[k]]
    # Pesaran CD from the pairwise correlations of the log levels
    R <- cor(z)
    n <- ncol(z)
    cd <- sqrt(2 * nrow(z) / (n * (n - 1))) * sum(R[upper.tri(R)])
    cd_p <- 2 * pnorm(-abs(cd))
    # CIPS*: mean of CADF statistics (one lag) truncated to [-6.42, 1.70]
    cips <- mean(pmax(-6.42, pmin(1.70, cadf_statistics(z, 1))))
    cips_st <- ""
    if (cips < -2.51) cips_st <- "*"
    if (cips < -2.55) cips_st <- "**"
    if (cips < -2.62) cips_st <- "***"
    rows <- c(rows, tex_row(c(labels[k], sprintf("%.3f", c(mean(z), sd(as.vector(z)), min(z),
        max(z))), sprintf("$%.3f^{%s}$", cd, stars(cd_p)), sprintf("$%.3f^{%s}$", cips,
        cips_st))))
  }
  write_table("table01_data", "\\begin{tabular}{@{}lrrrrrr@{}}",
      "Variable & Mean & Std & Min & Max & CD & CIPS\\\\", rows, "\\end{tabular}")
}

# -----------------------------------------------------------------------------
# table02_regional_composition()   Table 2: regional groups and coverage
#   Called by: run_all.R
# -----------------------------------------------------------------------------
table02_regional_composition <- function(panel) {
  parts <- list(c("North East", "North West", "Yorkshire and The Humber", "Scotland"),
      c("East Midlands", "West Midlands", "Wales"), "East of England", "London",
      c("South East", "South West"))
  names5 <- c("North and Scotland", "Midlands and Wales", "East", "London", "South")
  gi <- region_index(panel)
  rows <- c()
  for (g in 1:5) {
    text <- c()
    for (r in parts[[g]]) {
      count <- sum(panel$units$region == r)
      text <- c(text, paste0(sub("The Humber", "the Humber", r), " (", count, ")"))
    }
    # "A, B, and C" / "A and B" / "A"
    if (length(text) == 2) text <- paste(text, collapse = " and ")
    if (length(text) > 2) text <- paste0(paste(text[-length(text)], collapse = ", "),
        ", and ", text[length(text)])
    n <- length(gi[[g]])
    rows <- c(rows, tex_row(c(names5[g], text, n, sprintf("%.1f", 100 * n / panel$N),
        formatC(n * panel$T, format = "d", big.mark = ",")), if (g <= 4) "[3pt]" else ""))
  }
  total <- tex_row(c("All LADs", "England, Scotland and Wales", panel$N, sprintf("%.1f", 100),
      formatC(panel$N * panel$T, format = "d", big.mark = ",")))
  body <- c(rows[1], "", rows[2], "", rows[3], "", rows[4], "", rows[5], "\\midrule", total)
  write_table("table02_regional_composition",
      c("\\begin{tabularx}{\\linewidth}{@{}>{\\raggedright\\arraybackslash}p{3.2cm}",
      "    >{\\raggedright\\arraybackslash}X rrr@{}}"),
      "Region & Constituent areas (sample LADs) & LADs & Share (\\%) & Observations\\\\",
      body, "\\end{tabularx}")
}

# -----------------------------------------------------------------------------
# regional_log_means()
#   What it does : Monthly equally weighted regional mean of a T x N series
#                  (all LADs). T x 5 matrix.
#   Called by    : table03_prices(); figure01_regional(); figure03_events() [08]
#   Feeds        : Table 3; Figures 1, 3
# -----------------------------------------------------------------------------
regional_log_means <- function(z, panel) {
  gi <- region_index(panel)
  out <- matrix(NA, panel$T, 5)
  for (g in 1:5) out[, g] <- rowMeans(z[, gi[[g]], drop = FALSE])
  out
}

# -----------------------------------------------------------------------------
# table03_prices()   Table 3: regional real-price changes over five episodes
#   Called by: run_all.R
# -----------------------------------------------------------------------------
table03_prices <- function(panel) {
  starts <- as.Date(c("2007-06-01", "2012-12-01", "2016-06-01", "2020-02-01", "2022-02-01"))
  ends <- as.Date(c("2009-06-01", "2015-12-01", "2019-06-01", "2022-02-01", "2025-02-01"))
  labels <- c("2007-09", "2012-15", "2016-19", "2020-22", "2022-25")
  m <- regional_log_means(panel$Y, panel)
  rows <- c()
  for (g in 1:5) {
    change <- 100 * expm1(m[match(ends, panel$dates), g] - m[match(starts, panel$dates), g])
    rows <- c(rows, tex_row(c(short_names[g], sprintf("%+.1f", change))))
  }
  write_table("table03_prices", "\\begin{tabular}{lrrrrr}",
      paste0("Region & ", paste(labels, collapse = " & "), "\\\\"), rows, "\\end{tabular}")
}

# -----------------------------------------------------------------------------
# table04_common()   Table 4: common dates, their 95% intervals from the
#                    common-date simulation, and the
#                    coverage of the baseline local breaks (four per LAD)
#                    within six months, with the two benchmarks, also with
#                    four breaks per LAD (mean and 2.5th-97.5th percentiles
#                    over the 199 draws).
#   Called by: run_all.R
# -----------------------------------------------------------------------------
table04_common <- function(panel, common_breaks, observed, random, sim_coverage, sim_dates) {
  pct <- function(x) sprintf("%.1f", 100 * x)
  reps <- dim(sim_coverage)[1]
  rows <- c()
  for (j in 1:length(common_breaks)) {
    ci <- quantile(sim_dates[, j], c(0.025, 0.975), type = 1)
    cells <- c(format(panel$dates[common_breaks[j] + 1], "%b %Y"),
        two_line(paste0(format(panel$dates[ci[1] + 1], "%b %Y"), "--"),
        format(panel$dates[ci[2] + 1], "%b %Y")))
    p_tail <- (1 + sum(sim_coverage[, j] <= observed[j])) / (reps + 1)
    cells <- c(cells, two_line(pct(observed[j]), sprintf("(%.3f)", p_tail)),
        two_line(pct(mean(random[, j])), paste0("[",
        pct(quantile(random[, j], 0.025, type = 8)), ", ",
        pct(quantile(random[, j], 0.975, type = 8)), "]")),
        two_line(pct(mean(sim_coverage[, j])), paste0("[",
        pct(quantile(sim_coverage[, j], 0.025, type = 8)), ", ",
        pct(quantile(sim_coverage[, j], 0.975, type = 8)), "]")))
    rows <- c(rows, tex_row(cells, "[3pt]"))
  }
  sim_head <- "\\begin{tabular}[b]{@{}c@{}}Common-date\\\\simulation\\end{tabular}"
  write_table("table04_common", "\\begin{tabular}{llccc}",
      paste0("New regime & 95\\% interval & Observed & Random & ", sim_head, "\\\\"),
      rows, "\\end{tabular}")
}

# -----------------------------------------------------------------------------
# panel_block()
#   What it does : Rows of several panels with an italic title above each and
#                  a mid rule between panels.
#   Called by    : table05_coefficients(); table07_09_common()
#   Feeds        : Tables 5, 7-9
# -----------------------------------------------------------------------------
panel_block <- function(rows, titles, size, ncol = 4) {
  out <- c()
  for (g in 1:length(titles)) {
    if (g > 1) out <- c(out, "\\midrule")
    out <- c(out, paste0("\\multicolumn{", ncol, "}{l}{\\textit{", titles[g], "}}\\\\"),
        rows[(g - 1) * size + (1:size)])
  }
  out
}

# -----------------------------------------------------------------------------
# table05_coefficients()   Table 5: regional means (Panel A) and medians
#                          (Panel B, the baseline) of the local-date coefficients
#   Called by: run_all.R
# -----------------------------------------------------------------------------
table05_coefficients <- function(bl) {
  rows <- c()
  for (g in 1:5) {
    cells <- c()
    for (k in 1:3) cells <- c(cells, stacked_cell(bl$period_mean[g, k], bl$period_se[g, k],
        bl$period_p[g, k]))
    rows <- c(rows, tex_row(c(short_names[g], cells), "[4pt]"))
  }
  for (g in 1:5) {
    cells <- c()
    for (k in 1:3) cells <- c(cells, stacked_cell(bl$median_point[g, k], bl$median_se[g, k],
        bl$median_p[g, k]))
    rows <- c(rows, tex_row(c(short_names[g], cells), "[4pt]"))
  }
  write_table("table05_coefficients", "\\begin{tabular}{lccc}",
      "Region & Spatial coefficient & Sales slope & Employment slope\\\\",
      panel_block(rows, c("Panel A: regional mean-group estimates",
      "Panel B: regional medians"), 5), "\\end{tabular}")
}

# -----------------------------------------------------------------------------
# table06_fullcoef()   Table 6: distribution of local-date coefficients across
#                      LAD-months
#   Called by: run_all.R
# -----------------------------------------------------------------------------
table06_fullcoef <- function(bl) {
  rows <- c()
  names3 <- c("Spatial", "Sales", "Employment")
  for (k in 1:3) {
    cells <- c()
    for (j in 1:6) cells <- c(cells, stacked_cell(bl$pooled_point[k, j], bl$pooled_se[k, j],
        bl$pooled_p[k, j]))
    rows <- c(rows, tex_row(c(names3[k], cells), "[4pt]"))
  }
  write_table("table06_fullcoef",
      "\\begin{tabular*}{\\linewidth}{@{\\extracolsep{\\fill}}lcccccc@{}}",
      "Coefficient & Mean & Median & P25 & P75 & P10 & P90\\\\", rows, "\\end{tabular*}")
}

# -----------------------------------------------------------------------------
# table07_09_common()   Tables 7 (regional median coefficients), 8 (effects of
#                       sales) and 9 (effects of employment) within the common
#                       regimes; effects are medians over the LADs of each
#                       region, with every LAD given its region's median
#                       spatial coefficient
#   Called by: run_all.R
# -----------------------------------------------------------------------------
table07_09_common <- function(panel, bc) {
  n_reg <- length(bc$starts)
  regime <- paste(format(panel$dates[bc$starts], "%Y-%m"), format(panel$dates[bc$ends], "%Y-%m"),
      sep = " to ")
  titles <- paste0("Panel ", LETTERS[1:5], ": ", c("North / Scotland", "Midlands / Wales",
      "East of England", "London", "South"))
  wide <- "\\begin{tabular*}{\\linewidth}{@{\\extracolsep{\\fill}}lccc@{}}"
  # Table 7
  rows <- c()
  for (g in 1:5) for (j in 1:n_reg) {
    cells <- c()
    for (k in 1:3) cells <- c(cells, coef_cell(bc$regime_point[j, g, k], bc$regime_se[j, g, k],
        bc$regime_p[j, g, k]))
    rows <- c(rows, tex_row(c(regime[j], cells)))
  }
  write_table("table07_commoncoef", wide, "Common regime & Spatial & Sales & Employment\\\\",
      panel_block(rows, titles, n_reg), "\\end{tabular*}")
  # Tables 8 and 9
  for (k in 1:2) {
    rows <- c()
    for (g in 1:5) for (j in 1:n_reg) {
      cells <- c()
      for (h in 1:3) cells <- c(cells, stacked_cell(bc$effect_point[j, g, k, h],
          bc$effect_se[j, g, k, h], bc$effect_p[j, g, k, h]))
      rows <- c(rows, tex_row(c(regime[j], cells)))
    }
    name <- if (k == 1) "table08_common_sales" else "table09_common_employment"
    write_table(name, wide, "Common regime & Direct & Indirect & Total\\\\",
        panel_block(rows, titles, n_reg), "\\end{tabular*}")
  }
}

# -----------------------------------------------------------------------------
# table10_spatial_weights()   Table 10: regional medians under three spatial
#                             weights matrices
#   Called by: run_all.R
# -----------------------------------------------------------------------------
table10_spatial_weights <- function(bl, wa) {
  rows <- c()
  for (g in 1:5) {
    cells <- c()
    for (k in 1:3) cells <- c(cells, stacked_cell(bl$median_point[g, k], bl$median_se[g, k],
        bl$median_p[g, k]))
    for (w in c("contiguity", "knn")) {
      for (k in 1:3) cells <- c(cells, stacked_cell(wa[[w]]$point[g, k], wa[[w]]$se[g, k],
          wa[[w]]$p[g, k]))
    }
    rows <- c(rows, tex_row(c(short_names[g], cells), "[4pt]"))
  }
  write_table("table10_spatial_weights",
      "\\begin{tabular*}{\\linewidth}{@{\\extracolsep{\\fill}}lccccccccc@{}}",
      c(" & \\multicolumn{3}{c}{Inverse distance} & \\multicolumn{3}{c}{Queen contiguity} & \\multicolumn{3}{c}{5 nearest neighbours}\\\\",
      "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\\cmidrule(lr){8-10}",
      "Region & Spatial & Sales & Empl. & Spatial & Sales & Empl. & Spatial & Sales & Empl.\\\\"),
      rows, "\\end{tabular*}")
}

# -----------------------------------------------------------------------------
# table11_min_regime()   Table 11 (Appendix A.1): pooled sequential tests and
#                        common dates with 24-, 12- and 6-month minimum regimes
#   Called by: run_all.R
# -----------------------------------------------------------------------------
table11_min_regime <- function(panel, settings, selection24, common_dp24, checks) {
  runs <- list(list(sel = selection24, dp = common_dp24),
      list(sel = checks$m12$selection, dp = checks$m12$common_dp),
      list(sel = checks$m6$selection, dp = checks$m6$common_dp))
  steps <- 0:selection24$m                       # up to the first non-rejection
  rows <- c()
  for (l in steps) {
    cells <- c(sprintf("%d vs %d", l, l + 1),
        sprintf("%.3f", selection24$tests$critical_value[l + 1]))
    for (r in runs) cells <- c(cells, sprintf("%.3f", r$sel$tests$statistic[l + 1]),
        fmt_p(r$sel$tests$p_value[l + 1]))
    rows <- c(rows, tex_row(cells))
  }
  # Summary rows: trimming at the last step shown, selected count, dates, SSR
  span <- function(x) paste0("\\multicolumn{2}{c}{", x, "}")
  last <- max(steps) + 1
  trim <- c()
  selected <- c()
  dates_top <- c()
  dates_bottom <- c()
  ssr <- c()
  for (r in runs) {
    t <- r$sel$tests[last, ]
    trim <- c(trim, span(sprintf("%.0f--%.0f", 100 * t$trim_min, 100 * t$trim_max)))
    selected <- c(selected, span(r$sel$m))
    d <- format(panel$dates[r$dp$breaks[[selection24$m + 1]] + 1], "%b %Y")
    dates_top <- c(dates_top, span(paste0(paste(d[1:2], collapse = ", "), ",")))
    dates_bottom <- c(dates_bottom, span(paste(d[3:4], collapse = ", ")))
    ssr <- c(ssr, span(sprintf("%.3f", r$dp$loss[selection24$m + 1])))
  }
  body <- c(rows, "\\midrule",
      tex_row(c(sprintf("Trimming at %d vs %d (\\%%)", last - 1, last), "", trim)),
      tex_row(c("Selected common breaks", "", selected)),
      tex_row(c(sprintf("Dates with %d breaks", selection24$m), "", dates_top)),
      tex_row(c("", "", dates_bottom)),
      tex_row(c(sprintf("Dating SSR with %d breaks", selection24$m), "", ssr)))
  write_table("table11_min_regime",
      "\\begin{tabular*}{\\linewidth}{@{\\extracolsep{\\fill}}lccccccc@{}}",
      c(sprintf(" & & \\multicolumn{2}{c}{%d-month minimum} & \\multicolumn{2}{c}{12-month minimum} & \\multicolumn{2}{c}{6-month minimum}\\\\",
      settings$min_length), "\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}\\cmidrule(lr){7-8}",
      "Test & 5\\% CV & Wald & $p$ & Wald & $p$ & Wald & $p$\\\\"),
      body, "\\end{tabular*}")
}

# -----------------------------------------------------------------------------
# table12_counts()   Table 12 (Appendix A.2): local break counts that the LAD-level
#                    sequential tests and the information criteria would select
#                    (the baseline gives every LAD four breaks)
#   Called by: run_all.R
# -----------------------------------------------------------------------------
table12_counts <- function(panel, search, cap) {
  gi <- region_index(panel)
  rows <- c()
  for (g in 1:5) {
    m <- search$m_seq[gi[[g]]]
    rows <- c(rows, tex_row(c(short_names[g], length(m), sum(m), sprintf("%.2f", mean(m)),
        sprintf("%.1f", 100 * mean(m == 0)), sprintf("%.1f", 100 * mean(m == cap)),
        sprintf("%.2f", mean(search$m_ic[gi[[g]]])), sprintf("%.2f", mean(search$m_ic5[gi[[g]]])))))
  }
  write_table("table12_counts", "\\begin{tabular}{lrrrrrrr}",
      c(" & & \\multicolumn{4}{c}{Sequential tests} & \\multicolumn{2}{c}{Information criterion}\\\\",
      "\\cmidrule(lr){3-6}\\cmidrule(lr){7-8}",
      "Region & LADs & Breaks & Mean & None (\\%) & At cap (\\%) & Mean & Stronger\\\\"),
      rows, "\\end{tabular}")
}

# -----------------------------------------------------------------------------
# table13_fit()   Table 13 (Appendix B): dating and final SSRs (Panel A) and
#                 the penalised comparison of the chronologies (Panel B)
#   Called by: run_all.R
# -----------------------------------------------------------------------------
table13_fit <- function(ssr_dating, ssr_final, pen) {
  labels <- c("Local CCEX-IV: four breaks (baseline)", "Common CCEX-IV: four breaks (baseline)",
      "Common nonspatial CCE: selected", "Local CCEX-IV: five breaks",
      "Common CCEX-IV: five breaks", "Common nonspatial CCE: five breaks")
  rows <- c()
  for (j in 1:6) rows <- c(rows, tex_row(c(labels[j], sprintf("%.3f", ssr_dating[j]),
      sprintf("%.3f", ssr_final[j]))))
  signed <- function(x) ifelse(x < 0, paste0("$-$", formatC(abs(round(x)), format = "d",
      big.mark = ",")), formatC(round(x), format = "d", big.mark = ","))
  pen_rows <- c(tex_row(c("Baseline chronologies (four breaks)", signed(pen$selected))),
      tex_row(c("Five breaks in both", signed(pen$five_breaks))),
      tex_row(c("Break counts re-selected by the criterion", signed(pen$reselected))),
      tex_row(c("\\quad Re-selected local breaks (mean)", sprintf("%.2f", pen$reselected_local_mean))),
      tex_row(c("\\quad Re-selected common breaks", sprintf("%d", pen$reselected_common))),
      tex_row(c("LADs favouring local dates (baseline)", sprintf("%d", pen$lads_favouring_local))))
  write_table("table13_fit", "\\begin{tabular}{lrr}",
      "Specification & OLS dating SSR & Final SSR\\\\",
      c("\\multicolumn{3}{l}{\\textit{Panel A: residual sums of squares}}\\\\", rows,
      "\\midrule",
      "\\multicolumn{3}{l}{\\textit{Panel B: penalised comparison, IC(local) $-$ IC(common)}}\\\\",
      "Comparison & BIC & LWZ\\\\", "\\midrule", pen_rows), "\\end{tabular}")
}
