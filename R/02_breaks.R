# =============================================================================
# 02_breaks.R
# Step 2. Break-date search and sequential sup-Wald tests (Section 3).
#
# Dating regression for LAD i: y_it on q_it = (Wy_t, s_it, e_it, sbar_t, ebar_t),
# all five coefficients fitted by OLS within each candidate regime.
#   * segment_costs() + run_dp(): the SSR-minimising partition for 0-9 breaks
#     (Bai and Perron 2003), for each LAD (local dates) and for the sum over
#     LADs (common dates).
#   * select_common_breaks(): pooled sequential sup-Wald tests of a change in
#     the average (pooled) slopes, simulated 15%-trimming critical values
#     (Dai et al. 2026) -> the four common dates.
#   * sequential_local_breaks(): the same tests LAD by LAD with the tabulated
#     5% critical values of Bai and Perron (2003b), q = 3 (a check: the
#     baseline gives every LAD the number of common breaks, four).
# Feeds: all local and common dates, hence every table and figure.
# Table and figure numbers in this package are those of the paper
# (Tables 1-13, Figures 1-7).
# =============================================================================

# -----------------------------------------------------------------------------
# safe_solve()
#   What it does : solve(A) (the inverse) or solve(A, b), but returns NULL
#                  instead of stopping when A is singular.
#   Called by    : pooled_split_F(); regime_wald_side()
#   Feeds        : the sequential tests (number of common breaks; Tables 11, 12)
# -----------------------------------------------------------------------------
safe_solve <- function(A, b = NULL) {
  if (is.null(b)) out <- try(solve(A), silent = TRUE) else out <- try(solve(A, b), silent = TRUE)
  if (inherits(out, "try-error")) return(NULL)
  out
}

# -----------------------------------------------------------------------------
# segment_costs()
#   What it does : OLS residual sum of squares of the dating regression on every
#                  interval [a, b] of at least min_length months. All intervals
#                  are done at once with cumulative cross-products: after an
#                  orthonormal rotation, the SSR of an interval is the last
#                  Schur complement of its cross-product matrix. Returns a T x T
#                  matrix (Inf where an interval is too short).
#   Called by    : search_breaks() [02]; simulate_one_draw(),
#                  min_regime_checks() [06]
#   Feeds        : all local and common dates (every table and figure)
# -----------------------------------------------------------------------------
segment_costs <- function(y, design, min_length) {
  nt <- length(y)
  qx <- qr(design, tol = 1e-10)
  basis <- qr.Q(qx)[, 1:qx$rank, drop = FALSE]
  residual <- y - basis %*% crossprod(basis, y)
  Z <- cbind(basis, as.numeric(residual))
  p <- ncol(Z)
  # All admissible intervals (start, end)
  n_start <- nt - min_length + 1
  starts <- rep(1:n_start, nt - min_length + 2 - (1:n_start))
  ends <- c()
  for (s in 1:n_start) ends <- c(ends, (s + min_length - 1):nt)
  count <- length(starts)
  # Cross-products of every interval from cumulative sums
  G <- array(0, c(count, p, p))
  for (j in 1:p) for (k in j:p) {
    cs <- c(0, cumsum(Z[, j] * Z[, k]))
    G[, j, k] <- cs[ends + 1] - cs[starts]
    G[, k, j] <- G[, j, k]
  }
  # Gaussian elimination, all intervals in parallel; G[, p, p] ends as the SSR
  bad <- rep(FALSE, count)
  for (j in 1:(p - 1)) {
    pivot <- G[, j, j]
    bad <- bad | !is.finite(pivot) | pivot < 1e-12
    pivot[pivot < 1e-12 | !is.finite(pivot)] <- 1
    for (k in (j + 1):p) for (l in k:p) {
      G[, k, l] <- G[, k, l] - G[, j, k] * G[, j, l] / pivot
      G[, l, k] <- G[, k, l]
    }
  }
  loss <- G[, p, p]
  bad <- bad | !is.finite(loss) | loss < -1e-9
  # Rare numerically awkward intervals: direct QR regression instead
  for (j in which(bad)) {
    rows <- starts[j]:ends[j]
    loss[j] <- sum(qr.resid(qr(design[rows, , drop = FALSE], tol = 1e-10), y[rows])^2)
  }
  out <- matrix(Inf, nt, nt)
  out[cbind(starts, ends)] <- pmax(loss, 0)
  out
}

# -----------------------------------------------------------------------------
# run_dp()
#   What it does : Dynamic programme of Bai and Perron (2003): given the
#                  interval losses, the globally SSR-minimising partition for
#                  each number of breaks 0, ..., max_breaks. Breaks are stored
#                  as the last month of the old regime.
#   Called by    : search_breaks() [02]; simulate_one_draw(),
#                  min_regime_checks() [06]
#   Feeds        : all local and common dates (every table and figure)
# -----------------------------------------------------------------------------
run_dp <- function(cost, max_breaks, min_length) {
  nt <- nrow(cost)
  max_breaks <- min(max_breaks, nt %/% min_length - 1)
  value <- matrix(Inf, max_breaks + 1, nt)
  previous <- matrix(NA_integer_, max_breaks + 1, nt)
  value[1, ] <- cost[1, ]
  if (max_breaks > 0) {
    for (r in 2:(max_breaks + 1)) {
      for (end in (r * min_length):nt) {
        candidates <- ((r - 1) * min_length):(end - min_length)
        z <- value[r - 1, candidates] + cost[cbind(candidates + 1, end)]
        best <- which.min(z)
        value[r, end] <- z[best]
        previous[r, end] <- candidates[best]
      }
    }
  }
  # Trace back the break dates for every number of breaks
  paths <- list()
  for (m in 0:max_breaks) {
    b <- integer(m)
    end <- nt
    if (m > 0) {
      for (r in (m + 1):2) {
        b[r - 1] <- previous[r, end]
        end <- b[r - 1]
      }
    }
    paths[[m + 1]] <- b
  }
  list(loss = value[, nt], breaks = paths)
}

# -----------------------------------------------------------------------------
# ic_break_count()
#   What it does : Number of breaks chosen by the information criterion of
#                  equation (4): log(SSR/T) + k (m + 1) log(T) / T, with k = 3
#                  penalised coefficients (or k = 5 for the stronger version).
#                  With k = 2 it also picks the break count of the nonspatial
#                  CCE comparator.
#   Called by    : search_breaks()
#   Feeds        : Table 12 (information-criterion columns); Table 13 (CCE row);
#                  Appendix A.1 text
# -----------------------------------------------------------------------------
ic_break_count <- function(loss, nt, k = 3) {
  m <- 0:(length(loss) - 1)
  which.min(log(loss / nt) + (m + 1) * k * log(nt) / nt) - 1
}

# -----------------------------------------------------------------------------
# bp_critical_values()
#   What it does : 5% critical values of the sup F(l+1|l) test with 15%
#                  trimming from Bai and Perron (2003b, Econometrics Journal 6,
#                  72-78), for q = 1-10 changing coefficients and l = 0-9
#                  breaks under the null (the table shipped with the R package
#                  mbreaks). Returned divided by q, because every statistic in
#                  this package is Wald / q.
#   Called by    : search_breaks(); run_all.R (footnote value 6.047)
#   Feeds        : the sequential local break counts (Table 12; Figure 7;
#                  Section 4.2 and Appendix A text)
# -----------------------------------------------------------------------------
bp_critical_values <- function(l, q = 3) {
  # Rows q = 1, ..., 10; columns l = 0, ..., 9
  table_05 <- matrix(c(
      8.58, 10.13, 11.14, 11.83, 12.25, 12.66, 13.08, 13.35, 13.75, 13.89,
      11.47, 12.95, 14.03, 14.85, 15.29, 15.80, 16.16, 16.44, 16.77, 16.84,
      13.98, 15.72, 16.83, 17.61, 18.14, 18.74, 19.09, 19.41, 19.68, 19.77,
      16.19, 18.11, 18.93, 19.64, 20.19, 20.54, 21.21, 21.42, 21.72, 21.97,
      18.23, 19.91, 20.99, 21.71, 22.37, 22.77, 23.15, 23.42, 24.04, 24.42,
      20.08, 22.11, 23.04, 23.77, 24.43, 24.75, 24.96, 25.22, 25.61, 25.93,
      21.87, 24.17, 25.13, 26.03, 26.65, 27.06, 27.37, 27.90, 28.18, 28.36,
      23.70, 25.75, 26.81, 27.65, 28.48, 28.80, 29.08, 29.30, 29.50, 29.69,
      25.65, 27.66, 28.91, 29.67, 30.52, 30.96, 31.48, 31.77, 31.94, 32.33,
      27.03, 29.24, 30.45, 31.45, 32.12, 32.50, 32.84, 33.12, 33.22, 33.85),
      10, 10, byrow = TRUE)
  table_05[q, l + 1] / q
}

# -----------------------------------------------------------------------------
# simulate_supF()
#   What it does : 10,000 draws of the sup-Wald limit: the supremum of a
#                  normalised three-dimensional Brownian bridge over [0.15, 0.85]
#                  (fixed seed).
#   Called by    : select_common_breaks()
#   Feeds        : number of common breaks; CV and p columns of Table 11
# -----------------------------------------------------------------------------
simulate_supF <- function(draws = 10000, grid = 600, trim = 0.15, seed = 20260908) {
  set.seed(seed)
  time <- (1:grid) / grid
  use <- which(time >= trim & time <= 1 - trim)
  maxima <- numeric(draws)
  for (b in 1:draws) {
    brown <- apply(matrix(rnorm(grid * 3), grid, 3), 2, cumsum) / sqrt(grid)
    bridge <- brown - time %o% brown[grid, ]
    maxima[b] <- max(rowSums(bridge[use, , drop = FALSE]^2) /
        (3 * time[use] * (1 - time[use])))
  }
  maxima
}

# -----------------------------------------------------------------------------
# pooled_split_F()
#   What it does : Pooled Wald statistic (divided by three) for splitting one
#                  regime [left, right] after month `split`. On each side of
#                  the split, the factor proxies (sbar, ebar) are partialled out
#                  LAD by LAD, so their coefficients (the loadings) are
#                  LAD-specific and free on each side. One common slope vector
#                  on (Wy, s, e) is then estimated on each side by pooled OLS
#                  over all LADs and months; it is NOT LAD-specific. Under the
#                  random-coefficient assumption of Dai et al. (2026, Section
#                  2.3) this pooled slope estimates the average slope, so the
#                  statistic tests equal average slopes on the two sides. The
#                  variance is clustered by LAD (robust to heterogeneous slopes
#                  and to serial correlation within a LAD, including the
#                  covariance between the two sides), and the statistic is
#                  scaled by a degrees-of-freedom factor.
#   Called by    : select_common_breaks()
#   Feeds        : the common dates (every table and figure); Table 11
# -----------------------------------------------------------------------------
pooled_split_F <- function(panel, WY, left, split, right, total_regimes) {
  blocks <- list(left:split, (split + 1):right)
  score <- matrix(0, panel$N, 6)
  beta <- numeric(6)
  bread <- matrix(0, 6, 6)
  for (j in 1:2) {
    rows <- blocks[[j]]
    qf <- qr(panel$F[rows, , drop = FALSE])
    yy <- qr.resid(qf, panel$Y[rows, , drop = FALSE])
    d <- list(qr.resid(qf, WY[rows, , drop = FALSE]),
        qr.resid(qf, panel$X[rows, , 1]), qr.resid(qf, panel$X[rows, , 2]))
    G <- matrix(0, 3, 3)
    v <- numeric(3)
    for (k in 1:3) {
      v[k] <- sum(d[[k]] * yy)
      for (l in 1:3) G[k, l] <- sum(d[[k]] * d[[l]])
    }
    inv <- safe_solve(G)
    if (is.null(inv)) return(NA)
    b <- as.numeric(inv %*% v)
    ee <- yy
    for (k in 1:3) ee <- ee - d[[k]] * b[k]
    jj <- (3 * (j - 1) + 1):(3 * j)
    beta[jj] <- b
    bread[jj, jj] <- inv
    for (k in 1:3) score[, jj[k]] <- colSums(d[[k]] * ee)
  }
  V <- bread %*% crossprod(score) %*% bread
  R <- cbind(diag(3), -diag(3))
  delta <- R %*% beta
  VR <- R %*% V %*% t(R)
  VR_delta <- safe_solve(VR, delta)
  if (is.null(VR_delta)) return(NA)
  w <- as.numeric(crossprod(delta, VR_delta))
  df_adjust <- (panel$N * (panel$T - total_regimes * 3) - total_regimes * 3) /
    (panel$N * panel$T)
  df_adjust * w / 3
}

# -----------------------------------------------------------------------------
# select_common_breaks()
#   What it does : Pooled sequential tests of l against l + 1 common breaks on
#                  the global l-break common partitions (l = 0, ..., 8). A
#                  regime of n months is split only at dates leaving at least
#                  h = max(min_length, ceiling(0.15 n)) months on each side.
#                  Critical value at step l: the (0.95)^(1/(l+1)) quantile of
#                  the simulated sup-Wald draws. The selected number is the
#                  first non-rejection. Also returns, for every step, the
#                  number of regimes searched and the effective trimming.
#   Called by    : run_all.R (Step 2); min_regime_checks() [06]
#   Feeds        : the four common dates (every table and figure); Table 11;
#                  footnote of Section 4.2 and Appendix A.1 text
# -----------------------------------------------------------------------------
select_common_breaks <- function(panel, common_dp, WY, min_length, max_breaks, alpha = 0.05) {
  draws <- simulate_supF()
  tests <- data.frame()
  for (m in 0:(max_breaks - 1)) {
    bounds <- c(0, common_dp$breaks[[m + 1]], panel$T)
    stats <- c()
    trims <- c()
    for (j in 1:(m + 1)) {
      s <- bounds[j] + 1
      e <- bounds[j + 1]
      n <- e - s + 1
      h <- max(min_length, ceiling(0.15 * n))
      if (n < 2 * h) next                      # regime too short to be split
      trims <- c(trims, h / n)
      for (tau in (s + h - 1):(e - h)) {
        stats <- c(stats, pooled_split_F(panel, WY, s, tau, e, m + 2))
      }
    }
    cv <- as.numeric(quantile(draws, (1 - alpha)^(1 / (m + 1)), type = 8))
    stat <- NA
    if (any(is.finite(stats))) stat <- max(stats, na.rm = TRUE)
    p_one <- (1 + sum(draws > stat)) / (length(draws) + 1)
    tests <- rbind(tests, data.frame(null_breaks = m, statistic = stat,
        critical_value = cv, p_value = 1 - (1 - p_one)^(m + 1),
        reject = is.finite(stat) && stat > cv,
        regimes = m + 1, regimes_tested = length(trims),
        trim_min = if (length(trims) > 0) min(trims) else NA,
        trim_max = if (length(trims) > 0) max(trims) else NA))
    cat("  Common sup-Wald", m, "vs", m + 1, ": statistic", round(stat, 3),
        "critical value", round(cv, 3), "\n")
  }
  if (all(tests$reject)) selected <- max_breaks else selected <- which(!tests$reject)[1] - 1
  list(m = selected, breaks = common_dp$breaks[[selected + 1]], tests = tests)
}

# -----------------------------------------------------------------------------
# regime_wald_side()
#   What it does : One side of every candidate split of a regime. For each
#                  left-block size k (left = TRUE: months 1..k; FALSE: k+1..n)
#                  and each outcome column: OLS coefficients on Z, residuals,
#                  and the Bartlett HAC variance (four lags) of the three
#                  tested slopes.
#   Called by    : regime_wald()
#   Feeds        : local break counts (see sequential_local_breaks())
# -----------------------------------------------------------------------------
regime_wald_side <- function(y, Z, k, left, zz, lags) {
  n <- nrow(Z)
  q <- ncol(Z)
  B <- ncol(y)
  K <- length(k)
  bartlett <- 1 - (1:lags) / (lags + 1)
  pair_t <- as.vector(t(matrix(1:(q * q), q)))  # position of (c, a) for pair (a, c)
  theta <- array(NA, c(q, K, B))
  A <- array(NA, c(3, q, K))
  in_block <- matrix(0, n, K)
  for (j in 1:K) {
    if (left) rows <- 1:k[j] else rows <- (k[j] + 1):n
    in_block[rows, j] <- 1
    Gi <- safe_solve(crossprod(Z[rows, , drop = FALSE]))
    if (is.null(Gi)) next
    A[, , j] <- Gi[1:3, ]
    theta[, j, ] <- Gi %*% crossprod(Z[rows, , drop = FALSE], y[rows, , drop = FALSE])
  }
  # Residuals of every (candidate, series) pair; zero outside the block
  u <- (y[, rep(1:B, each = K), drop = FALSE] - Z %*% matrix(theta, q)) *
    as.numeric(in_block)
  u[is.na(u)] <- 0
  # HAC "meat" for every column: sum over t of z_t z_{t-l}' u_t u_{t-l}
  omega <- crossprod(zz[[1]], u * u)
  for (l in 1:lags) {
    g <- crossprod(zz[[l + 1]], u[(l + 1):n, , drop = FALSE] * u[1:(n - l), , drop = FALSE])
    omega <- omega + bartlett[l] * (g + g[pair_t, , drop = FALSE])
  }
  # Variance of the three tested slopes: A Omega A'
  Ab <- A[, , rep(1:K, B), drop = FALSE]
  V <- array(0, c(3, 3, K * B))
  for (r in 1:3) for (c0 in r:3) {
    v <- 0
    for (a in 1:q) for (c1 in 1:q) v <- v + Ab[r, a, ] * omega[(c1 - 1) * q + a, ] * Ab[c0, c1, ]
    V[r, c0, ] <- v
    V[c0, r, ] <- v
  }
  list(b = matrix(theta[1:3, , , drop = FALSE], 3), V = V)
}

# -----------------------------------------------------------------------------
# regime_wald()
#   What it does : LAD-level Wald statistic W(k) for every admissible
#                  split k of one regime and every outcome column of y:
#                  W(k) = (b_L - b_R)' (V_L + V_R)^(-1) (b_L - b_R) / 3, where
#                  b are the slopes on (Wy, s, e) with the factor-proxy
#                  coefficients free on each side, and V their HAC variances.
#                  Z holds the regressors with the three tested ones first.
#   Called by    : lad_sup_wald()
#   Feeds        : local break counts (see sequential_local_breaks())
# -----------------------------------------------------------------------------
regime_wald <- function(y, Z, k, lags = 4) {
  y <- as.matrix(y)
  n <- nrow(Z)
  q <- ncol(Z)
  K <- length(k)
  # Products z_{t,a} z_{t-l,c} for all pairs (a, c) and lags l = 0, ..., 4
  zz <- list()
  for (l in 0:lags) {
    lead <- Z[(l + 1):n, , drop = FALSE]
    lagged <- Z[1:(n - l), , drop = FALSE]
    zz[[l + 1]] <- lead[, rep(1:q, q), drop = FALSE] * lagged[, rep(1:q, each = q), drop = FALSE]
  }
  L <- regime_wald_side(y, Z, k, TRUE, zz, lags)
  R <- regime_wald_side(y, Z, k, FALSE, zz, lags)
  d <- L$b - R$b
  S <- L$V + R$V
  # d' S^(-1) d through the cofactors of the 3 x 3 matrix S
  c11 <- S[2, 2, ] * S[3, 3, ] - S[2, 3, ]^2
  c12 <- S[1, 3, ] * S[2, 3, ] - S[1, 2, ] * S[3, 3, ]
  c13 <- S[1, 2, ] * S[2, 3, ] - S[1, 3, ] * S[2, 2, ]
  c22 <- S[1, 1, ] * S[3, 3, ] - S[1, 3, ]^2
  c23 <- S[1, 2, ] * S[1, 3, ] - S[1, 1, ] * S[2, 3, ]
  c33 <- S[1, 1, ] * S[2, 2, ] - S[1, 2, ]^2
  det <- S[1, 1, ] * c11 + S[1, 2, ] * c12 + S[1, 3, ] * c13
  w <- (d[1, ]^2 * c11 + d[2, ]^2 * c22 + d[3, ]^2 * c33 + 2 * (d[1, ] * d[2, ] * c12 +
      d[1, ] * d[3, ] * c13 + d[2, ] * d[3, ] * c23)) / det / 3
  w[!is.finite(w) | !(det > 0)] <- NA
  matrix(w, K, ncol(y))
}

# -----------------------------------------------------------------------------
# lad_sup_wald()
#   What it does : sup-F(l+1|l) statistic for one LAD: the largest W(k) over the
#                  admissible splits of all regimes of an l-break partition.
#                  Same splitting rule as the pooled test. y may hold several
#                  series in columns (one statistic per column).
#   Called by    : sequential_local_breaks()
#   Feeds        : local break counts (Table 12; Figure 7)
# -----------------------------------------------------------------------------
lad_sup_wald <- function(y, X, Fm, breaks, min_length) {
  y <- as.matrix(y)
  bounds <- c(0, breaks, nrow(y))
  stat <- rep(-Inf, ncol(y))
  for (j in 1:(length(bounds) - 1)) {
    s <- bounds[j] + 1
    e <- bounds[j + 1]
    n <- e - s + 1
    h <- max(min_length, ceiling(0.15 * n))
    if (n < 2 * h) next
    # Orthonormal basis of [F, X] on the regime, tested columns first. The
    # statistic does not depend on this rotation; it only helps numerically.
    Z <- qr.Q(qr(cbind(Fm, X)[s:e, , drop = FALSE]))[, c(3:5, 1:2), drop = FALSE]
    w <- regime_wald(y[s:e, , drop = FALSE], Z, h:(n - h))
    w[is.na(w)] <- -Inf
    stat <- pmax(stat, apply(w, 2, max))
  }
  stat
}

# -----------------------------------------------------------------------------
# sequential_local_breaks()
#   What it does : Number of local breaks for one LAD: starting from zero, test
#                  l against l + 1 breaks on the LAD's own global partitions
#                  (dp from run_dp()) and stop at the first non-rejection or
#                  at max_breaks. cv[l + 1] is the critical value of step l.
#   Called by    : search_breaks() [02]; simulate_one_draw(),
#                  min_regime_checks() [06]
#   Feeds        : the break counts of the sequential tests, used as a check
#                  on the baseline of four breaks per LAD: Table 12;
#                  Figure 7; Section 4.2 and Appendices A.1, A.2 text
# -----------------------------------------------------------------------------
sequential_local_breaks <- function(y, X, Fm, dp, min_length, max_breaks, cv) {
  max_breaks <- min(max_breaks, length(dp$breaks) - 1)
  m <- 0
  while (m < max_breaks) {
    stat <- lad_sup_wald(y, X, Fm, dp$breaks[[m + 1]], min_length)
    if (!is.finite(stat) || stat <= cv[m + 1]) break
    m <- m + 1
  }
  m
}

# -----------------------------------------------------------------------------
# search_breaks()
#   What it does : Step 2 for the observed data. For every LAD: interval losses,
#                  global partitions with 0-9 breaks, and the break count from
#                  the sequential tests (and from the two information
#                  criteria). Sums the losses over LADs for the common search,
#                  and runs the date search of the nonspatial CCE comparator.
#   Called by    : run_all.R (Step 2)
#   Feeds        : every table and figure (Table 13 for the CCE comparator)
# -----------------------------------------------------------------------------
search_breaks <- function(panel, settings) {
  ml <- settings$min_length
  mb <- settings$max_breaks
  WY <- panel$Y %*% t(panel$W)
  common_cost <- matrix(0, panel$T, panel$T)
  cce_cost <- matrix(0, panel$T, panel$T)
  costs <- list()
  local_dp <- list()
  for (i in 1:panel$N) {
    design <- cbind(WY[, i], panel$X[, i, ], panel$F)
    costs[[i]] <- segment_costs(panel$Y[, i], design, ml)
    common_cost <- common_cost + costs[[i]]
    local_dp[[i]] <- run_dp(costs[[i]], mb, ml)
    cce_design <- cbind(panel$X[, i, ], panel$F, rowMeans(panel$Y))
    cce_cost <- cce_cost + segment_costs(panel$Y[, i], cce_design, ml)
    if (i %% 50 == 0) cat("  Date search:", i, "of", panel$N, "LADs\n")
  }
  common_dp <- run_dp(common_cost, mb, ml)
  cce_dp <- run_dp(cce_cost, mb, ml)
  # Break counts by the sequential tests (a check) and by equation (4)
  cv <- bp_critical_values(0:(mb - 1))
  m_seq <- parallel::mclapply(1:panel$N, function(i)
    sequential_local_breaks(panel$Y[, i], cbind(WY[, i], panel$X[, i, ]), panel$F,
        local_dp[[i]], ml, mb, cv), mc.cores = settings$cores)
  m_ic <- m_ic5 <- numeric(panel$N)
  for (i in 1:panel$N) {
    m_ic[i] <- ic_break_count(local_dp[[i]]$loss, panel$T, 3)
    m_ic5[i] <- ic_break_count(local_dp[[i]]$loss, panel$T, 5)
  }
  list(WY = WY, costs = costs, local_dp = local_dp, common_dp = common_dp, cce_dp = cce_dp,
      m_seq = unlist(m_seq), m_ic = m_ic, m_ic5 = m_ic5, local_cv = cv,
      cce_m = ic_break_count(cce_dp$loss / panel$N, panel$T, 2))
}

# -----------------------------------------------------------------------------
# breaks_with_m()
#   What it does : Each LAD's break dates from its own global partition with
#                  m[i] breaks (m may be a single number for all LADs).
#   Called by    : run_all.R; min_regime_checks() [06]
#   Feeds        : the baseline local dates (four breaks per LAD): Tables 4-6,
#                  10, 13, Figures 2, 4-6; and the sequential and
#                  information-criterion dates (Figure 7, Appendix A text)
# -----------------------------------------------------------------------------
breaks_with_m <- function(local_dp, m) {
  if (length(m) == 1) m <- rep(m, length(local_dp))
  out <- list()
  for (i in 1:length(local_dp)) out[[i]] <- local_dp[[i]]$breaks[[m[i] + 1]]
  out
}
