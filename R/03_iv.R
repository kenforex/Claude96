# =============================================================================
# 03_iv.R
# Step 3. CCEX-IV estimation of equation (1) within each regime (Section 3).
#
# Within a regime the two cross-sectional averages (sbar, ebar) are partialled
# out, and the coefficients on the spatial lag Wy, sales and employment are
# estimated by 2SLS with instruments [X, WX, W^2X].
# Feeds: Tables 5-10 and 13; Figures 4-6.
# =============================================================================

# -----------------------------------------------------------------------------
# residualize()
#   What it does : Residuals of a vector or matrix after an OLS projection on
#                  the columns of F (QR decomposition).
#   Called by    : fit_iv_regime(); fit_panel_cce()
#   Feeds        : Tables 5-10 and 13; Figures 4-6
# -----------------------------------------------------------------------------
residualize <- function(a, F) qr.resid(qr(F, tol = 1e-10), as.matrix(a))

# -----------------------------------------------------------------------------
# spatial_instruments()
#   What it does : Instrument blocks [X, WX, W^2X] for the spatial lag, as a
#                  list of three T x N x 2 arrays (X = sales and employment).
#   Called by    : fit_panel_iv()
#   Feeds        : Tables 5-10 and 13; Figures 4-6
# -----------------------------------------------------------------------------
spatial_instruments <- function(panel) {
  WX <- panel$X
  W2X <- panel$X
  for (k in 1:2) {
    WX[, , k] <- panel$X[, , k] %*% t(panel$W)
    W2X[, , k] <- WX[, , k] %*% t(panel$W)
  }
  list(panel$X, WX, W2X)
}

# -----------------------------------------------------------------------------
# fit_iv_regime()
#   What it does : CCEX-IV for one LAD and one regime: partial out F, then 2SLS
#                  of y on d = (Wy, s, e) with instruments z. Returns the three
#                  coefficients, the residuals and an orthonormal basis of the
#                  instruments (reused by the bootstrap). Returns NA
#                  coefficients if the first-stage fitted values are (nearly)
#                  collinear, which does not happen with 24-month regimes.
#   Called by    : fit_panel_iv()
#   Feeds        : Tables 5-10 and 13; Figures 4-6
# -----------------------------------------------------------------------------
fit_iv_regime <- function(y, d, z, F) {
  yf <- residualize(y, F)
  df <- residualize(d, F)
  zf <- residualize(z, F)
  qz <- qr(zf, tol = 1e-10)
  qb <- qr.Q(qz)[, 1:qz$rank, drop = FALSE]    # orthonormal basis of the instruments
  dhat <- qb %*% crossprod(qb, df)             # first-stage fitted values
  qd <- qr(dhat, tol = 1e-10)
  if (qd$rank < ncol(d) || kappa(dhat) > 1e6) {
    return(list(coef = rep(NA, 3),
        residual = as.numeric(qr.resid(qr(df, tol = 1e-10), yf)), z_basis = NULL))
  }
  bhat <- as.numeric(qr.coef(qd, yf))
  e <- as.numeric(yf - df %*% bhat)
  list(coef = bhat, residual = e, z_basis = qb)
}

# -----------------------------------------------------------------------------
# fit_panel_iv()
#   What it does : fit_iv_regime() for every LAD and every regime of a given set
#                  of break dates (breaks[[i]] for LAD i). Returns monthly
#                  coefficient paths (T x N x 3: spatial, sales, employment),
#                  residuals, and for each regime its months and instrument
#                  basis (reused by the bootstrap).
#   Called by    : run_all.R (Step 3); weights_analysis() [05]
#   Feeds        : Tables 5-10 and 13; Figures 4-6
# -----------------------------------------------------------------------------
fit_panel_iv <- function(panel, breaks) {
  Z <- spatial_instruments(panel)
  WY <- panel$Y %*% t(panel$W)
  coef <- array(NA, c(panel$T, panel$N, 3))
  residual <- matrix(NA, panel$T, panel$N)
  details <- list()
  for (i in 1:panel$N) {
    bounds <- c(0, breaks[[i]], panel$T)
    details[[i]] <- list()
    for (j in 1:(length(bounds) - 1)) {
      rows <- (bounds[j] + 1):bounds[j + 1]
      z <- cbind(Z[[1]][rows, i, ], Z[[2]][rows, i, ], Z[[3]][rows, i, ])
      d <- cbind(WY[rows, i], panel$X[rows, i, ])
      fit <- fit_iv_regime(panel$Y[rows, i], d, z, panel$F[rows, , drop = FALSE])
      for (k in 1:3) coef[rows, i, k] <- fit$coef[k]
      residual[rows, i] <- fit$residual
      details[[i]][[j]] <- list(rows = rows, z_basis = fit$z_basis)
    }
  }
  list(coef = coef, residual = residual, details = details)
}

# -----------------------------------------------------------------------------
# fit_panel_cce()
#   What it does : Nonspatial CCE comparator: y on sales and employment with
#                  sbar, ebar and ybar as cross-sectional averages, on a common
#                  partition. Returns the residuals.
#   Called by    : run_all.R (Step 3)
#   Feeds        : Table 13, Panel A (CCE rows)
# -----------------------------------------------------------------------------
fit_panel_cce <- function(panel, breaks) {
  residual <- matrix(NA, panel$T, panel$N)
  bounds <- c(0, breaks, panel$T)
  ybar <- rowMeans(panel$Y)
  for (i in 1:panel$N) {
    for (j in 1:(length(bounds) - 1)) {
      rows <- (bounds[j] + 1):bounds[j + 1]
      F <- cbind(panel$F[rows, , drop = FALSE], ybar[rows])
      y <- residualize(panel$Y[rows, i], F)
      x <- residualize(panel$X[rows, i, ], F)
      residual[rows, i] <- lm.fit(x, y)$residuals
    }
  }
  list(residual = residual)
}
