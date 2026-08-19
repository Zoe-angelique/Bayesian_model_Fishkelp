# tdr_sampler_v7.R
# Vendored from FishKelp commit 8c46ebb0f3d31857e468376778833e933f6e4f5d.
# Source SHA-256: 49e53b33f7e83afecf67213c8b1cbc78fe91ab431b27e3f4dd5b5751540584e9.
# Keep this file aligned with data/fishkelp_v7/model/tdr_model_data_v7.rds.
# Shared functions for the dependency model (writing/methods/Kelp_Dependency_Model_v7.tex), on the evidentiary-ladder schema,
# the evidentiary-ladder response. The response is the highest ladder rung a source supports for a
# species on a route, with C ordered categories (rung 0 = no link, the informative negative).
# Sourced by code/model_v7/02_fit_m1_v7.R (the fit), code/model_v7/03_fit_m0_v7.R and
# code/model_v7/04_sensitivity_v7.R.
# (the region-only supplement), and code/model_v7/04_sensitivity_v7.R (sensitivity).
#
# Hierarchical ordered-PROBIT model fitted by Albert-Chib Gibbs sampling (compilation-free path; the
# canonical ordered-probit Stan reference is stan/dependency_v7.stan). The model has, on the latent
# scale,
#     z_j ~ Normal( X_j' beta_{d(j)} + u_{i(j),d(j)} + a_{s(j)} , sigma_j^2 ),  y_j = sum_c 1{z_j > kappa_{d,c}}
# with
#   * ROUTE-SPECIFIC fixed effects beta_d (region + standardised/contrast-coded species traits), so a
#     trait may act differently on the trophic and habitat routes; per-column regularising prior.
#   * correlated trophic/habitat species effects u_i ~ N2(0, Sigma) (inverse-Wishart), giving the
#     residual trophic-habitat correlation rho_TH AFTER the shared traits are in the location.
#   * a route-specific BASELINE INTERCEPT (a column of the design X) plus ordered cut-points kappa_d
#     with the first anchored at zero for identification (kappa[d,1] = 0). The intercept carries the
#     baseline location, so the reference profile is not forced to a fixed None probability, and the
#     anchor removes the location/cut-point ridge and keeps the sampler mixing; the free cut-points
#     use a weakly-informative, UNBIASED N(0, cut_sd) prior.
#   * OPTIONALLY a source random intercept a_s ~ N(0, tau_s^2) (use_source = TRUE) to absorb the
#     non-independence of multiple records from one study. It is not separately identifiable from the
#     species effect given the sparse positive signal (it degrades mixing without shifting the
#     reaches), so it is OFF by default and used only as a robustness check.
#   * Knowledge Specificity as the fixed per-row latent scale sigma (High most precise), i.e. a
#     per-record precision weight 1/sigma^2 on the location; it never enters the location itself.
# Convergence is assessed with rank-normalised split-Rhat and bulk effective sample size (Vehtari et
# al. 2021), implemented in mcmc_diag() below.

# Vectorised truncated-normal draw on (lo, hi], robust to +/-Inf bounds.
rtn <- function(mean, sd, lo, hi) {
  pl <- pnorm(lo, mean, sd); ph <- pnorm(hi, mean, sd)
  u  <- runif(length(mean), pl, ph)
  u  <- pmin(pmax(u, 1e-12), 1 - 1e-12)
  q  <- qnorm(u, mean, sd)
  pmin(pmax(q, lo), hi)
}

# Map a latent z to its ordinal category 1..C, given the per-route cut-point matrix kappa (D x C-1).
z_to_cat <- function(z, dimension, kappa, C) {
  cat <- rep(1L, length(z))
  for (m in 1:(C - 1L)) cat <- cat + as.integer(z > kappa[cbind(dimension, m)])
  cat
}

# One MCMC chain.
#   prior_scale : per-column SD of the N(0, prior_scale) prior on each route's beta
#                 (current M1: intercept, region, and invertebrate 2.5; traits 0.5)
#   nu0, tau0   : inverse-Wishart df and target between-species SD (S0 = tau0^2 (nu0-D-1) I)
#   cut_sd      : weakly-informative N(0, cut_sd) prior on the free cut-points
#   as0, bs0    : inverse-gamma prior on the source variance tau_s^2 (shape, rate)
gibbs_chain <- function(dat, y, dimension, species, sigma, X, source = dat$source,
                        niter = 6000, nburn = 2500, seed = 1,
                        nu0 = 40, tau0 = 0.7, prior_scale = NULL,
                        cut_sd = 5, as0 = 3, bs0 = 0.5, use_source = FALSE) {
  set.seed(seed)
  N <- length(y); I <- dat$I; D <- dat$D; P <- ncol(X); C <- dat$C
  S <- dat$S; if (is.null(source)) source <- rep(1L, N)
  if (is.null(prior_scale)) prior_scale <- rep(1, P)
  s2 <- sigma^2; w <- 1 / s2

  P0     <- diag(1 / prior_scale^2, P)         # per-column prior precision on beta
  S0     <- diag(D) * tau0^2 * (nu0 - D - 1)
  cut_mu <- 0                                  # free cut-point prior N(0, cut_sd), weakly-informative

  rows_d <- lapply(1:D, function(d) which(dimension == d))
  # per-route posterior precision of beta is fixed (X, w, P0 fixed): precompute the covariance.
  Vbeta  <- lapply(1:D, function(d) {
    Xd <- X[rows_d[[d]], , drop = FALSE]; wd <- w[rows_d[[d]]]
    solve(crossprod(Xd * wd, Xd) + P0) })
  cVbeta <- lapply(Vbeta, function(V) t(chol(V)))

  # state. Cut-points are anchored at kappa[d,1] = 0 (no separate intercept; the mean-zero random
  # effects carry the rest of the location), which removes the location/cut-point ridge and keeps the
  # sampler mixing; the remaining C-2 cut-points per route are free under a weakly-informative,
  # UNBIASED N(0, cut_sd) prior (the earlier N(1,1) prior biased the lower cut-points).
  beta  <- matrix(0, D, P)                     # route-specific fixed effects
  u     <- matrix(0, I, D); Sigma <- diag(D)
  a     <- rep(0, S); tau_s2 <- 0.25           # source effects and their variance
  kappa <- matrix(rep(0:(C - 2L), each = D), D, C - 1L)        # anchored start 0,1,2,...

  keep <- niter - nburn
  out_beta  <- array(NA, c(keep, D, P))
  out_kappa <- array(NA, c(keep, D, C - 1L))
  out_Sigma <- array(NA, c(keep, D, D))
  out_u     <- array(NA, c(keep, I, D))
  out_a     <- matrix(NA, keep, S)
  out_taus  <- numeric(keep)

  for (it in 1:niter) {
    Xbeta <- rowSums(X * beta[dimension, , drop = FALSE])     # X_j' beta_{d(j)}
    eta   <- Xbeta + u[cbind(species, dimension)] + a[source]

    # (1) latent z | rest, truncated to the y-category interval
    lo <- ifelse(y == 1, -Inf, kappa[cbind(dimension, pmax(y - 1L, 1L))])
    hi <- ifelse(y == C,  Inf, kappa[cbind(dimension, pmin(y, C - 1L))])
    z  <- rtn(eta, sigma, lo, hi)

    # (2) route-specific beta_d | rest
    for (d in 1:D) {
      R  <- rows_d[[d]]
      r  <- z[R] - u[species[R], d] - a[source[R]]
      m  <- Vbeta[[d]] %*% crossprod(X[R, , drop = FALSE] * w[R], r)
      beta[d, ] <- as.vector(m + cVbeta[[d]] %*% rnorm(P))
    }
    Xbeta <- rowSums(X * beta[dimension, , drop = FALSE])

    # (3) correlated species effects u_i | rest
    e_u   <- z - Xbeta - a[source]
    precD <- sapply(1:D, function(d) {
      rs <- rep(0, I); g <- rowsum(w[rows_d[[d]]], species[rows_d[[d]]])
      rs[as.integer(rownames(g))] <- g[, 1]; rs })
    bD <- sapply(1:D, function(d) {
      rs <- rep(0, I); g <- rowsum((e_u * w)[rows_d[[d]]], species[rows_d[[d]]])
      rs[as.integer(rownames(g))] <- g[, 1]; rs })
    Sinv <- solve(Sigma)
    for (i in 1:I) {
      V <- solve(Sinv + diag(precD[i, ], D))
      u[i, ] <- V %*% bD[i, ] + t(chol(V)) %*% rnorm(D)
    }

    # (4) Sigma | rest (inverse-Wishart)
    Sigma <- solve(rWishart(1, nu0 + I, solve(S0 + crossprod(u)))[, , 1])

    # (5) source random intercept a_s | rest (regularised toward modest source clustering)
    if (use_source) {
      e_a  <- z - Xbeta - u[cbind(species, dimension)]
      sw   <- rowsum(w, source);        prec_s <- rep(0, S); prec_s[as.integer(rownames(sw))] <- sw[, 1]
      swe  <- rowsum(e_a * w, source);  bb_s   <- rep(0, S); bb_s[as.integer(rownames(swe))]  <- swe[, 1]
      v_s  <- 1 / (prec_s + 1 / tau_s2)
      a    <- v_s * bb_s + sqrt(v_s) * rnorm(S)
      # (6) source variance tau_s^2 | rest (inverse-gamma, weakly informative)
      tau_s2 <- 1 / rgamma(1, as0 + S / 2, bs0 + 0.5 * sum(a^2))
    }

    # (7) free ordered cut-points kappa_d | rest (kappa[,1] fixed at 0; sample m = 2..C-1 ascending)
    for (d in 1:D) for (m in 2:(C - 1L)) {
      upper_nb <- if (m == C - 1L) Inf else kappa[d, m + 1L]
      zc  <- z[dimension == d & y == m]         # category m  (latent below the cut)
      zc1 <- z[dimension == d & y == (m + 1L)]  # category m+1 (latent above the cut)
      L <- max(c(kappa[d, m - 1L], zc)); U <- min(c(upper_nb, zc1))
      kappa[d, m] <- rtn(cut_mu, cut_sd, L, U)
    }

    if (it > nburn) {
      k <- it - nburn
      out_beta[k, , ]  <- beta
      out_kappa[k, , ] <- kappa
      out_Sigma[k, , ] <- Sigma
      out_u[k, , ]     <- u
      out_a[k, ]       <- a
      out_taus[k]      <- sqrt(tau_s2)
    }
  }
  list(beta = out_beta, kappa = out_kappa, Sigma = out_Sigma, u = out_u,
       a = out_a, tau_s = out_taus)
}

run_model <- function(dat, y, dimension, species, sigma, X, source = dat$source,
                      nchains = 4, niter = 6000, nburn = 2500, nu0 = 40, tau0 = 0.7,
                      prior_scale = dat$prior_scale, cut_sd = 5, use_source = FALSE) {
  lapply(1:nchains, function(ch)
    gibbs_chain(dat, y, dimension, species, sigma, X, source = source,
                niter = niter, nburn = nburn, seed = 100 + ch,
                nu0 = nu0, tau0 = tau0, prior_scale = prior_scale, cut_sd = cut_sd,
                use_source = use_source))
}

rho_of <- function(ch) ch$Sigma[, 1, 2] / sqrt(ch$Sigma[, 1, 1] * ch$Sigma[, 2, 2])

pool_draws <- function(chains, dat) {
  keep <- dim(chains[[1]]$beta)[1]; S <- keep * length(chains)
  D <- dat$D; P <- dat$P; I <- dat$I; Cm1 <- dat$C - 1L; Ns <- dat$S
  beta_d  <- array(0, c(S, D, P)); kappa_d <- array(0, c(S, D, Cm1))
  u_d     <- array(0, c(S, I, D)); a_d <- matrix(0, S, Ns)
  for (ch in seq_along(chains)) {
    idx <- ((ch - 1) * keep + 1):(ch * keep)
    beta_d[idx, , ]  <- chains[[ch]]$beta
    kappa_d[idx, , ] <- chains[[ch]]$kappa
    u_d[idx, , ]     <- chains[[ch]]$u
    a_d[idx, ]       <- chains[[ch]]$a
  }
  list(beta = beta_d, kappa = kappa_d, u = u_d, a = a_d, S = S)
}

# Per-species evidentiary reach at the species' reference design X_ref (West/kelp context, own
# traits), a population-typical source (a = 0), and High specificity (sigma = 1). Returns the expected
# rung (category - 1) on each route, with a 95% credible interval. Uses the ROUTE-SPECIFIC beta_d.
derive_ranks <- function(dat, draws) {
  beta_d <- draws$beta; kappa_d <- draws$kappa; u_d <- draws$u; C <- dat$C
  Xref <- dat$X_ref
  n_id <- table(factor(dat$species, levels = 1:dat$I), factor(dat$dimension, levels = 1:dat$D))
  est  <- data.frame(species = dat$species_lab, stringsAsFactors = FALSE)
  for (d in 1:dat$D) {
    ERm <- RCm <- RClo <- RChi <- numeric(dat$I)
    for (i in 1:dat$I) {
      loc <- as.vector(beta_d[, d, ] %*% Xref[i, ]) + u_d[, i, d]   # S draws, West context, route d
      kap <- matrix(kappa_d[, d, ], nrow = length(loc))            # S x (C-1)
      er  <- rep(0, length(loc)); prev <- rep(0, length(loc))
      for (cc in 1:C) {
        cur <- if (cc == C) rep(1, length(loc)) else pnorm(kap[, cc] - loc)
        er  <- er + cc * (cur - prev); prev <- cur
      }
      reach <- er - 1
      ERm[i] <- mean(er); RCm[i] <- mean(reach)
      RClo[i] <- quantile(reach, .025); RChi[i] <- quantile(reach, .975)
    }
    tag <- dat$dim_lab[d]
    est[[paste0("n_", tag)]]        <- as.integer(n_id[, d])
    est[[paste0(tag, "_erank")]]    <- round(ERm, 3)
    est[[paste0(tag, "_reach")]]    <- round(RCm, 3)
    est[[paste0(tag, "_reach_lo")]] <- round(RClo, 3)
    est[[paste0(tag, "_reach_hi")]] <- round(RChi, 3)
  }
  est$T_evidence <- est$n_T > 0
  est$H_evidence <- est$n_H > 0
  est[order(-est$T_reach), ]
}

# --------------------------------------------------------------------------------
# Convergence diagnostics: rank-normalised split-Rhat and bulk-ESS (Vehtari et al. 2021, Bayesian
# Analysis 16:667-718). M is an iterations x chains matrix of draws for one scalar quantity.
# --------------------------------------------------------------------------------
.rank_normalise <- function(x) {
  r <- rank(x, ties.method = "average")
  qnorm((r - 3/8) / (length(x) + 1/4))            # Blom rankits
}
.rhat_rank <- function(M) {                        # classic Rhat on a matrix iter x chain
  N <- nrow(M); K <- ncol(M)
  cm <- colMeans(M); gm <- mean(cm)
  B <- N / (K - 1) * sum((cm - gm)^2)
  W <- mean(apply(M, 2, var))
  if (W <= 0) return(NA_real_)
  sqrt(((N - 1) / N * W + B / N) / W)
}
.fft_acov <- function(x) {                         # biased autocovariance via FFT
  n <- length(x); x <- x - mean(x)
  m <- 2^ceiling(log2(2 * n))
  f <- fft(c(x, rep(0, m - n)))
  Re(fft(Mod(f)^2, inverse = TRUE))[1:n] / (m * n)
}
.ess <- function(M) {                              # bulk-ESS from rank-normalised draws
  N <- nrow(M); K <- ncol(M)
  if (N < 8) return(NA_real_)
  acov <- sapply(1:K, function(k) .fft_acov(M[, k]))   # N x K
  cvar <- apply(M, 2, var)                              # within-chain variance (n-1)
  W <- mean(cvar)
  B <- if (K > 1) N * var(colMeans(M)) else 0
  var_plus <- (N - 1) / N * W + (K > 1) * B / N
  if (var_plus <= 0) return(NA_real_)
  rho <- 1 - (W - rowMeans(acov)) / var_plus           # rho[1] = lag 0
  # Geyer initial positive + monotone on paired lags
  t <- 1; Gamma <- numeric(0)
  while (2 * t + 1 <= N) {
    g <- rho[2 * t] + rho[2 * t + 1]                   # lags (2t-1, 2t)
    if (g < 0) break
    Gamma <- c(Gamma, g); t <- t + 1
  }
  Gamma <- cummin(pmax(Gamma, 0))                      # enforce monotone non-increasing
  tau <- -1 + 2 * sum(Gamma)
  if (!is.finite(tau) || tau < 1) tau <- 1
  (N * K) / tau
}
mcmc_diag <- function(M) {                         # split-Rhat + bulk-ESS for one scalar
  N <- nrow(M); K <- ncol(M); half <- floor(N / 2)
  Msplit <- cbind(M[1:half, , drop = FALSE], M[(N - half + 1):N, , drop = FALSE])
  z <- matrix(.rank_normalise(as.vector(Msplit)), nrow = half)
  c(Rhat = .rhat_rank(z), ESS = round(.ess(z)))
}
# Backwards-compatible name used by the fit/sensitivity scripts.
rhat_ess <- function(M) mcmc_diag(M)

gate_convergence <- function(rhats, warn_at = 1.01, stop_at = 1.05) {
  if (anyNA(rhats)) stop("Convergence gate: one or more split-Rhat values are NA.")
  mx <- max(rhats)
  if (mx > stop_at) stop(sprintf("Convergence gate: max split-Rhat = %.3f exceeds %.2f; not converged.", mx, stop_at))
  if (mx > warn_at) warning(sprintf("Convergence gate: max split-Rhat = %.3f exceeds %.2f; mixing marginal.", mx, warn_at))
  invisible(mx)
}
