# methods.R
# Subsampling methods for optimal subsampling simulation
#
# Methods:
#   run_full                    : full data OLS
#   run_uni_ipw                 : uniform inclusion probability
#   run_lev_ipw                 : leverage-proportional inclusion probability
#   run_opt_homo                : OPT under homoscedasticity
#   run_opt_hetero_oracle       : OPT-hetero with true sigma2
#   run_opt_hetero_plugin       : OPT-hetero with pilot-estimated sigma2 (raw)
#   run_opt_hetero_tempering    : OPT-hetero with score tempering (pi ∝ score^gamma)
#   run_opt_hetero_plugin_stab  : OPT-hetero with stabilized pilot-estimated sigma2
#                                 (normalize → clip → shrink → tempering)
#
# Note: pi_i = inclusion probability, sum(pi) ≈ k (Bernoulli subsampling)

# ── 공통 유틸리티 ────────────────────────────────────────────────────────────

compute_pi <- function(score, k, max_iter = 100, tol = 1e-10) {

  N <- length(score)

  if (any(!is.finite(score))) stop("score must be finite.")
  if (any(score < 0))         stop("score must be non-negative.")
  if (sum(score) == 0)        stop("sum(score) must be positive.")
  if (k <= 0 || k >= N)       stop("k must be in (0, N).")

  pi    <- score / sum(score) * k
  fixed <- rep(FALSE, N)

  for (iter in 1:max_iter) {
    new_clip <- !fixed & (pi > 1)
    if (!any(new_clip)) break

    fixed[new_clip] <- TRUE
    pi[fixed] <- 1

    budget_remain <- k - sum(fixed)
    free <- !fixed

    if (budget_remain <= 0 || !any(free)) {
      pi[free] <- 0
      break
    }

    pi[free] <- score[free] / sum(score[free]) * budget_remain
  }

  pi <- pmax(pmin(pi, 1), 0)

  if (abs(sum(pi) - k) > tol * k) {
    warning(sprintf("sum(pi) = %.6f, target k = %.6f, diff = %.2e",
                    sum(pi), k, abs(sum(pi) - k)))
  }

  return(pi)
}


fit_wls <- function(X, y, idx, pi) {
  X_sub <- X[idx, , drop = FALSE]
  y_sub <- y[idx]
  w     <- 1 / pi[idx]

  XtW  <- t(X_sub) * w
  XtWX <- XtW %*% X_sub
  XtWy <- XtW %*% y_sub

  beta_hat <- solve(XtWX, XtWy)
  as.vector(beta_hat)
}

# ── Methods ──────────────────────────────────────────────────────────────────

run_full <- function(dat) {
  N   <- nrow(dat$X)
  idx <- 1:N
  pi  <- rep(1, N)

  list(
    method     = "FULL",
    beta_hat   = fit_wls(dat$X, dat$y, idx = idx, pi = pi),
    pi         = pi,
    idx        = idx,
    n_realized = N,
    clip_rate  = 0
  )
}


run_uni_ipw <- function(dat, k) {
  N  <- nrow(dat$X)
  pi <- rep(k / N, N)
  idx <- which(rbinom(N, 1, pi) == 1)

  list(
    method     = "UNI-IPW",
    beta_hat   = fit_wls(dat$X, dat$y, idx = idx, pi = pi),
    pi         = pi,
    idx        = idx,
    n_realized = length(idx),
    clip_rate  = 0
  )
}


run_lev_ipw <- function(dat, k) {
  N  <- nrow(dat$X)
  pi <- compute_pi(score = dat$ell, k = k)
  idx <- which(rbinom(N, 1, pi) == 1)

  list(
    method     = "LEV-IPW",
    beta_hat   = fit_wls(dat$X, dat$y, idx = idx, pi = pi),
    pi         = pi,
    idx        = idx,
    n_realized = length(idx),
    clip_rate  = mean(pi >= 1 - 1e-10)
  )
}


run_opt_homo <- function(dat, k) {
  N  <- nrow(dat$X)
  pi <- compute_pi(score = sqrt(dat$ell), k = k)
  idx <- which(rbinom(N, 1, pi) == 1)

  list(
    method     = "OPT-homo",
    beta_hat   = fit_wls(dat$X, dat$y, idx = idx, pi = pi),
    pi         = pi,
    idx        = idx,
    n_realized = length(idx),
    clip_rate  = mean(pi >= 1 - 1e-10)
  )
}


run_opt_hetero_oracle <- function(dat, k) {
  N  <- nrow(dat$X)
  pi <- compute_pi(score = sqrt(dat$sigma2 * dat$ell), k = k)
  idx <- which(rbinom(N, 1, pi) == 1)

  list(
    method     = "OPT-hetero-oracle",
    beta_hat   = fit_wls(dat$X, dat$y, idx = idx, pi = pi),
    pi         = pi,
    idx        = idx,
    n_realized = length(idx),
    clip_rate  = mean(pi >= 1 - 1e-10)
  )
}


run_opt_hetero_plugin <- function(dat, k, n0 = 200) {
  N <- nrow(dat$X)
  p <- ncol(dat$X)

  X_df <- as.data.frame(dat$X)
  colnames(X_df) <- paste0("V", 1:p)

  # Step 1. pilot sample
  idx_pilot  <- sample(N, n0, replace = FALSE)
  pi_pilot   <- rep(n0 / N, N)
  beta_pilot <- fit_wls(dat$X, dat$y, idx = idx_pilot, pi = pi_pilot)

  # Step 2. sigma2 추정 (raw)
  resid_pilot <- as.vector(dat$y[idx_pilot] - dat$X[idx_pilot, ] %*% beta_pilot)
  log_resid2  <- log(resid_pilot^2 + 1e-10)

  var_model  <- lm(log_resid2 ~ ., data = X_df[idx_pilot, , drop = FALSE])
  sigma2_hat <- exp(predict(var_model, newdata = X_df))

  bad <- !is.finite(sigma2_hat)
  if (any(bad)) sigma2_hat[bad] <- median(sigma2_hat[!bad])

  sigma2_hat <- pmin(sigma2_hat, quantile(sigma2_hat, 0.99))
  sigma2_hat <- pmax(sigma2_hat, 1e-6)
  sigma2_hat <- sigma2_hat / mean(sigma2_hat)

  # Step 3. plugin score
  pi  <- compute_pi(score = sqrt(sigma2_hat * dat$ell), k = k)
  idx <- which(rbinom(N, 1, pi) == 1)

  list(
    method     = "OPT-hetero-plugin",
    beta_hat   = fit_wls(dat$X, dat$y, idx = idx, pi = pi),
    pi         = pi,
    idx        = idx,
    n_realized = length(idx),
    clip_rate  = mean(pi >= 1 - 1e-10),
    sigma2_hat = sigma2_hat,
    n_pilot    = n0
  )
}


run_opt_hetero_tempering <- function(dat, k, gamma = 0.5) {
  # Tempering: pi ∝ score^gamma, 0 < gamma < 1
  # gamma = 1 → oracle, gamma = 0 → uniform (연속적 보간)
  N  <- nrow(dat$X)

  score_raw <- sqrt(dat$sigma2 * dat$ell)
  pi <- compute_pi(score = score_raw^gamma, k = k)
  idx <- which(rbinom(N, 1, pi) == 1)

  list(
    method     = paste0("OPT-hetero-tempering-", gamma),
    beta_hat   = fit_wls(dat$X, dat$y, idx = idx, pi = pi),
    pi         = pi,
    idx        = idx,
    n_realized = length(idx),
    clip_rate  = mean(pi >= 1 - 1e-10),
    gamma      = gamma
  )
}


run_opt_hetero_plugin_stab <- function(dat, k, n0 = 200,
                                        lambda = 0.5,
                                        gamma  = 0.5,
                                        clip_L = 0.25,
                                        clip_U = 4.0) {
  # Stabilized plug-in: normalize → clip → shrink → score tempering
  # lambda: shrinkage toward homoskedastic baseline (0 = no shrink, 1 = full shrink)
  # gamma:  score tempering (1 = raw plugin score, 0 = uniform)
  # clip_L, clip_U: clipping bounds after normalization
  N <- nrow(dat$X)
  p <- ncol(dat$X)

  X_df <- as.data.frame(dat$X)
  colnames(X_df) <- paste0("V", 1:p)

  # X + X^2 feature matrix (nonlinear variance signal 포착)
  X_df_q <- X_df
  for (j in 1:p) X_df_q[[paste0("V", j, "sq")]] <- X_df[[paste0("V", j)]]^2

  # Step 1. Pilot → beta_hat
  idx_pilot  <- sample(N, n0, replace = FALSE)
  pi_pilot   <- rep(n0 / N, N)
  beta_pilot <- fit_wls(dat$X, dat$y, idx_pilot, pi_pilot)

  # Step 2. Residual
  resid <- as.vector(dat$y[idx_pilot] - dat$X[idx_pilot, ] %*% beta_pilot)

  # Step 3. Log residual squared model (X + X^2)
  # centered Bernoulli linear term이 S3 group effect 포착의 핵심
  c_eps  <- max(1e-8, 0.1 * median(resid^2))
  log_r2 <- log(resid^2 + c_eps)

  var_mod    <- lm(log_r2 ~ ., data = X_df_q[idx_pilot, , drop = FALSE])
  sigma2_hat <- exp(predict(var_mod, newdata = X_df_q))

  # NA/Inf 방어
  bad <- !is.finite(sigma2_hat)
  if (any(bad)) sigma2_hat[bad] <- median(sigma2_hat[!bad])

  # Step 4. Normalize (평균 1)
  sigma2_hat <- sigma2_hat / mean(sigma2_hat)

  # Step 5. Clip
  sigma2_hat <- pmin(pmax(sigma2_hat, clip_L), clip_U)

  # Step 6. Shrink toward 1 (homoskedastic baseline)
  sigma2_stab <- (1 - lambda) * sigma2_hat + lambda

  # Step 7. Score tempering
  score <- sqrt(sigma2_stab * dat$ell)
  pi    <- compute_pi(score^gamma, k = k)

  idx <- which(rbinom(N, 1, pi) == 1)

  list(
    method     = sprintf("OPT-plugin-stab(λ=%.1f,γ=%.1f)", lambda, gamma),
    beta_hat   = fit_wls(dat$X, dat$y, idx, pi),
    pi         = pi,
    idx        = idx,
    n_realized = length(idx),
    clip_rate  = mean(pi >= 1 - 1e-10),
    sigma2_hat = sigma2_stab,
    lambda     = lambda,
    gamma      = gamma
  )
}