# methods.R
# Subsampling methods — Phase 1 (저차원 검증)
#
# Methods:
#   run_full                 : full data OLS
#   run_uni_ipw              : uniform inclusion probability
#   run_lev_ipw              : leverage-proportional inclusion probability
#   run_opt_homo             : OPT under homoscedasticity  (pi ∝ sqrt(ell))
#   run_opt_hetero_oracle    : OPT-hetero with true sigma2 (pi ∝ sqrt(sigma2 * ell))
#   run_opt_hetero_plugin    : OPT-hetero with pilot-estimated sigma2
#                              pilot: uniform n0 draw → lm(log(resid^2) ~ X) → exp(predict)
#
# 공통:
#   - Bernoulli subsampling: delta_i ~ Bernoulli(pi_i)
#   - sum(pi) = k (expected subsample size)
#   - weight = 1/pi_i (IPW)

# ── 유틸리티 ──────────────────────────────────────────────────────────────────

compute_pi <- function(score, k, max_iter = 100, tol = 1e-10) {
  # score 벡터를 받아 sum(pi) = k가 되도록 정규화
  # pi_i > 1인 관측치는 1로 clip하고 나머지 budget을 재분배 (iterative)

  N <- length(score)

  if (any(!is.finite(score))) stop("score must be finite.")
  if (any(score < 0))         stop("score must be non-negative.")
  if (sum(score) == 0)        stop("sum(score) must be positive.")
  if (k <= 0 || k >= N)       stop("k must be in (0, N).")

  pi    <- score / sum(score) * k
  fixed <- rep(FALSE, N)

  for (iter in seq_len(max_iter)) {
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

  pmax(pmin(pi, 1), 0)
}


fit_wls <- function(X, y, idx, pi) {
  # 선택된 idx의 subsample로 weighted OLS
  # beta_hat = (X'WX)^{-1} X'Wy,  W = diag(1/pi_i)

  X_sub <- X[idx, , drop = FALSE]
  y_sub <- y[idx]
  w     <- 1 / pi[idx]

  XtW  <- t(X_sub) * w
  XtWX <- XtW %*% X_sub
  XtWy <- XtW %*% y_sub

  as.vector(solve(XtWX, XtWy))
}

# ── Methods ───────────────────────────────────────────────────────────────────

run_full <- function(dat) {
  N   <- nrow(dat$X)
  idx <- seq_len(N)
  pi  <- rep(1, N)

  list(
    method     = "FULL",
    beta_hat   = fit_wls(dat$X, dat$y, idx, pi),
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
    beta_hat   = fit_wls(dat$X, dat$y, idx, pi),
    pi         = pi,
    idx        = idx,
    n_realized = length(idx),
    clip_rate  = 0
  )
}


run_lev_ipw <- function(dat, k) {
  # pi ∝ ell_i (leverage에 비례)
  N  <- nrow(dat$X)
  pi <- compute_pi(score = dat$ell, k = k)
  idx <- which(rbinom(N, 1, pi) == 1)

  list(
    method     = "LEV-IPW",
    beta_hat   = fit_wls(dat$X, dat$y, idx, pi),
    pi         = pi,
    idx        = idx,
    n_realized = length(idx),
    clip_rate  = mean(pi >= 1 - 1e-10)
  )
}


run_opt_homo <- function(dat, k) {
  # pi ∝ sqrt(ell_i)
  # Cauchy-Schwarz로 유도된 homoscedastic 하의 최적 샘플링 확률
  N  <- nrow(dat$X)
  pi <- compute_pi(score = sqrt(dat$ell), k = k)
  idx <- which(rbinom(N, 1, pi) == 1)

  list(
    method     = "OPT-homo",
    beta_hat   = fit_wls(dat$X, dat$y, idx, pi),
    pi         = pi,
    idx        = idx,
    n_realized = length(idx),
    clip_rate  = mean(pi >= 1 - 1e-10)
  )
}


run_opt_hetero_oracle <- function(dat, k) {
  # pi ∝ sqrt(sigma2_i * ell_i)
  # 진짜 sigma2 사용 → asymptotic optimal (finite sample에서 불안정 가능)
  N  <- nrow(dat$X)
  pi <- compute_pi(score = sqrt(dat$sigma2 * dat$ell), k = k)
  idx <- which(rbinom(N, 1, pi) == 1)

  list(
    method     = "OPT-hetero-oracle",
    beta_hat   = fit_wls(dat$X, dat$y, idx, pi),
    pi         = pi,
    idx        = idx,
    n_realized = length(idx),
    clip_rate  = mean(pi >= 1 - 1e-10)
  )
}


run_opt_hetero_plugin <- function(dat, k, n0 = 200) {
  # pi ∝ sqrt(sigma2_hat * ell_i)
  # sigma2 추정:
  #   1. uniform pilot draw (n0개)
  #   2. pilot residual 계산
  #   3. lm(log(resid^2 + eps) ~ X) 적합 (1차 선형)
  #   4. 전체 N개에 predict → exp(·)로 복원
  #   5. 평균 1로 정규화

  N <- nrow(dat$X)
  p <- ncol(dat$X)

  # Step 1. pilot
  idx_pilot  <- sample(N, n0, replace = FALSE)
  pi_pilot   <- rep(n0 / N, N)
  beta_pilot <- fit_wls(dat$X, dat$y, idx_pilot, pi_pilot)

  # Step 2. log residual squared model
  resid      <- dat$y[idx_pilot] - dat$X[idx_pilot, , drop = FALSE] %*% beta_pilot
  log_resid2 <- log(as.vector(resid)^2 + 1e-10)

  X_pilot_df <- as.data.frame(dat$X[idx_pilot, , drop = FALSE])
  colnames(X_pilot_df) <- paste0("V", seq_len(p))

  var_model <- lm(log_resid2 ~ ., data = X_pilot_df)

  # Step 3. 전체 N에 predict
  X_full_df <- as.data.frame(dat$X)
  colnames(X_full_df) <- paste0("V", seq_len(p))

  sigma2_hat <- exp(predict(var_model, newdata = X_full_df))

  # NA/Inf 방어
  bad <- !is.finite(sigma2_hat)
  if (any(bad)) sigma2_hat[bad] <- median(sigma2_hat[!bad])

  # Step 4. 정규화
  sigma2_hat <- sigma2_hat / mean(sigma2_hat)

  # Step 5. plugin score → pi
  pi  <- compute_pi(score = sqrt(sigma2_hat * dat$ell), k = k)
  idx <- which(rbinom(N, 1, pi) == 1)

  list(
    method     = "OPT-hetero-plugin",
    beta_hat   = fit_wls(dat$X, dat$y, idx, pi),
    pi         = pi,
    idx        = idx,
    n_realized = length(idx),
    clip_rate  = mean(pi >= 1 - 1e-10),
    sigma2_hat = sigma2_hat
  )
}