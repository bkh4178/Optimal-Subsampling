compute_pi <- function(score, k, max_iter = 100, tol = 1e-10) {
  
  N <- length(score)
  
  # 입력 검사
  if (any(!is.finite(score))) stop("score must be finite.")
  if (any(score < 0))         stop("score must be non-negative.")
  if (sum(score) == 0)        stop("sum(score) must be positive.")
  if (k <= 0 || k >= N)       stop("k must be in (0, N).")
  
  # 초기 정규화
  pi    <- score / sum(score) * k
  fixed <- rep(FALSE, N)
  
  for (iter in 1:max_iter) {
    
    # 이번 반복에서 새로 clip될 항목 (free set 중 pi > 1)
    new_clip <- !fixed & (pi > 1)
    
    if (!any(new_clip)) break
    
    # fixed set에 누적
    fixed[new_clip] <- TRUE
    pi[fixed] <- 1
    
    # 남은 budget을 free set에 score 비례 재분배
    budget_remain <- k - sum(fixed)
    free <- !fixed
    
    if (budget_remain <= 0 || !any(free)) {
      pi[free] <- 0
      break
    }
    
    pi[free] <- score[free] / sum(score[free]) * budget_remain
  }
  
  # 부동소수점 오차 방지
  pi <- pmax(pmin(pi, 1), 0)
  
  # sum 검증
  if (abs(sum(pi) - k) > tol * k) {
    warning(sprintf("sum(pi) = %.6f, target k = %.6f, diff = %.2e",
                    sum(pi), k, abs(sum(pi) - k)))
  }
  
  return(pi)
}


fit_wls <- function(X, y, idx, pi) {
  # idx: 선택된 관측치 인덱스
  # pi: 해당 관측치의 inclusion probability
  # weight = 1 / pi_i (IPW)
  
  X_sub <- X[idx, , drop = FALSE]
  y_sub <- y[idx]
  w     <- 1 / pi[idx]
  
  # weighted OLS: beta = (X^T W X)^{-1} X^T W y
  XtW  <- t(X_sub) * w          # p × n_sub
  XtWX <- XtW %*% X_sub         # p × p
  XtWy <- XtW %*% y_sub         # p × 1
  
  beta_hat <- solve(XtWX, XtWy)
  as.vector(beta_hat)
}


run_full <- function(dat) {
  N <- nrow(dat$X)
  idx <- 1:N
  pi  <- rep(1, N)
  
  beta_hat <- fit_wls(dat$X, dat$y, idx = idx, pi = pi)
  
  list(
    method     = "FULL",
    beta_hat   = beta_hat,
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
  N   <- nrow(dat$X)
  
  # score = ell_i = x_i^T A_hat^{-1} x_i
  pi  <- compute_pi(score = dat$ell, k = k)
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
  N   <- nrow(dat$X)
  
  # score = sqrt(ell_i)
  pi  <- compute_pi(score = sqrt(dat$ell), k = k)
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
  N   <- nrow(dat$X)
  
  # score = sqrt(sigma2_i * ell_i), 진짜 sigma2 사용
  pi  <- compute_pi(score = sqrt(dat$sigma2 * dat$ell), k = k)
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
  
  # 열 이름 통일 (lm 변수명 매칭 문제 방지)
  X_df <- as.data.frame(dat$X)
  colnames(X_df) <- paste0("V", 1:p)
  
  # Step 1. pilot sample (uniform, size n0)
  idx_pilot  <- sample(N, n0, replace = FALSE)
  pi_pilot   <- rep(n0 / N, N)
  beta_pilot <- fit_wls(dat$X, dat$y, idx = idx_pilot, pi = pi_pilot)
  
  # Step 2. pilot residual로 sigma2 추정
  resid_pilot <- as.vector(dat$y[idx_pilot] - dat$X[idx_pilot, ] %*% beta_pilot)
  log_resid2  <- log(resid_pilot^2 + 1e-10)
  
  var_model  <- lm(log_resid2 ~ ., data = X_df[idx_pilot, , drop = FALSE])
  sigma2_hat <- exp(predict(var_model, newdata = X_df))
  
  # NA/Inf 방어
  bad <- !is.finite(sigma2_hat)
  if (any(bad)) {
    sigma2_hat[bad] <- median(sigma2_hat[!bad])
  }
  
  # quantile clipping + 평균 1 정규화
  sigma2_hat <- pmin(sigma2_hat, quantile(sigma2_hat, 0.99))
  sigma2_hat <- pmax(sigma2_hat, 1e-6)
  sigma2_hat <- sigma2_hat / mean(sigma2_hat)
  
  # Step 3. plugin score
  pi  <- compute_pi(score = sqrt(sigma2_hat * dat$ell), k = k)
  idx <- which(rbinom(N, 1, pi) == 1)
  
  list(
    method                = "OPT-hetero-plugin",
    beta_hat              = fit_wls(dat$X, dat$y, idx = idx, pi = pi),
    pi                    = pi,
    idx                   = idx,
    n_realized            = length(idx),
    clip_rate             = mean(pi >= 1 - 1e-10),
    sigma2_hat            = sigma2_hat,
    n_pilot               = n0,
    idx_pilot             = idx_pilot,
    pilot_second_overlap  = length(intersect(idx_pilot, idx))
  )
}
