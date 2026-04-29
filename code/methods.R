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


set.seed(42)
N <- 10000

# 1. uniform score → pi_i = k/N
score_uni <- rep(1, N)
pi_uni <- compute_pi(score_uni, k = 500)
cat("uniform: sum =", sum(pi_uni), "| range =", range(pi_uni), "\n")
# 기대: sum = 500, all pi = 0.05

# 2. 일반적인 score (leverage-like)
score_lev <- rchisq(N, df = 12)
pi_lev <- compute_pi(score_lev, k = 500)
cat("leverage: sum =", round(sum(pi_lev), 6),
    "| range =", round(range(pi_lev), 4),
    "| clip rate =", round(mean(pi_lev >= 1 - 1e-10), 4), "\n")

# 3. 극단값 포함 score
score_ext <- c(rep(1, N - 5), rep(1e6, 5))
pi_ext <- compute_pi(score_ext, k = 500)
cat("extreme: sum =", round(sum(pi_ext), 6),
    "| clip rate =", round(mean(pi_ext >= 1 - 1e-10), 4), "\n")
# 기대: 극단값 5개 pi=1, 나머지에 495 재분배

# 4. k = 200, 500, 1000 전부 통과하는지
for (k_val in c(200, 500, 1000)) {
  pi_test <- compute_pi(score_lev, k = k_val)
  cat("k =", k_val, "| sum =", round(sum(pi_test), 6),
      "| clip rate =", round(mean(pi_test >= 1 - 1e-10), 4), "\n")
}


dat <- generate_data(case = 1)
N   <- nrow(dat$X)

# FULL: 일반 OLS와 동일해야 함
beta_wls  <- fit_wls(dat$X, dat$y, idx = 1:N, pi = rep(1, N))
beta_ols  <- as.vector(solve(crossprod(dat$X), crossprod(dat$X, dat$y)))
cat("max diff (FULL vs OLS):", max(abs(beta_wls - beta_ols)), "\n")
# 기대: 거의 0


dat <- generate_data(case = 1)

res_full <- run_full(dat)
res_uni  <- run_uni_ipw(dat, k = 500)

res_full$method
res_uni$method
res_uni$n_realized
mean(res_uni$pi)
sum(res_uni$pi)
res_uni$clip_rate


set.seed(1)
dat <- generate_data(case = 2)
res_plugin <- run_opt_hetero_plugin(dat, k = 500, n0 = 200)

cat("n_realized:", res_plugin$n_realized, "\n")
cat("clip_rate:", res_plugin$clip_rate, "\n")
cat("sum(pi):", round(sum(res_plugin$pi), 6), "\n")
cat("sigma2_hat summary:\n")
print(summary(res_plugin$sigma2_hat))
cat("cor(sigma2_hat, true sigma2):", 
    round(cor(res_plugin$sigma2_hat, dat$sigma2), 4), "\n")


cat("cor(sigma2_hat, ell):",
    round(cor(res_plugin$sigma2_hat, dat$ell), 4), "\n")

cat("pi range:\n")
print(range(res_plugin$pi))

cat("pi summary:\n")
print(summary(res_plugin$pi))

cat("true sigma2 summary:\n")
print(summary(dat$sigma2))

for (cc in c(1, 2, 3)) {
  set.seed(1)
  dat <- generate_data(case = cc)
  res <- run_opt_hetero_plugin(dat, k = 500, n0 = 200)
  
  cat("\nCase", cc, "\n")
  cat("n_realized:", res$n_realized, "\n")
  cat("sum(pi):", round(sum(res$pi), 6), "\n")
  cat("clip_rate:", res$clip_rate, "\n")
  cat("sigma2_hat summary:\n")
  print(summary(res$sigma2_hat))
  cat("cor(sigma2_hat, true sigma2):",
      round(cor(res$sigma2_hat, dat$sigma2), 4), "\n")
  cat("cor(sigma2_hat, ell):",
      round(cor(res$sigma2_hat, dat$ell), 4), "\n")
}