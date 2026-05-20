# metrics.R
# 지표 계산 함수 — Phase 1
#
# 함수:
#   compute_param_mse   : parameter MSE = ||beta_hat - beta0||^2
#   compute_excess_risk : excess prediction risk = mean((X_test %*% (beta0 - beta_hat))^2)
#   compute_test_mse    : test MSE = mean((y_test - X_test %*% beta_hat)^2)
#   compute_c1          : 이론적 목적함수 sample analogue = (1/N) * sum(sigma2 * ell / pi)
#   compute_metrics     : 위 네 개를 한 번에 계산 → data.frame 한 행 반환

# 1. Parameter MSE
compute_param_mse <- function(beta_hat, beta0) {
  sum((beta_hat - beta0)^2)
}

# 2. Excess prediction risk
# test point x*가 train과 같은 모집단에서 온다고 가정
# E[(x*'beta0 - x*'beta_hat)^2] ≈ mean((X_test %*% (beta0 - beta_hat))^2)
compute_excess_risk <- function(beta_hat, beta0, X_test) {
  diff <- X_test %*% (beta0 - beta_hat)
  mean(diff^2)
}

# 3. Test MSE (noisy response 기준)
# y_test에 irreducible noise epsilon이 포함 → excess risk보다 항상 크거나 같음
compute_test_mse <- function(beta_hat, X_test, y_test) {
  resid <- y_test - X_test %*% beta_hat
  mean(resid^2)
}

# 4. 이론적 목적함수 sample analogue
# leading term of excess risk: (1/N) * sum(sigma2_i * ell_i / pi_i)
# 주의: FULL은 pi_i = 1이라 다른 방법들과 단위가 달라 → c1 비교 시 FULL 제외
compute_c1 <- function(pi, sigma2, ell, N) {
  sum(sigma2 * ell / pi) / N
}

# 5. 한 번에 모든 지표 계산 → data.frame 한 행 반환
compute_metrics <- function(result, dat, X_test, y_test) {
  N <- nrow(dat$X)

  data.frame(
    method      = result$method,
    param_mse   = compute_param_mse(result$beta_hat, dat$beta0),
    excess_risk = compute_excess_risk(result$beta_hat, dat$beta0, X_test),
    test_mse    = compute_test_mse(result$beta_hat, X_test, y_test),
    c1_pi       = compute_c1(result$pi, dat$sigma2, dat$ell, N),
    n_realized  = result$n_realized,
    clip_rate   = result$clip_rate
  )
}