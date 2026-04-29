# metrics.R
# 지표 계산 함수들
#
# 입력: method 결과 (run_* 반환값), dat (generate_data 반환값), X_test / y_test
# 출력: data.frame (한 행)

# 1. Parameter MSE
compute_param_mse <- function(beta_hat, beta0) {
  sum((beta_hat - beta0)^2)
}

# 2. Test MSE
compute_test_mse <- function(beta_hat, X_test, y_test) {
  resid <- y_test - X_test %*% beta_hat
  mean(resid^2)
}

# 3. Excess risk
compute_excess_risk <- function(beta_hat, beta0, X_test) {
  diff <- X_test %*% (beta0 - beta_hat)
  mean(diff^2)
}

# 4. 한 번에 모든 지표 계산 → data.frame 반환
compute_metrics <- function(result, dat, X_test, y_test) {
  data.frame(
    method      = result$method,
    param_mse   = compute_param_mse(result$beta_hat, dat$beta0),
    test_mse    = compute_test_mse(result$beta_hat, X_test, y_test),
    excess_risk = compute_excess_risk(result$beta_hat, dat$beta0, X_test),
    n_realized  = result$n_realized,
    clip_rate   = result$clip_rate
  )
}

set.seed(1)
dat    <- generate_data(case = 1)
dat_te <- generate_data(case = 1)

res_full   <- run_full(dat)
res_oracle <- run_opt_hetero_oracle(dat, k = 500)

m_full   <- compute_metrics(res_full,   dat, dat_te$X, dat_te$y)
m_oracle <- compute_metrics(res_oracle, dat, dat_te$X, dat_te$y)

# FULL이 가장 낮은 excess risk, test MSE여야 함
print(m_full)
print(m_oracle)