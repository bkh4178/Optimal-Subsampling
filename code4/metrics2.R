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

# 4. Theoretical leading term C1(pi) = sum(sigma2_i * ell_i / pi_i)
# 이론에서 최소화하는 목적함수의 sample analogue
# 주의: FULL은 pi=1이라 C1 = sum(sigma2*ell), 다른 방법들과 단위 다름 → 비교 시 FULL 제외
compute_c1 <- function(pi, sigma2, ell) {
  sum(sigma2 * ell / pi)
}

# 5. A_hat_pi^{-1} 근사 오차 및 condition number
# 이론: (X_s' W X_s / N)^{-1} →p A^{-1} as k,N → ∞
# finite sample에서 이 근사가 얼마나 불안정한지 측정
compute_Ainv_error <- function(dat, idx, pi) {
  X_s      <- dat$X[idx, , drop = FALSE]
  w        <- 1 / pi[idx]
  A_hat_pi <- crossprod(X_s * sqrt(w)) / nrow(dat$X)

  out <- tryCatch({
    A_hat_pi_inv <- solve(A_hat_pi)
    list(
      Ainv_frob_error = norm(A_hat_pi_inv - dat$A_hat_inv, type = "F"),
      A_pi_cond       = kappa(A_hat_pi),
      A_solve_fail    = FALSE
    )
  }, error = function(e) {
    list(
      Ainv_frob_error = NA_real_,
      A_pi_cond       = kappa(A_hat_pi),
      A_solve_fail    = TRUE
    )
  })
  out
}

# 6. 한 번에 모든 지표 계산 → data.frame 반환
compute_metrics <- function(result, dat, X_test, y_test) {
  ainv <- compute_Ainv_error(dat, result$idx, result$pi)

  data.frame(
    method          = result$method,
    param_mse       = compute_param_mse(result$beta_hat, dat$beta0),
    test_mse        = compute_test_mse(result$beta_hat, X_test, y_test),
    excess_risk     = compute_excess_risk(result$beta_hat, dat$beta0, X_test),
    c1_pi           = compute_c1(result$pi, dat$sigma2, dat$ell),
    Ainv_frob_error = ainv$Ainv_frob_error,
    A_pi_cond       = ainv$A_pi_cond,
    A_solve_fail    = ainv$A_solve_fail,
    n_realized      = result$n_realized,
    clip_rate       = result$clip_rate
  )
}