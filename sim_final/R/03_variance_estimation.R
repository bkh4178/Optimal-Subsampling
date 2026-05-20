# 03_variance_estimation.R
# sigma^2(Z) plugin 추정
#
# 참 분산 모형: sigma(Z) = 0.5 + |0.5*Z1 + 0.5*Z2 + Z6 + Z8|
# Plugin 모형: abs_resid ~ abs_s  (sigma-scale 선형 근사)
#              sigma2_hat = (sigma_hat / sqrt(2/pi))^2
#              abs_s  = |0.5*Z1 + 0.5*Z2 + Z6 + Z8|
#
# X 열 순서: (intercept, Z1, Z2, Z3, Z4, Z5, Z6, Z7, Z8, Z9)
#   col:          1       2   3   4   5   6   7   8   9   10

make_abs_s <- function(X) {
  Z1 <- X[, 2L]
  Z2 <- X[, 3L]
  Z6 <- X[, 7L]
  Z8 <- X[, 9L]
  abs(0.5 * Z1 + 0.5 * Z2 + Z6 + Z8)
}

# 파라미터:
#   dat        : generate_data_final() 반환값
#   idx_pilot  : pilot 인덱스 (정수 벡터)
#   beta_pilot : pilot에서 추정한 beta (길이 10)
# 반환: 길이 N 의 양수 벡터 (sigma^2 hat)
estimate_sigma2_plugin <- function(dat, idx_pilot, beta_pilot) {
  resid_pilot     <- dat$y[idx_pilot] -
                     as.numeric(dat$X[idx_pilot, , drop = FALSE] %*% beta_pilot)
  abs_resid_pilot <- abs(resid_pilot)

  abs_s_pilot <- make_abs_s(dat$X[idx_pilot, , drop = FALSE])
  pilot_df    <- data.frame(abs_resid = abs_resid_pilot, abs_s = abs_s_pilot)
  sig_model   <- lm(abs_resid ~ abs_s, data = pilot_df)

  pred_df   <- data.frame(abs_s = make_abs_s(dat$X))
  sigma_hat <- as.numeric(predict(sig_model, newdata = pred_df))

  sigma2_hat <- (sigma_hat / sqrt(2 / pi))^2

  bad <- !is.finite(sigma2_hat) | sigma2_hat <= 0
  if (any(bad)) {
    good_vals       <- sigma2_hat[!bad]
    fallback        <- if (length(good_vals) > 0L) median(good_vals) else 1
    sigma2_hat[bad] <- fallback
  }
  pmax(sigma2_hat, 1e-10)
}
