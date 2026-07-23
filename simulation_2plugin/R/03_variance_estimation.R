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
estimate_sigma2_plugin_param <- function(dat, idx_pilot, beta_pilot) {
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


# TODO(상준이 형 코드 반영, 07/21): Random Forest 기반 nonparametric plugin.
# param과 동일하게 beta_pilot을 공유받아 잔차를 계산 (pilot lm 재적합 안 함 —
# param/nonparam 비교가 순수하게 "분산 추정 방식" 차이만 반영하도록 하기 위함).
# 예측변수는 Z1~Z9 전체 (intercept 제외).

estimate_sigma2_plugin_nonparam <- function(dat, idx_pilot, beta_pilot,
                                             ntree    = 100,
                                             mtry     = RF_MTRY_FIXED,
                                             nodesize = RF_NODESIZE_FIXED) {
  X_pilot <- dat$X[idx_pilot, , drop = FALSE]

  resid_pilot     <- dat$y[idx_pilot] - as.numeric(X_pilot %*% beta_pilot)
  abs_resid_pilot <- abs(resid_pilot)

  Z_pilot  <- as.data.frame(X_pilot[, -1, drop = FALSE])
  pilot_df <- data.frame(abs_resid = abs_resid_pilot, Z_pilot)

  rf_model <- randomForest::randomForest(abs_resid ~ ., data = pilot_df,
                                          ntree = ntree, mtry = mtry, nodesize = nodesize)

  Z_full     <- as.data.frame(dat$X[, -1, drop = FALSE])
  pred_abs   <- as.numeric(predict(rf_model, newdata = Z_full))

  sigma_hat  <- pmax(pred_abs, 1e-10) / sqrt(2 / pi)   # param과 동일한 정규분포 보정
  pred_var_full <- sigma_hat^2

  bad <- !is.finite(pred_var_full) | pred_var_full <= 0
  if (any(bad)) {
    good_vals <- pred_var_full[!bad]
    pred_var_full[bad] <- if (length(good_vals) > 0L) median(good_vals) else 1
  }
  pmax(pred_var_full, 1e-10)
}