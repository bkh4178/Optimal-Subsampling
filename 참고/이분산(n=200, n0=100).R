# ============================================================
# 이분산 선형회귀 시뮬레이션
# 비교 방법: OPT vs SRS vs FULL
# 최종 표본의 기대 크기는 200명이다.
# ============================================================

# ------------------------------------------------------------
# 1. 환경 설정 및 시드 고정
# ------------------------------------------------------------

set.seed(123)

N <- 10000
p <- 6


# ------------------------------------------------------------
# 2. 훈련용 모집단 생성
# ------------------------------------------------------------

x1 <- rnorm(N, 0, 1)
x2 <- rnorm(N, 0, 1)

x3 <- sample(
  1:5,
  N,
  replace = TRUE
)

x4 <- sample(
  1:5,
  N,
  replace = TRUE
)

x5 <- rbinom(
  N,
  1,
  0.5
)

x6 <- rbinom(
  N,
  1,
  0.1
)


# 절편을 포함한 설계행렬
X_pop <- cbind(
  intercept = 1,
  x1 = x1,
  x2 = x2,
  x3 = x3,
  x4 = x4,
  x5 = x5,
  x6 = x6
)


# 참 회귀계수
beta_true <- c(
  intercept = 1.0,
  x1 = 1.5,
  x2 = 1.0,
  x3 = 0.5,
  x4 = 0.1,
  x5 = -1.0,
  x6 = 2.0
)


# 이분산 적용
#
# sigma_vec은 조건부 분산이 아니라
# 조건부 표준편차 sigma(x)를 나타낸다.
#
#   sigma(x) = |x1 + x2| + 0.5

sigma_vec <- (
  abs(X_pop[, "x1"] + X_pop[, "x2"]) +
    0.5
)


# 이분산 오차항
error_hetero <- rnorm(
  N,
  mean = 0,
  sd = sigma_vec
)


# 반응변수
y_hetero <- drop(
  X_pop %*% beta_true +
    error_hetero
)


# 훈련용 데이터프레임
full_data_hetero <- data.frame(
  y = y_hetero,
  X_pop
)


# ------------------------------------------------------------
# 3. 독립 테스트 자료 생성
# ------------------------------------------------------------

N_test <- 10000


x1_t <- rnorm(N_test, 0, 1)
x2_t <- rnorm(N_test, 0, 1)

x3_t <- sample(
  1:5,
  N_test,
  replace = TRUE
)

x4_t <- sample(
  1:5,
  N_test,
  replace = TRUE
)

x5_t <- rbinom(
  N_test,
  1,
  0.5
)

x6_t <- rbinom(
  N_test,
  1,
  0.1
)


# 테스트 설계행렬
X_test <- cbind(
  intercept = 1,
  x1 = x1_t,
  x2 = x2_t,
  x3 = x3_t,
  x4 = x4_t,
  x5 = x5_t,
  x6 = x6_t
)


# 테스트 자료에도 동일한 이분산 구조 적용
sigma_vec_test <- (
  abs(X_test[, "x1"] + X_test[, "x2"]) +
    0.5
)


error_test <- rnorm(
  N_test,
  mean = 0,
  sd = sigma_vec_test
)


y_test <- drop(
  X_test %*% beta_true +
    error_test
)


test_data_hetero <- data.frame(
  y = y_test,
  X_test
)


# ------------------------------------------------------------
# 4. 레버리지 계산을 위한 A_hat
# ------------------------------------------------------------

# A_hat = X^T X / N
A_hat <- crossprod(X_pop) / N


# A_hat의 역행렬
A_inv <- solve(A_hat)


# 각 관측치의 leverage 부분
#
#   x_i^T A_hat^{-1} x_i

leverage_part <- rowSums(
  (X_pop %*% A_inv) * X_pop
)

# ------------------------------------------------------------
# 5. 시뮬레이션 설정
# ------------------------------------------------------------

# 파일럿 표본크기와 최종 표본 기대 크기는 별도이다.
#
# 파일럿 표본: 정확히 100명
# 최종 OPT/SRS 표본: 기대 크기 200명

n_total <- 200
n_pilot <- 100
n_reps <- 1000



X_test_mat <- as.matrix(
  test_data_hetero[, -1]
)

y_test_vec <- test_data_hetero$y


library(randomForest)


# ------------------------------------------------------------
# 6. 반복 시뮬레이션
# ------------------------------------------------------------

for (r in seq_len(n_reps)) {
  
  
  # ----------------------------------------------------------
  # Step 2. Pilot sampling
  # ----------------------------------------------------------
  
  idx_pilot <- sample(
    x = seq_len(N),
    size = n_pilot,
    replace = FALSE
  )
  
  
  d_pilot <- full_data_hetero[
    idx_pilot,
  ]
  
  
  # 파일럿 회귀모형
  #
  # intercept 열이 이미 데이터에 있으므로
  # R의 자동 절편은 제거한다.
  
  fit_pilot <- lm(
    y ~ . - 1,
    data = d_pilot
  )
  
  
  # ----------------------------------------------------------
  # Step 3. Random Forest를 이용한 조건부 분산 추정
  # ----------------------------------------------------------
  
  # 파일럿 표본의 잔차제곱
  pilot_res_sq <- (
    d_pilot$y -
      predict(
        fit_pilot,
        newdata = d_pilot
      )
  )^2
  
  
  # Random Forest 학습자료
  #
  # d_pilot의 첫 번째 열은 y,
  # 두 번째 열은 intercept이므로 둘 다 제외한다.
  
  d_pilot_var <- data.frame(
    res_sq = pilot_res_sq,
    d_pilot[, c(-1, -2)]
  )
  
  
  # 잔차제곱을 이용하여 조건부 분산함수 추정
  var_rf_model <- randomForest(
    res_sq ~ .,
    data = d_pilot_var,
    ntree = 100
  )
  
  
  # 모집단 전체에 대해 조건부 분산 예측
  pred_var_full <- predict(
    var_rf_model,
    newdata = full_data_hetero[, c(-1, -2)]
  )
  
  
  # 음수 또는 지나치게 작은 추정분산 방지
  sigma_sq_hat <- pmax(
    pred_var_full,
    0.01
  )
  
  
  # ----------------------------------------------------------
  # Step 4-5. OPT 포함확률 계산
  # ----------------------------------------------------------
  
  # omega_i
  #
  #   omega_i
  #     = sigma_hat^2(x_i)
  #       * x_i^T A_hat^{-1} x_i
  
  omega_i <- (
    sigma_sq_hat *
      leverage_part
  )
  
  
  # 이분산 조건에서의 최적 포함확률
  #
  #   pi_i^OPT ∝ sqrt(omega_i)
  #
  # 파일럿 표본도 최종 표본에 다시 포함될 수 있으므로
  # 모집단 전체 N명을 대상으로 포함확률을 계산한다.
  #
  # sum(pi_opt) = n_total이므로
  # 최종 OPT 표본의 기대 크기는 200명이다.
  
  pi_opt <- (
    n_total * sqrt(omega_i) /
      sum(sqrt(omega_i))
  )
  
  
  # ----------------------------------------------------------
  # Step 6-1. OPT final sampling
  # ----------------------------------------------------------
  
  # 모집단 전체에서 최종 표본을 추출하므로
  # 파일럿 표본이 다시 선택될 수도 있다.
  
  delta_opt <- rbinom(
    n = N,
    size = 1,
    prob = pi_opt
  )
  
  
  d_opt <- full_data_hetero[
    delta_opt == 1,
  ]
  
  
  pi_opt_selected <- pi_opt[
    delta_opt == 1
  ]
  
  # ----------------------------------------------------------
  # Step 6-2. SRS 비교군
  # ----------------------------------------------------------
  
  # 모집단 전체 N명에게 동일한 포함확률을 부여한다.
  #
  # 파일럿 표본도 SRS 최종 표본에 다시 포함될 수 있다.
  #
  # sum(pi_srs) = n_total이므로
  # SRS의 기대 표본크기도 200명이다.
  
  pi_srs <- rep(
    n_total / N,
    N
  )
  
  
  delta_srs <- rbinom(
    n = N,
    size = 1,
    prob = pi_srs
  )
  
  
  d_srs <- full_data_hetero[
    delta_srs == 1,
  ]
  
  # ----------------------------------------------------------
  # Step 7. 최종 회귀계수 추정
  # ----------------------------------------------------------
  
  # OPT: IPW 가중회귀
  fit_opt <- lm(
    y ~ . - 1,
    data = d_opt,
    weights = 1 / pi_opt_selected
  )
  
  
  beta_opt <- coef(fit_opt)
  
  
  # SRS는 포함확률이 모두 같으므로
  # 일반 최소제곱법과 IPW 결과가 동일하다.
  
  fit_srs <- lm(
    y ~ . - 1,
    data = d_srs
  )
  
  
  beta_srs <- coef(fit_srs)
  
  
  
  # ----------------------------------------------------------
  # Parameter squared error
  # ----------------------------------------------------------
  #
  # 한 반복에서 7개 회귀계수의 제곱오차를 모두 더한다.
  #
  # 이후 반복 결과를 평균내어 최종 Parameter MSE를 구한다.
  
  results_opt$param_mse[r] <- sum(
    (beta_opt - beta_true)^2
  )
  
  
  results_srs$param_mse[r] <- sum(
    (beta_srs - beta_true)^2
  )
  
  
  # ----------------------------------------------------------
  # Prediction MSE
  # ----------------------------------------------------------
  
  y_pred_opt <- drop(
    X_test_mat %*% beta_opt
  )
  
  
  y_pred_srs <- drop(
    X_test_mat %*% beta_srs
  )
  
  
  results_opt$pred_mse[r] <- mean(
    (y_test_vec - y_pred_opt)^2
  )
  
  
  results_srs$pred_mse[r] <- mean(
    (y_test_vec - y_pred_srs)^2
  )
  
  
  # 진행 상황 출력
  if (r %% 100L == 0L) {
    cat(
      r,
      "번째 반복 완료...\n"
    )
  }
}


# ------------------------------------------------------------
# 7. OPT 및 SRS 결과 요약
# ------------------------------------------------------------

final_table <- data.frame(
  
  Method = c(
    "OPT",
    "SRS"
  ),
  
  Param_MSE = c(
    mean(
      results_opt$param_mse,
      na.rm = TRUE
    ),
    mean(
      results_srs$param_mse,
      na.rm = TRUE
    )
  ),
  
  Pred_MSE = c(
    mean(
      results_opt$pred_mse,
      na.rm = TRUE
    ),
    mean(
      results_srs$pred_mse,
      na.rm = TRUE
    )
  )
)


# ------------------------------------------------------------
# 8. FULL model benchmark
# ------------------------------------------------------------

fit_full <- lm(
  y ~ . - 1,
  data = full_data_hetero
)


beta_full <- coef(fit_full)


# FULL parameter squared error
param_mse_full <- sum(
  (beta_full - beta_true)^2
)


# FULL 테스트 예측값
pred_full <- drop(
  X_test_mat %*% beta_full
)


# FULL prediction MSE
pred_mse_full <- mean(
  (y_test_vec - pred_full)^2
)


# FULL 결과 추가
final_table <- rbind(
  final_table,
  data.frame(
    Method = "FULL",
    Param_MSE = param_mse_full,
    Pred_MSE = pred_mse_full
  )
)


# ------------------------------------------------------------
# 9. 최종 결과 출력
# ------------------------------------------------------------

print(final_table)

