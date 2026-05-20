# --- 0. 환경 설정 및 시드 고정 ---
set.seed(123) 
N <- 10000     # 모집단 크기
p <- 6         # 변수 차원

# --- 1. 훈련용 모집단(Population) 생성 ---
# 독립변수 X 생성 (교수님 논문 3.1절 분포 응용)
x1 <- rnorm(N, 0, 1)                  # 연속형 (정규분포)
x2 <- rnorm(N, 0, 1)                  # 연속형 (정규분포)
x3 <- sample(1:5, N, replace = TRUE)  # 이산형 (1~5 균등)
x4 <- sample(1:5, N, replace = TRUE)  # 이산형 (1~5 균등)
x5 <- rbinom(N, 1, 0.5)               # 범주형 (흔함 p=0.5)
x6 <- rbinom(N, 1, 0.1)               # 범주형 (희귀 p=0.1, low prevalence)

# 설계행렬 X 구성 (절편항 포함 총 7개 컬럼)
X_pop <- cbind(1, x1, x2, x3, x4, x5, x6)

# 참값(Beta_0) 설정: 영향력이 작은 0.1 포함
beta_true <- c(1.0, 1.5, 1.0, 0.5, 0.1, -1.0, 2.0) 

# 등분산 오차항 생성 및 Y 합성
error_pop <- rnorm(N, 0, 1) 
y_pop <- X_pop %*% beta_true + error_pop

# 훈련용 데이터프레임
full_data <- data.frame(y = y_pop, X_pop)
colnames(full_data) <- c("y", "intercept", paste0("x", 1:p))




# --- 2. 테스트 데이터(Test Set) 생성 ---
# 훈련 데이터와 동일한 메커니즘으로 생성하되, 별도의 샘플임
N_test <- 10000
x1_t <- rnorm(N_test, 0, 1); x2_t <- rnorm(N_test, 0, 1)
x3_t <- sample(1:5, N_test, replace = TRUE); x4_t <- sample(1:5, N_test, replace = TRUE)
x5_t <- rbinom(N_test, 1, 0.5); x6_t <- rbinom(N_test, 1, 0.1)

X_test <- cbind(1, x1_t, x2_t, x3_t, x4_t, x5_t, x6_t)
error_test <- rnorm(N_test, 0, 1)
y_test <- X_test %*% beta_true + error_test

# 테스트용 데이터프레임
test_data <- data.frame(y = y_test, X_test)
colnames(test_data) <- c("y", "intercept", paste0("x", 1:p))

# 생성 확인
cat("훈련용 모집단 생성 완료:", nrow(full_data), "행\n")
cat("테스트 데이터 생성 완료:", nrow(test_data), "행\n")


# --- 3. 샘플링 확률 계산 (n = 200 기준으로 예시) ---
n <- 200  

# (1) SRS & MLE 확률: 모든 데이터에 동일한 확률 부여
pi_srs <- rep(n / N, N)

# (2) OPT(최적 샘플링) 확률: 레버리지 스코어 기반
# 레버리지 스코어 h_i = x_i^T * (X^T * X)^-1 * x_i 계산
# 여기서 X_pop는 아까 만든 10000x7 행렬
X_mat <- as.matrix(X_pop)
N <- nrow(X_mat)
A_hat <- (t(X_mat) %*% X_mat) / N 
A_inv <- solve(A_hat)

# 각 관측치별 레버리지 스코어 계산
leverage_scores <- rowSums((X_mat %*% A_inv) * X_mat)

# 논문의 Algorithm 1: pi_i는 레버리지 스코어에 비례함
# pi_i = n * (h_i / sum(h_i))
pi_opt <- n * (sqrt(leverage_scores) / sum(sqrt(leverage_scores)))

# --- 확률값 보정 ---
# 확률은 1을 넘을 수 없으므로, 1보다 큰 값은 1로 고정(Capping)해줘야 에러가 안 남
pi_opt <- pmin(pi_opt, 1)
summary(pi_opt)

# --- 결과 저장용 행렬 생성 ---
n_reps <- 1000
results <- data.frame(
  Method = rep(c("OPT", "SRS", "NAIVE_OPT", "FULL"), each = n_reps),
  Param_MSE = 0,
  Pred_MSE = 0
)

# --- 반복 시뮬레이션 시작 ---
for (i in 1:n_reps) {
  
  # (A) 샘플링 지시변수(delta) 생성: 베르누이 시행
  delta_opt <- rbinom(N, 1, pi_opt)
  delta_srs <- rbinom(N, 1, pi_srs)
  
  # (B) 데이터 추출
  # OPT와 SRS는 추출 확률(pi)이 다르므로 각각의 인덱스로 추출
  d_opt <- full_data[delta_opt == 1, ]
  d_srs <- full_data[delta_srs == 1, ]
  
  # (C) 파라미터 추정 (Beta_hat)
  
  # 1. OPT: IPW 가중치 적용 (1/pi_opt)
  # 가중치 계산 시 추출된 데이터의 확률만 가져옴
  w_opt <- 1 / pi_opt[delta_opt == 1]
  fit_opt <- lm(y ~ . - 1 - intercept + intercept, data = d_opt, weights = w_opt)
  beta_opt <- coef(fit_opt)
  
  # 2. SRS: IPW 가중치 적용 (사실상 동일 가중치라 평범한 lm과 같음)
  w_srs <- 1 / pi_srs[delta_srs == 1]
  fit_srs <- lm(y ~ . - 1 - intercept + intercept, data = d_srs, weights = w_srs)
  beta_srs <- coef(fit_srs)
  
  # 3. naive: 최적 샘플링으로 뽑았으나 가중치를 무시한 경우 (Naive)
  fit_naive <- lm(y ~ . - 1 - intercept + intercept, data = d_opt)
  beta_naive <- coef(fit_naive)
  
  # 4. FULL: 모집단 전체를 사용 (Benchmark)
  if(i == 1) { # FULL은 한 번만 계산하면 됨
    fit_full <- lm(y ~ . - 1 - intercept + intercept, data = full_data)
    beta_full <- coef(fit_full)
  }
  
  # (D) MSE 계산 및 저장
  
  # 1) Parameter MSE: ||beta_hat - beta_true||^2
  results$Param_MSE[i] <- sum((beta_opt - beta_true)^2)
  results$Param_MSE[i + n_reps] <- sum((beta_srs - beta_true)^2)
  results$Param_MSE[i + 2*n_reps] <- sum((beta_naive - beta_true)^2)
  results$Param_MSE[i + 3*n_reps] <- sum((beta_full - beta_true)^2)
  
  # 2) Prediction MSE: Test 데이터에 적용
  # y_hat = X_test %*% beta_hat
  pred_opt <- as.matrix(test_data[,-1]) %*% beta_opt
  pred_srs <- as.matrix(test_data[,-1]) %*% beta_srs
  pred_naive <- as.matrix(test_data[,-1]) %*% beta_naive
  pred_full <- as.matrix(test_data[,-1]) %*% beta_full
  
  results$Pred_MSE[i] <- mean((test_data$y - pred_opt)^2)
  results$Pred_MSE[i + n_reps] <- mean((test_data$y - pred_srs)^2)
  results$Pred_MSE[i + 2*n_reps] <- mean((test_data$y - pred_naive)^2)
  results$Pred_MSE[i + 3*n_reps] <- mean((test_data$y - pred_full)^2)
  
  # 진행 상황 출력 (100번마다)
  if(i %% 100 == 0) cat(i, "번째 반복 완료...\n")
}

# 방법론별 평균 MSE 확인
summary_table <- aggregate(cbind(Param_MSE, Pred_MSE) ~ Method, data = results, FUN = mean)
print(summary_table)

# --- 4. 샘플링 확률 계산 (n = 400) ---
n <- 400  

# (1) SRS  확률: 모든 데이터에 동일한 확률 부여
pi_srs <- rep(n / N, N)

# (2) OPT(최적 샘플링) 확률: 레버리지 스코어 기반
# 레버리지 스코어 h_i = x_i^T * (X^T * X)^-1 * x_i 계산
# 여기서 X_pop는 아까 만든 10000x7 행렬
X_mat <- as.matrix(X_pop)
N <- nrow(X_mat)
A_hat <- (t(X_mat) %*% X_mat) / N 
A_inv <- solve(A_hat)

# 각 관측치별 레버리지 스코어 계산
leverage_scores <- rowSums((X_mat %*% A_inv) * X_mat)

# 네 논문의 Algorithm 1: pi_i는 레버리지 스코어에 비례함
# pi_i = n * (h_i / sum(h_i))
pi_opt <- n * (sqrt(leverage_scores) / sum(sqrt(leverage_scores)))

# --- 확률값 보정 ---
# 확률은 1을 넘을 수 없으므로, 1보다 큰 값은 1로 고정(Capping)해줘야 에러가 안 남
pi_opt <- pmin(pi_opt, 1)
summary(pi_opt)

# --- 결과 저장용 행렬 생성 ---
n_reps <- 1000
results <- data.frame(
  Method = rep(c("OPT", "SRS", "NAIVE_OPT", "FULL"), each = n_reps),
  Param_MSE = 0,
  Pred_MSE = 0
)

# --- 반복 시뮬레이션 시작 ---
for (i in 1:n_reps) {
  
  # (A) 샘플링 지시변수(delta) 생성: 베르누이 시행
  delta_opt <- rbinom(N, 1, pi_opt)
  delta_srs <- rbinom(N, 1, pi_srs)
  
  # (B) 데이터 추출
  # OPT와 SRS는 추출 확률(pi)이 다르므로 각각의 인덱스로 추출
  d_opt <- full_data[delta_opt == 1, ]
  d_srs <- full_data[delta_srs == 1, ]
  
  # (C) 파라미터 추정 (Beta_hat)
  
  # 1. OPT: IPW 가중치 적용 (1/pi_opt)
  # 가중치 계산 시 추출된 데이터의 확률만 가져옴
  w_opt <- 1 / pi_opt[delta_opt == 1]
  fit_opt <- lm(y ~ . - 1 - intercept + intercept, data = d_opt, weights = w_opt)
  beta_opt <- coef(fit_opt)
  
  # 2. SRS: IPW 가중치 적용 (사실상 동일 가중치라 평범한 lm과 같음)
  w_srs <- 1 / pi_srs[delta_srs == 1]
  fit_srs <- lm(y ~ . - 1 - intercept + intercept, data = d_srs, weights = w_srs)
  beta_srs <- coef(fit_srs)
  
  # 3. NAIVE: 최적 샘플링으로 뽑았으나 가중치를 무시한 경우 (Naive)
  fit_naive <- lm(y ~ . - 1 - intercept + intercept, data = d_opt)
  beta_naive <- coef(fit_naive)
  
  # 4. FULL: 모집단 전체를 사용 (Benchmark)
  if(i == 1) { # FULL은 한 번만 계산하면 됨
    fit_full <- lm(y ~ . - 1 - intercept + intercept, data = full_data)
    beta_full <- coef(fit_full)
  }
  
  # (D) MSE 계산 및 저장
  
  # 1) Parameter MSE: ||beta_hat - beta_true||^2
  results$Param_MSE[i] <- sum((beta_opt - beta_true)^2)
  results$Param_MSE[i + n_reps] <- sum((beta_srs - beta_true)^2)
  results$Param_MSE[i + 2*n_reps] <- sum((beta_naive - beta_true)^2)
  results$Param_MSE[i + 3*n_reps] <- sum((beta_full - beta_true)^2)
  
  # 2) Prediction MSE: Test 데이터에 적용
  # y_hat = X_test %*% beta_hat
  pred_opt <- as.matrix(test_data[,-1]) %*% beta_opt
  pred_srs <- as.matrix(test_data[,-1]) %*% beta_srs
  pred_naive <- as.matrix(test_data[,-1]) %*% beta_naive
  pred_full <- as.matrix(test_data[,-1]) %*% beta_full
  
  results$Pred_MSE[i] <- mean((test_data$y - pred_opt)^2)
  results$Pred_MSE[i + n_reps] <- mean((test_data$y - pred_srs)^2)
  results$Pred_MSE[i + 2*n_reps] <- mean((test_data$y - pred_naive)^2)
  results$Pred_MSE[i + 3*n_reps] <- mean((test_data$y - pred_full)^2)
  
  # 진행 상황 출력 (100번마다)
  if(i %% 100 == 0) cat(i, "번째 반복 완료...\n")
}

# 방법론별 평균 MSE 확인
summary_table <- aggregate(cbind(Param_MSE, Pred_MSE) ~ Method, data = results, FUN = mean)
print(summary_table)


# --- 5. 샘플링 확률 계산 (n = 800) ---
n <- 800

# (1) SRS 확률: 모든 데이터에 동일한 확률 부여
pi_srs <- rep(n / N, N)

# (2) OPT(최적 샘플링) 확률: 레버리지 스코어 기반
# 레버리지 스코어 h_i = x_i^T * (X^T * X)^-1 * x_i 계산
# 여기서 X_pop는 아까 만든 10000x7 행렬
X_mat <- as.matrix(X_pop)
N <- nrow(X_mat)
A_hat <- (t(X_mat) %*% X_mat) / N 
A_inv <- solve(A_hat)

# 각 관측치별 레버리지 스코어 계산
leverage_scores <- rowSums((X_mat %*% A_inv) * X_mat)

# 네 논문의 Algorithm 1: pi_i는 레버리지 스코어에 비례함
# pi_i = n * (h_i / sum(h_i))
pi_opt <- n * (sqrt(leverage_scores) / sum(sqrt(leverage_scores)))

# --- 확률값 보정 (중요!) ---
# 확률은 1을 넘을 수 없으므로, 1보다 큰 값은 1로 고정(Capping)해줘야 에러가 안 남
pi_opt <- pmin(pi_opt, 1)
summary(pi_opt)

# --- 결과 저장용 행렬 생성 ---
n_reps <- 1000
results <- data.frame(
  Method = rep(c("OPT", "SRS", "NAIVE_OPT", "FULL"), each = n_reps),
  Param_MSE = 0,
  Pred_MSE = 0
)

# ---반복 시뮬레이션 시작 ---
for (i in 1:n_reps) {
  
  # (A) 샘플링 지시변수(delta) 생성: 베르누이 시행
  delta_opt <- rbinom(N, 1, pi_opt)
  delta_srs <- rbinom(N, 1, pi_srs)
  
  # (B) 데이터 추출
  # OPT와 SRS는 추출 확률(pi)이 다르므로 각각의 인덱스로 추출
  d_opt <- full_data[delta_opt == 1, ]
  d_srs <- full_data[delta_srs == 1, ]
  
  # (C) 파라미터 추정 (Beta_hat)
  
  # 1. OPT: IPW 가중치 적용 (1/pi_opt)
  # 가중치 계산 시 추출된 데이터의 확률만 가져옴
  w_opt <- 1 / pi_opt[delta_opt == 1]
  fit_opt <- lm(y ~ . - 1 - intercept + intercept, data = d_opt, weights = w_opt)
  beta_opt <- coef(fit_opt)
  
  # 2. SRS: IPW 가중치 적용 (사실상 동일 가중치라 평범한 lm과 같음)
  w_srs <- 1 / pi_srs[delta_srs == 1]
  fit_srs <- lm(y ~ . - 1 - intercept + intercept, data = d_srs, weights = w_srs)
  beta_srs <- coef(fit_srs)
  
  # 3. MLE: 최적 샘플링으로 뽑았으나 가중치를 무시한 경우 (Naive)
  fit_naive <- lm(y ~ . - 1 - intercept + intercept, data = d_opt)
  beta_naive <- coef(fit_naive)
  
  # 4. FULL: 모집단 전체를 사용 (Benchmark)
  if(i == 1) { # FULL은 한 번만 계산하면 됨
    fit_full <- lm(y ~ . - 1 - intercept + intercept, data = full_data)
    beta_full <- coef(fit_full)
  }
  
  # (D) MSE 계산 및 저장
  
  # 1) Parameter MSE: ||beta_hat - beta_true||^2
  results$Param_MSE[i] <- sum((beta_opt - beta_true)^2)
  results$Param_MSE[i + n_reps] <- sum((beta_srs - beta_true)^2)
  results$Param_MSE[i + 2*n_reps] <- sum((beta_naive - beta_true)^2)
  results$Param_MSE[i + 3*n_reps] <- sum((beta_full - beta_true)^2)
  
  # 2) Prediction MSE: Test 데이터에 적용
  # y_hat = X_test %*% beta_hat
  pred_opt <- as.matrix(test_data[,-1]) %*% beta_opt
  pred_srs <- as.matrix(test_data[,-1]) %*% beta_srs
  pred_naive <- as.matrix(test_data[,-1]) %*% beta_naive
  pred_full <- as.matrix(test_data[,-1]) %*% beta_full
  
  results$Pred_MSE[i] <- mean((test_data$y - pred_opt)^2)
  results$Pred_MSE[i + n_reps] <- mean((test_data$y - pred_srs)^2)
  results$Pred_MSE[i + 2*n_reps] <- mean((test_data$y - pred_naive)^2)
  results$Pred_MSE[i + 3*n_reps] <- mean((test_data$y - pred_full)^2)
  
  # 진행 상황 출력 (100번마다)
  if(i %% 100 == 0) cat(i, "번째 반복 완료...\n")
}

# 방법론별 평균 MSE 확인
summary_table <- aggregate(cbind(Param_MSE, Pred_MSE) ~ Method, data = results, FUN = mean)
print(summary_table)


#이분산 case
# --- 1. 환경 설정 및 시드 고정 ---
set.seed(456)
N <- 10000
p <- 6
n_total <- 200
n_pilot <- 100

# --- 2. 훈련용 모집단(Population) 생성 ---
x1 <- rnorm(N, 0, 1); x2 <- rnorm(N, 0, 1)
x3 <- sample(1:5, N, replace = TRUE); x4 <- sample(1:5, N, replace = TRUE)
x5 <- rbinom(N, 1, 0.5); x6 <- rbinom(N, 1, 0.1)
X_pop <- cbind(1, x1, x2, x3, x4, x5, x6)
beta_true <- c(1.0, 1.5, 1.0, 0.5, 0.1, -1.0, 2.0)

# [이분산 적용] x1, x2에 따라 오차 분산이 변함
sigma_vec <- abs(X_pop[,2] + X_pop[,3]) + 0.5 
error_hetero <- rnorm(N, mean = 0, sd = sigma_vec)
y_hetero <- X_pop %*% beta_true + error_hetero

full_data_hetero <- data.frame(y = y_hetero, X_pop)
colnames(full_data_hetero) <- c("y", "intercept", paste0("x", 1:p))


# --- 3. 테스트 데이터(Test Set) 생성 (이분산성 동일 적용) ---
N_test <- 10000
x1_t <- rnorm(N_test, 0, 1); x2_t <- rnorm(N_test, 0, 1)
x3_t <- sample(1:5, N_test, replace = TRUE); x4_t <- sample(1:5, N_test, replace = TRUE)
x5_t <- rbinom(N_test, 1, 0.5); x6_t <- rbinom(N_test, 1, 0.1)

X_test <- cbind(1, x1_t, x2_t, x3_t, x4_t, x5_t, x6_t)
# 테스트 셋도 동일한 이분산 구조를 가져야 함
sigma_vec_test <- abs(X_test[,2] + X_test[,3]) + 0.5
error_test <- rnorm(N_test, 0, sd = sigma_vec_test)
y_test <- X_test %*% beta_true + error_test

test_data_hetero <- data.frame(y = y_test, X_test)
colnames(test_data_hetero) <- c("y", "intercept", paste0("x", 1:p))


# --- 4. 레버리지 계산을 위한 A_hat (Step 1) ---
X_mat <- as.matrix(X_pop)
A_hat <- (t(X_mat) %*% X_mat) / N 
A_inv <- solve(A_hat)
leverage_part <- rowSums((X_mat %*% A_inv) * X_mat)

# --- 1. (n=200) ---
n_total <- 200   
n_pilot <- 100   
n_reps <- 1000   

# 결과 저장용 데이터프레임 
results_opt <- data.frame(param_mse = numeric(n_reps), pred_mse = numeric(n_reps))
results_srs <- data.frame(param_mse = numeric(n_reps), pred_mse = numeric(n_reps))
results_naive <- data.frame(param_mse = numeric(n_reps), pred_mse = numeric(n_reps))

X_test_mat <- as.matrix(test_data_hetero[, -1])
y_test_vec <- test_data_hetero$y

# --- 반복 시뮬레이션 루프 ---
for (r in 1:n_reps) {
  
  # [Step 2] Pilot Sampling
  idx_pilot <- sample(1:N, n_pilot)
  d_pilot <- full_data_hetero[idx_pilot, ]
  fit_pilot <- lm(y ~ . - 1, data = d_pilot)
  
  # [Step 3] 분산 추정
  y_hat_pilot <- predict(fit_pilot, full_data_hetero)
  sigma_sq_hat <- (full_data_hetero$y - y_hat_pilot)^2
  
  # [Step 4 & 5] 확률 계산
  omega_i <- sigma_sq_hat * leverage_part
  pi_opt <- n_total * (sqrt(omega_i) / sum(sqrt(omega_i)))
  pi_opt <- pmin(pi_opt, 1) 
  
  # [Step 6] Final Sampling
  delta_opt <- rbinom(N, 1, pi_opt)
  d_opt <- full_data_hetero[delta_opt == 1, ]
  pi_opt_selected <- pi_opt[delta_opt == 1]
  
  # SRS 비교군
  pi_srs_val <- n_total / N
  delta_srs <- rbinom(N, 1, rep(pi_srs_val, N))
  d_srs <- full_data_hetero[delta_srs == 1, ]
  
  # [Step 7] 최종 추정
  # 1. OPT (IPW 가중치 적용)
  fit_opt <- lm(y ~ . - 1, data = d_opt, weights = 1/pi_opt_selected)
  beta_opt <- coef(fit_opt)
  
  # 2. NAIVE (OPT 샘플을 사용하되 가중치 1/pi를 무시함)
  fit_naive <- lm(y ~ . - 1, data = d_opt)
  beta_naive <- coef(fit_naive)
  
  # 3. SRS
  fit_srs <- lm(y ~ . - 1, data = d_srs)
  beta_srs <- coef(fit_srs)
  
  # --- MSE 계산 및 저장 ---
  # Param_MSE
  results_opt$param_mse[r] <- sum((beta_opt - beta_true)^2)
  results_naive$param_mse[r] <- sum((beta_naive - beta_true)^2)
  results_srs$param_mse[r] <- sum((beta_srs - beta_true)^2)
  
  # Pred_MSE
  y_pred_opt <- X_test_mat %*% beta_opt
  y_pred_naive <- X_test_mat %*% beta_naive
  y_pred_srs <- X_test_mat %*% beta_srs
  
  results_opt$pred_mse[r] <- mean((y_test_vec - y_pred_opt)^2)
  results_naive$pred_mse[r] <- mean((y_test_vec - y_pred_naive)^2)
  results_srs$pred_mse[r] <- mean((y_test_vec - y_pred_srs)^2)
}

# --- 최종 결과 요약 ---
final_table <- data.frame(
  Method = c("OPT", "NAIVE_OPT", "SRS"),
  Param_MSE = c(mean(results_opt$param_mse), mean(results_naive$param_mse), mean(results_srs$param_mse)),
  Pred_MSE = c(mean(results_opt$pred_mse), mean(results_naive$pred_mse), mean(results_srs$pred_mse))
)

# FULL 모델 Benchmark
fit_full <- lm(y ~ . - 1, data = full_data_hetero)
beta_full <- coef(fit_full)
param_mse_full <- sum((beta_full - beta_true)^2)
pred_mse_full <- mean((y_test_vec - X_test_mat %*% beta_full)^2)

final_table <- rbind(final_table, data.frame(Method = "FULL", Param_MSE = param_mse_full, Pred_MSE = pred_mse_full))

print(final_table)


# ---2.(n=400) ---
n_total <- 400   
n_pilot <- 100   
n_reps <- 1000   

# 결과 저장용 데이터프레임 
results_opt <- data.frame(param_mse = numeric(n_reps), pred_mse = numeric(n_reps))
results_srs <- data.frame(param_mse = numeric(n_reps), pred_mse = numeric(n_reps))
results_naive <- data.frame(param_mse = numeric(n_reps), pred_mse = numeric(n_reps))

X_test_mat <- as.matrix(test_data_hetero[, -1])
y_test_vec <- test_data_hetero$y

# ---반복 시뮬레이션 루프 ---
for (r in 1:n_reps) {
  
  # [Step 2] Pilot Sampling
  idx_pilot <- sample(1:N, n_pilot)
  d_pilot <- full_data_hetero[idx_pilot, ]
  fit_pilot <- lm(y ~ . - 1, data = d_pilot)
  
  # [Step 3] 분산 추정
  y_hat_pilot <- predict(fit_pilot, full_data_hetero)
  sigma_sq_hat <- (full_data_hetero$y - y_hat_pilot)^2
  
  # [Step 4 & 5] 확률 계산
  omega_i <- sigma_sq_hat * leverage_part
  pi_opt <- n_total * (sqrt(omega_i) / sum(sqrt(omega_i)))
  pi_opt <- pmin(pi_opt, 1) 
  
  # [Step 6] Final Sampling
  delta_opt <- rbinom(N, 1, pi_opt)
  d_opt <- full_data_hetero[delta_opt == 1, ]
  pi_opt_selected <- pi_opt[delta_opt == 1]
  
  # SRS 비교군
  pi_srs_val <- n_total / N
  delta_srs <- rbinom(N, 1, rep(pi_srs_val, N))
  d_srs <- full_data_hetero[delta_srs == 1, ]
  
  # [Step 7] 최종 추정
  # 1. OPT (IPW 가중치 적용)
  fit_opt <- lm(y ~ . - 1, data = d_opt, weights = 1/pi_opt_selected)
  beta_opt <- coef(fit_opt)
  
  # 2. MLE (OPT 샘플을 사용하되 가중치 1/pi를 무시함)
  fit_naive <- lm(y ~ . - 1, data = d_opt)
  beta_naive <- coef(fit_naive)
  
  # 3. SRS
  fit_srs <- lm(y ~ . - 1, data = d_srs)
  beta_srs <- coef(fit_srs)
  
  # --- MSE 계산 및 저장 ---
  # Param_MSE
  results_opt$param_mse[r] <- sum((beta_opt - beta_true)^2)
  results_naive$param_mse[r] <- sum((beta_naive - beta_true)^2)
  results_srs$param_mse[r] <- sum((beta_srs - beta_true)^2)
  
  # Pred_MSE
  y_pred_opt <- X_test_mat %*% beta_opt
  y_pred_naive <- X_test_mat %*% beta_naive
  y_pred_srs <- X_test_mat %*% beta_srs
  
  results_opt$pred_mse[r] <- mean((y_test_vec - y_pred_opt)^2)
  results_naive$pred_mse[r] <- mean((y_test_vec - y_pred_naive)^2)
  results_srs$pred_mse[r] <- mean((y_test_vec - y_pred_srs)^2)
}

# ---최종 결과 요약 ---
final_table <- data.frame(
  Method = c("OPT", "NAIVE_OPT", "SRS"),
  Param_MSE = c(mean(results_opt$param_mse), mean(results_naive$param_mse), mean(results_srs$param_mse)),
  Pred_MSE = c(mean(results_opt$pred_mse), mean(results_naive$pred_mse), mean(results_srs$pred_mse))
)

# FULL 모델 Benchmark
fit_full <- lm(y ~ . - 1, data = full_data_hetero)
beta_full <- coef(fit_full)
param_mse_full <- sum((beta_full - beta_true)^2)
pred_mse_full <- mean((y_test_vec - X_test_mat %*% beta_full)^2)

final_table <- rbind(final_table, data.frame(Method = "FULL", Param_MSE = param_mse_full, Pred_MSE = pred_mse_full))

print(final_table)

# --- 3.(n=800) ---
n_total <- 800   
n_pilot <- 100   
n_reps <- 1000   

# 결과 저장용 데이터프레임
results_opt <- data.frame(param_mse = numeric(n_reps), pred_mse = numeric(n_reps))
results_srs <- data.frame(param_mse = numeric(n_reps), pred_mse = numeric(n_reps))
results_naive <- data.frame(param_mse = numeric(n_reps), pred_mse = numeric(n_reps))

X_test_mat <- as.matrix(test_data_hetero[, -1])
y_test_vec <- test_data_hetero$y

# --- 반복 시뮬레이션 루프 ---
for (r in 1:n_reps) {
  
  # [Step 2] Pilot Sampling
  idx_pilot <- sample(1:N, n_pilot)
  d_pilot <- full_data_hetero[idx_pilot, ]
  fit_pilot <- lm(y ~ . - 1, data = d_pilot)
  
  # [Step 3] 분산 추정
  y_hat_pilot <- predict(fit_pilot, full_data_hetero)
  sigma_sq_hat <- (full_data_hetero$y - y_hat_pilot)^2
  
  # [Step 4 & 5] 확률 계산
  omega_i <- sigma_sq_hat * leverage_part
  pi_opt <- n_total * (sqrt(omega_i) / sum(sqrt(omega_i)))
  pi_opt <- pmin(pi_opt, 1) 
  
  # [Step 6] Final Sampling
  delta_opt <- rbinom(N, 1, pi_opt)
  d_opt <- full_data_hetero[delta_opt == 1, ]
  pi_opt_selected <- pi_opt[delta_opt == 1]
  
  # SRS 비교군
  pi_srs_val <- n_total / N
  delta_srs <- rbinom(N, 1, rep(pi_srs_val, N))
  d_srs <- full_data_hetero[delta_srs == 1, ]
  
  # [Step 7] 최종 추정
  # 1. OPT (IPW 가중치 적용)
  fit_opt <- lm(y ~ . - 1, data = d_opt, weights = 1/pi_opt_selected)
  beta_opt <- coef(fit_opt)
  
  # 2. MLE (OPT 샘플을 사용하되 가중치 1/pi를 무시함)
  fit_naive <- lm(y ~ . - 1, data = d_opt)
  beta_naive <- coef(fit_naive)
  
  # 3. SRS
  fit_srs <- lm(y ~ . - 1, data = d_srs)
  beta_srs <- coef(fit_srs)
  
  # --- MSE 계산 및 저장 ---
  # Param_MSE
  results_opt$param_mse[r] <- sum((beta_opt - beta_true)^2)
  results_naive$param_mse[r] <- sum((beta_naive - beta_true)^2)
  results_srs$param_mse[r] <- sum((beta_srs - beta_true)^2)
  
  # Pred_MSE
  y_pred_opt <- X_test_mat %*% beta_opt
  y_pred_naive <- X_test_mat %*% beta_naive
  y_pred_srs <- X_test_mat %*% beta_srs
  
  results_opt$pred_mse[r] <- mean((y_test_vec - y_pred_opt)^2)
  results_naive$pred_mse[r] <- mean((y_test_vec - y_pred_naive)^2)
  results_srs$pred_mse[r] <- mean((y_test_vec - y_pred_srs)^2)
}

# --- 최종 결과 요약 ---
final_table <- data.frame(
  Method = c("OPT", "NAVIE_OPT", "SRS"),
  Param_MSE = c(mean(results_opt$param_mse), mean(results_naive$param_mse), mean(results_srs$param_mse)),
  Pred_MSE = c(mean(results_opt$pred_mse), mean(results_naive$pred_mse), mean(results_srs$pred_mse))
)

# FULL 모델 Benchmark
fit_full <- lm(y ~ . - 1, data = full_data_hetero)
beta_full <- coef(fit_full)
param_mse_full <- sum((beta_full - beta_true)^2)
pred_mse_full <- mean((y_test_vec - X_test_mat %*% beta_full)^2)

final_table <- rbind(final_table, data.frame(Method = "FULL", Param_MSE = param_mse_full, Pred_MSE = pred_mse_full))

print(final_table)
