# ============================================================
# California Housing Real Data Application
#
# 1. Homoscedastic design
#    - Final expected sample size: n = 200
#    - Leverage score만 이용
#
# 2. Heteroscedastic design
#    - Final expected sample size: n = 200
#    - Pilot sample size: n0 = 100
#    - Random Forest로 조건부 분산 추정

library(randomForest)
library(ggplot2)

set.seed(456)


# ------------------------------------------------------------
# 0. California Housing 데이터 불러오기
# ------------------------------------------------------------

url <- paste0(
  "https://raw.githubusercontent.com/ageron/",
  "handson-ml/master/datasets/housing/housing.csv"
)

housing <- read.csv(url)


# 결측값 제거
housing <- na.omit(housing)


# 분석에 사용할 변수만 선택
housing <- housing[, c(
  "median_house_value",
  "longitude",
  "latitude",
  "housing_median_age",
  "total_rooms",
  "total_bedrooms",
  "population",
  "households",
  "median_income"
)]


# 반응변수의 단위를 줄임
housing$y <- housing$median_house_value / 100000

housing$median_house_value <- NULL



# ------------------------------------------------------------
# 1. Train / Test split
# ------------------------------------------------------------

N_all <- nrow(housing)


# 훈련자료 10,000개 추출
idx_train <- sample(
  x = seq_len(N_all),
  size = 10000,
  replace = FALSE
)


train_data_raw <- housing[
  idx_train,
]

test_data_raw <- housing[
  -idx_train,
]


# 테스트 자료가 10,000개보다 많으면 10,000개만 사용
if (nrow(test_data_raw) > 10000) {
  
  idx_test <- sample(
    x = seq_len(nrow(test_data_raw)),
    size = 10000,
    replace = FALSE
  )
  
  test_data_raw <- test_data_raw[
    idx_test,
  ]
}



# ------------------------------------------------------------
# 2. 설명변수 표준화
# ------------------------------------------------------------

x_cols <- setdiff(
  colnames(train_data_raw),
  "y"
)


# 훈련자료의 평균과 표준편차를 사용
x_mean <- sapply(
  train_data_raw[, x_cols],
  mean
)

x_sd <- sapply(
  train_data_raw[, x_cols],
  sd
)


train_X_scaled <- scale(
  train_data_raw[, x_cols],
  center = x_mean,
  scale = x_sd
)

test_X_scaled <- scale(
  test_data_raw[, x_cols],
  center = x_mean,
  scale = x_sd
)


# 절편을 포함한 설계행렬
X_pop <- cbind(
  intercept = 1,
  train_X_scaled
)

X_test <- cbind(
  intercept = 1,
  test_X_scaled
)


p <- ncol(X_pop) - 1
N <- nrow(X_pop)
N_test <- nrow(X_test)


# 열 이름 통일
colnames(X_pop) <- c(
  "intercept",
  paste0("x", seq_len(p))
)

colnames(X_test) <- c(
  "intercept",
  paste0("x", seq_len(p))
)


x_names <- colnames(X_pop)


# 분석용 데이터프레임
full_data <- data.frame(
  y = train_data_raw$y,
  X_pop
)

test_data <- data.frame(
  y = test_data_raw$y,
  X_test
)


X_test_mat <- as.matrix(
  test_data[, -1]
)

y_test_vec <- test_data$y



# ------------------------------------------------------------
# 3. 공통 설정 및 레버리지 계산
# ------------------------------------------------------------

n_reps <- 100


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


# 수치오차로 인해 매우 작은 음수가 발생하는 경우 방지
leverage_part <- pmax(
  leverage_part,
  0
)



# ============================================================
# 4. Homoscedastic Design
# ============================================================
#
# 조건부 분산을 따로 추정하지 않고
#
#   pi_i^OPT ∝ sqrt(x_i^T A_hat^{-1} x_i)
#
# 를 사용한다.
#
# 최종 OPT와 SRS의 기대 표본크기는 모두 200이다.
# ============================================================

n_homo <- 200


# ------------------------------------------------------------
# 4-1. OPT 포함확률
# ------------------------------------------------------------

pi_opt_homo <- (
  n_homo * sqrt(leverage_part) /
    sum(sqrt(leverage_part))
)


# ------------------------------------------------------------
# 4-2. SRS 포함확률
# ------------------------------------------------------------

pi_srs_homo <- rep(
  n_homo / N,
  N
)


# ------------------------------------------------------------
# 4-3. 결과 저장 벡터
# ------------------------------------------------------------

pred_mse_homo_opt <- rep(
  NA_real_,
  n_reps
)

pred_mse_homo_srs <- rep(
  NA_real_,
  n_reps
)


sample_size_homo_opt <- rep(
  NA_integer_,
  n_reps
)

sample_size_homo_srs <- rep(
  NA_integer_,
  n_reps
)



# ------------------------------------------------------------
# 4-4. 등분산 설계 반복
# ------------------------------------------------------------

for (r in seq_len(n_reps)) {
  
  
  # ----------------------------------------------------------
  # OPT 표본추출
  # ----------------------------------------------------------
  
  delta_opt <- rbinom(
    n = N,
    size = 1,
    prob = pi_opt_homo
  )
  
  
  selected_opt <- which(
    delta_opt == 1
  )
  
  
  sample_size_homo_opt[r] <- length(
    selected_opt
  )
  
  
  # ----------------------------------------------------------
  # SRS 표본추출
  # ----------------------------------------------------------
  
  delta_srs <- rbinom(
    n = N,
    size = 1,
    prob = pi_srs_homo
  )
  
  
  selected_srs <- which(
    delta_srs == 1
  )
  
  
  sample_size_homo_srs[r] <- length(
    selected_srs
  )
  
  # 추출자료
  d_opt <- full_data[
    selected_opt,
  ]
  
  d_srs <- full_data[
    selected_srs,
  ]
  
  
  # ----------------------------------------------------------
  # 최종 회귀모형 적합
  # ----------------------------------------------------------
  
  # OPT: IPW 가중회귀
  fit_opt <- lm(
    y ~ . - 1,
    data = d_opt,
    weights = 1 / pi_opt_homo[selected_opt]
  )
  
  
  # SRS: 모든 포함확률이 같으므로 일반 OLS 사용
  fit_srs <- lm(
    y ~ . - 1,
    data = d_srs
  )
  
  
  beta_opt <- coef(fit_opt)[x_names]
  beta_srs <- coef(fit_srs)[x_names]
  
  
  # ----------------------------------------------------------
  # Test Prediction MSE
  # ----------------------------------------------------------
  
  y_pred_opt <- drop(
    X_test_mat %*% beta_opt
  )
  
  y_pred_srs <- drop(
    X_test_mat %*% beta_srs
  )
  
  
  pred_mse_homo_opt[r] <- mean(
    (y_test_vec - y_pred_opt)^2
  )
  
  pred_mse_homo_srs[r] <- mean(
    (y_test_vec - y_pred_srs)^2
  )
  
  
  if (r %% 20L == 0L) {
    
    cat(
      "[Homoscedastic]",
      r,
      "번째 반복 완료\n"
    )
  }
}



# ============================================================
# 5. Heteroscedastic Design
# ============================================================
#
# Pilot sample n0 = 100을 이용하여
# Random Forest로 sigma^2(x)를 추정한다.
#
# 이후
#
#   pi_i^OPT
#     ∝ sqrt{
#         sigma_hat^2(x_i)
#         x_i^T A_hat^{-1} x_i
#       }
#
# 를 사용한다.
#
# 최종 OPT와 SRS의 기대 표본크기는 모두 200이다.
#
# 파일럿 표본은 최종 표본에 다시 선택될 수 있다.
# ============================================================

n_hetero <- 200
n_pilot <- 100


# SRS 포함확률
pi_srs_hetero <- rep(
  n_hetero / N,
  N
)


# 결과 저장 벡터
pred_mse_hetero_opt <- rep(
  NA_real_,
  n_reps
)

pred_mse_hetero_srs <- rep(
  NA_real_,
  n_reps
)


sample_size_hetero_opt <- rep(
  NA_integer_,
  n_reps
)

sample_size_hetero_srs <- rep(
  NA_integer_,
  n_reps
)



# ------------------------------------------------------------
# 5-1. 이분산 설계 반복
# ------------------------------------------------------------

for (r in seq_len(n_reps)) {
  
  
  # ----------------------------------------------------------
  # Step 1. Pilot sampling
  # ----------------------------------------------------------
  
  idx_pilot <- sample(
    x = seq_len(N),
    size = n_pilot,
    replace = FALSE
  )
  
  
  d_pilot <- full_data[
    idx_pilot,
  ]
  
  
  # 파일럿 선형회귀모형
  fit_pilot <- lm(
    y ~ . - 1,
    data = d_pilot
  )
  
  
  # ----------------------------------------------------------
  # Step 2. Random Forest로 조건부 분산 추정
  # ----------------------------------------------------------
  
  # 파일럿 표본에서 잔차제곱 계산
  pilot_res_sq <- residuals(
    fit_pilot
  )^2
  
  
  # RF 학습자료
  #
  # y와 intercept는 제외하고
  # 실제 설명변수 x1~xp만 사용한다.
  
  d_pilot_var <- data.frame(
    res_sq = pilot_res_sq,
    d_pilot[, paste0("x", seq_len(p))]
  )
  
  
  # Random Forest로 조건부 분산함수 추정
  var_rf_model <- randomForest(
    res_sq ~ .,
    data = d_pilot_var,
    ntree = 100
  )
  
  
  # 모집단 전체에 대해 조건부 분산 예측
  pred_var_full <- predict(
    var_rf_model,
    newdata = full_data[
      ,
      paste0("x", seq_len(p))
    ]
  )
  
  
  # 음수 또는 지나치게 작은 분산 예측값 방지
  sigma_sq_hat <- pmax(
    pred_var_full,
    0.01
  )
  
  
  # ----------------------------------------------------------
  # Step 3. OPT 포함확률 계산
  # ----------------------------------------------------------
  
  omega_i <- (
    sigma_sq_hat *
      leverage_part
  )
  
  
  pi_opt_hetero <- (
    n_hetero * sqrt(omega_i) /
      sum(sqrt(omega_i))
  )
  
  # ----------------------------------------------------------
  # Step 4. 최종 OPT 표본추출
  # ----------------------------------------------------------
  
  delta_opt <- rbinom(
    n = N,
    size = 1,
    prob = pi_opt_hetero
  )
  
  
  selected_opt <- which(
    delta_opt == 1
  )
  
  
  sample_size_hetero_opt[r] <- length(
    selected_opt
  )
  
  
  # ----------------------------------------------------------
  # Step 5. 최종 SRS 표본추출
  # ----------------------------------------------------------
  
  delta_srs <- rbinom(
    n = N,
    size = 1,
    prob = pi_srs_hetero
  )
  
  
  selected_srs <- which(
    delta_srs == 1
  )
  
  
  sample_size_hetero_srs[r] <- length(
    selected_srs
  )
  
  d_opt <- full_data[
    selected_opt,
  ]
  
  d_srs <- full_data[
    selected_srs,
  ]
  
  
  # ----------------------------------------------------------
  # Step 6. 최종 회귀모형 적합
  # ----------------------------------------------------------
  
  # OPT: IPW 가중회귀
  fit_opt <- lm(
    y ~ . - 1,
    data = d_opt,
    weights = 1 / pi_opt_hetero[selected_opt]
  )
  
  
  # SRS: 일반 OLS
  fit_srs <- lm(
    y ~ . - 1,
    data = d_srs
  )
  
  
  beta_opt <- coef(fit_opt)[x_names]
  beta_srs <- coef(fit_srs)[x_names]
  
  
  # ----------------------------------------------------------
  # Step 7. Test Prediction MSE
  # ----------------------------------------------------------
  
  y_pred_opt <- drop(
    X_test_mat %*% beta_opt
  )
  
  y_pred_srs <- drop(
    X_test_mat %*% beta_srs
  )
  
  
  pred_mse_hetero_opt[r] <- mean(
    (y_test_vec - y_pred_opt)^2
  )
  
  pred_mse_hetero_srs[r] <- mean(
    (y_test_vec - y_pred_srs)^2
  )
  
  
  if (r %% 20L == 0L) {
    
    cat(
      "[Heteroscedastic]",
      r,
      "번째 반복 완료\n"
    )
  }
}



# ============================================================
# 6. 평균 Prediction MSE 및 평균 표본크기
# ============================================================

real_results <- data.frame(
  
  Scenario = c(
    "Homoscedastic",
    "Homoscedastic",
    "Heteroscedastic",
    "Heteroscedastic"
  ),
  
  Method = c(
    "OPT",
    "SRS",
    "OPT",
    "SRS"
  ),
  
  Pred_MSE = c(
    mean(pred_mse_homo_opt),
    mean(pred_mse_homo_srs),
    mean(pred_mse_hetero_opt),
    mean(pred_mse_hetero_srs)
  ),
  
  Mean_Sample_Size = c(
    mean(sample_size_homo_opt),
    mean(sample_size_homo_srs),
    mean(sample_size_hetero_opt),
    mean(sample_size_hetero_srs)
  )
)


print(real_results)



# ============================================================
# 7. Improvement 및 Win rate
# ============================================================

# ------------------------------------------------------------
# 7-1. Homoscedastic
# ------------------------------------------------------------

improvement_homo <- (
  100 *
    (
      pred_mse_homo_srs -
        pred_mse_homo_opt
    ) /
    pred_mse_homo_srs
)


wins_homo <- sum(
  pred_mse_homo_opt <
    pred_mse_homo_srs
)


win_rate_homo <- (
  100 *
    wins_homo /
    n_reps
)


cat("\n[Homoscedastic design]\n")


cat(
  "Mean improvement (%):",
  mean(improvement_homo),
  "\n"
)

cat(
  "Median improvement (%):",
  median(improvement_homo),
  "\n"
)

cat(
  "Wins:",
  wins_homo,
  "\n"
)

cat(
  "Win rate (%):",
  win_rate_homo,
  "\n"
)



# ------------------------------------------------------------
# 7-2. Heteroscedastic
# ------------------------------------------------------------

improvement_hetero <- (
  100 *
    (
      pred_mse_hetero_srs -
        pred_mse_hetero_opt
    ) /
    pred_mse_hetero_srs
)


wins_hetero <- sum(
  pred_mse_hetero_opt <
    pred_mse_hetero_srs
)


win_rate_hetero <- (
  100 *
    wins_hetero /
    n_reps
)


cat("\n[Heteroscedastic design]\n")

cat(
  "Mean improvement (%):",
  mean(improvement_hetero),
  "\n"
)

cat(
  "Median improvement (%):",
  median(improvement_hetero),
  "\n"
)

cat(
  "Wins:",
  wins_hetero,
  "\n"
)

cat(
  "Win rate (%):",
  win_rate_hetero,
  "\n"
)

