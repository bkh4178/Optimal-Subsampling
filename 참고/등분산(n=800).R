# ============================================================
# 등분산 선형회귀 시뮬레이션
# 비교 방법: OPT vs SRS vs FULL
# ============================================================
#
# [평가기준]
#
# 1. Parameter MSE
#
#    Param_MSE
#      = mean((beta_hat - beta_true)^2)
#      = 절편을 포함한 전체 회귀계수들의 평균제곱오차
#
# 2. Prediction MSE
#
#    Pred_MSE
#      = mean((y_test - X_test %*% beta_hat)^2)
#
#    새로운 독립 테스트 자료에서 실제 반응변수 y_test와
#    예측값 사이의 평균제곱오차를 계산한다.
#
#
# [시뮬레이션 구조]
#
# 모집단과 테스트 자료는 처음에 한 번 생성하여 고정한다.
# 각 반복에서는 동일한 모집단으로부터 표본만 새로 추출한다.
#
# 따라서 이 시뮬레이션은 고정된 유한모집단에서
# 표본추출 설계로 인해 발생하는 변동을 비교한다.
# ============================================================


# ------------------------------------------------------------
# 0. 기본 설정
# ------------------------------------------------------------

set.seed(123)

N <- 10000          # 훈련용 유한모집단 크기
N_test <- 10000     # 독립 테스트 자료 크기
p <- 6              # 설명변수 개수, 절편 제외
n <- 800            # 목표 기대 표본크기
n_reps <- 1000      # 시뮬레이션 반복 횟수


# 설계행렬의 열 이름
x_names <- c(
  "intercept",
  paste0("x", seq_len(p))
)


# 참 회귀계수
#
# 계수에 이름을 붙여두면 추정계수의 순서가 변하더라도
# 잘못된 계수끼리 비교하는 문제를 방지할 수 있다.

beta_true <- c(
  intercept = 1.0,
  x1 = 1.5,
  x2 = 1.0,
  x3 = 0.5,
  x4 = 0.1,
  x5 = -1.0,
  x6 = 2.0
)

# ------------------------------------------------------------
# 1. 훈련용 유한모집단 생성
# ------------------------------------------------------------

# 연속형 설명변수
x1 <- rnorm(
  n = N,
  mean = 0,
  sd = 1
)

x2 <- rnorm(
  n = N,
  mean = 0,
  sd = 1
)


# 이산형 수치변수
# x3와 x4는 factor가 아니라 1, 2, 3, 4, 5의 값을 갖는
# 수치형 변수로 회귀모형에 들어간다.

x3 <- sample(
  x = 1:5,
  size = N,
  replace = TRUE
)

x4 <- sample(
  x = 1:5,
  size = N,
  replace = TRUE
)


# 이진 설명변수
x5 <- rbinom(
  n = N,
  size = 1,
  prob = 0.5
)

x6 <- rbinom(
  n = N,
  size = 1,
  prob = 0.1
)


# 절편을 포함한 N x 7 설계행렬
X_pop <- cbind(
  intercept = 1,
  x1 = x1,
  x2 = x2,
  x3 = x3,
  x4 = x4,
  x5 = x5,
  x6 = x6
)


# 등분산 오차항
#
# epsilon_i ~ N(0, 1)

error_pop <- rnorm(
  n = N,
  mean = 0,
  sd = 1
)


# 반응변수
#
# y = X beta + epsilon

y_pop <- drop(
  X_pop %*% beta_true + error_pop
)

# ------------------------------------------------------------
# 2. 독립 테스트 자료 생성
# ------------------------------------------------------------

# 훈련용 모집단과 동일한 자료생성 메커니즘을 사용하되
# 훈련용 모집단과 독립적으로 새로 생성한다.

x1_test <- rnorm(
  n = N_test,
  mean = 0,
  sd = 1
)

x2_test <- rnorm(
  n = N_test,
  mean = 0,
  sd = 1
)

x3_test <- sample(
  x = 1:5,
  size = N_test,
  replace = TRUE
)

x4_test <- sample(
  x = 1:5,
  size = N_test,
  replace = TRUE
)

x5_test <- rbinom(
  n = N_test,
  size = 1,
  prob = 0.5
)

x6_test <- rbinom(
  n = N_test,
  size = 1,
  prob = 0.1
)


# 절편을 포함한 테스트 설계행렬
X_test <- cbind(
  intercept = 1,
  x1 = x1_test,
  x2 = x2_test,
  x3 = x3_test,
  x4 = x4_test,
  x5 = x5_test,
  x6 = x6_test
)


# 테스트 자료의 등분산 오차항
error_test <- rnorm(
  n = N_test,
  mean = 0,
  sd = 1
)


# 테스트 자료의 반응변수
y_test <- drop(
  X_test %*% beta_true + error_test
)



# ------------------------------------------------------------
# 3. SRS 포함확률 계산
# ------------------------------------------------------------
#
# 모든 관측치에 동일한 포함확률을 부여한다.
#
# 각 관측치는 서로 독립적으로
#
#   delta_i ~ Bernoulli(n / N)
#
# 에 따라 추출된다.
#
# 따라서 실제 추출 표본크기는 반복마다 정확히 n이 아닐 수 있지만,
#
#   E[표본크기] = sum(pi_srs) = n
#
# 이므로 평균적으로 n명 정도가 추출된다.

pi_srs <- rep(
  x = n / N,
  times = N
)


# ------------------------------------------------------------
# 4. OPT 포함확률 계산
# ------------------------------------------------------------

# 모집단 전체의 설명변수를 이용하여
#
#   A_hat = (1/N) X^T X
#
# 를 계산한다.

A_hat <- crossprod(X_pop) / N


# A_hat의 역행렬 계산
A_inv <- solve(A_hat)


# 각 관측치의 leverage score 계산
#
#   omega_i = x_i^T A_hat^{-1} x_i
#
# 아래 식은 각 행별로 x_i^T A_hat^{-1} x_i를 계산한다.

leverage_scores <- rowSums(
  (X_pop %*% A_inv) * X_pop
)

# 등분산 조건에서의 최적 포함확률
#
#   pi_i^opt ∝ sqrt(omega_i)
#
# 비례상수를 조정하여
#
#   sum(pi_opt) = n
#
# 이 되도록 한다.
#
# 따라서 OPT 역시 실제 표본크기가 항상 n인 것은 아니지만,
# 기대 표본크기는 n이다.

pi_opt <- (
  n * sqrt(leverage_scores) /
    sum(sqrt(leverage_scores))
)

# ------------------------------------------------------------
# 5. IPW 회귀계수 추정 함수
# ------------------------------------------------------------
#
# 추출된 관측치에 대하여 다음 IPW 제곱오차를 최소화한다.
#
#   sum_i [delta_i / pi_i] *
#     (y_i - x_i^T beta)^2
#
# 추출된 관측치만 사용하므로 lm.wfit()에는
#
#   weights = 1 / pi_i
#
# 를 입력한다.

fit_ipw <- function(
    X,
    y,
    delta,
    pi,
    coefficient_names
) {
  
  # 추출된 관측치의 인덱스
  selected <- which(delta == 1L)
  
  
  # 추출 표본의 설계행렬
  X_sample <- X[
    selected,
    ,
    drop = FALSE
  ]
  
  
  # 추출 표본의 반응변수
  y_sample <- y[selected]
  
  
  # IPW 가중치
  w_sample <- 1 / pi[selected]
  
  
  # IPW 가중최소제곱 적합
  fit <- lm.wfit(
    x = X_sample,
    y = y_sample,
    w = w_sample
  )
  
  
  # 회귀계수 추출
  beta_hat <- fit$coefficients
  
  
  # 회귀계수 이름 지정
  names(beta_hat) <- colnames(X_sample)
  
  
  # 참 회귀계수와 동일한 순서로 정렬
  beta_hat <- beta_hat[coefficient_names]
  
  
  return(beta_hat)
}


# ------------------------------------------------------------
# 6. FULL benchmark 계산
# ------------------------------------------------------------
#
# FULL은 모집단 전체의 X와 y를 모두 사용하여
# 일반 최소제곱법으로 회귀계수를 추정한다.
#
# FULL 결과는 반복할 때마다 변하지 않으므로
# 반복문 밖에서 한 번만 계산한다.

fit_full <- lm.fit(
  x = X_pop,
  y = y_pop
)


beta_full <- fit_full$coefficients


# 추정계수 이름 지정
names(beta_full) <- colnames(X_pop)


# 참 계수와 동일한 순서로 정렬
beta_full <- beta_full[x_names]


# FULL 테스트 예측값
pred_full <- drop(
  X_test %*% beta_full
)


# FULL parameter MSE
#
# 절편을 포함한 7개 계수의 오차 제곱의 합

full_param_mse <- sum(
  (beta_full - beta_true)^2
)


# FULL prediction MSE
#
# 테스트 자료의 실제 y와 예측값 사이의 평균제곱오차

full_pred_mse <- mean(
  (y_test - pred_full)^2
)



# ------------------------------------------------------------
# 7. 반복 결과 저장용 벡터
# ------------------------------------------------------------

# 실제 추출 표본크기
opt_sample_size <- integer(n_reps)
srs_sample_size <- integer(n_reps)


# Parameter MSE
opt_param_mse <- numeric(n_reps)
srs_param_mse <- numeric(n_reps)


# Prediction MSE
opt_pred_mse <- numeric(n_reps)
srs_pred_mse <- numeric(n_reps)



# ------------------------------------------------------------
# 8. 반복 시뮬레이션
# ------------------------------------------------------------

for (rep_id in seq_len(n_reps)) {
  
  
  # ----------------------------------------------------------
  # A. 베르누이 표본추출
  # ----------------------------------------------------------
  
  delta_opt <- rbinom(
    n = N,
    size = 1,
    prob = pi_opt
  )
  
  
  delta_srs <- rbinom(
    n = N,
    size = 1,
    prob = pi_srs
  )
  
  
  # 실제 추출 표본크기 저장
  opt_sample_size[rep_id] <- sum(delta_opt)
  srs_sample_size[rep_id] <- sum(delta_srs)
  
  
  
  # ----------------------------------------------------------
  # B. IPW 회귀계수 추정
  # ----------------------------------------------------------
  
  beta_opt <- fit_ipw(
    X = X_pop,
    y = y_pop,
    delta = delta_opt,
    pi = pi_opt,
    coefficient_names = x_names
  )
  
  
  beta_srs <- fit_ipw(
    X = X_pop,
    y = y_pop,
    delta = delta_srs,
    pi = pi_srs,
    coefficient_names = x_names
  )
  
  
  
  # ----------------------------------------------------------
  # C. Parameter MSE 계산
  # ----------------------------------------------------------
  #
  # 절편을 포함한 7개 회귀계수의 제곱오차를 더한다.
  #
  # Param_MSE
  #   = mean((beta_hat - beta_true)^2)
  
  opt_param_mse[rep_id] <- sum(
    (beta_opt - beta_true)^2
  )
  
  
  srs_param_mse[rep_id] <- sum(
    (beta_srs - beta_true)^2
  )
  
  
  
  # ----------------------------------------------------------
  # D. Prediction MSE 계산
  # ----------------------------------------------------------
  
  # 테스트 자료에 대한 예측값
  pred_opt <- drop(
    X_test %*% beta_opt
  )
  
  
  pred_srs <- drop(
    X_test %*% beta_srs
  )
  
  
  # 실제 테스트 반응변수와 예측값 사이의 MSE
  opt_pred_mse[rep_id] <- mean(
    (y_test - pred_opt)^2
  )
  
  
  srs_pred_mse[rep_id] <- mean(
    (y_test - pred_srs)^2
  )
  
  
  
  # ----------------------------------------------------------
  # E. 진행 상황 출력
  # ----------------------------------------------------------
  
  if (rep_id %% 100L == 0L) {
    
    cat(
      rep_id,
      "번째 반복 완료...\n"
    )
  }
}



# ------------------------------------------------------------
# 9. 반복별 결과 데이터 생성
# ------------------------------------------------------------
#
# OPT와 SRS는 각각 n_reps개의 반복 결과를 저장한다.
#
# FULL은 반복에 따라 변하지 않는 고정된 benchmark이므로
# 한 행만 저장한다.

results <- rbind(
  
  # OPT 결과
  data.frame(
    Rep = seq_len(n_reps),
    Method = "OPT",
    Sample_Size = opt_sample_size,
    Param_MSE = opt_param_mse,
    Pred_MSE = opt_pred_mse
  ),
  
  
  # SRS 결과
  data.frame(
    Rep = seq_len(n_reps),
    Method = "SRS",
    Sample_Size = srs_sample_size,
    Param_MSE = srs_param_mse,
    Pred_MSE = srs_pred_mse
  ),
  
  
  # FULL 결과
  data.frame(
    Rep = NA_integer_,
    Method = "FULL",
    Sample_Size = N,
    Param_MSE = full_param_mse,
    Pred_MSE = full_pred_mse
  )
)


# 방법 순서 지정
results$Method <- factor(
  results$Method,
  levels = c(
    "OPT",
    "SRS",
    "FULL"
  )
)



# ------------------------------------------------------------
# 10. 방법별 평균 결과
# ------------------------------------------------------------

summary_table <- aggregate(
  cbind(
    Sample_Size,
    Param_MSE,
    Pred_MSE
  ) ~ Method,
  data = results,
  FUN = mean
)


cat("\n방법별 평균 결과\n")

print(
  summary_table
)

