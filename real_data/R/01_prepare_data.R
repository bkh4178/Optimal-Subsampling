# =============================================================================
# 01_prepare_data.R
# NYC Yellow Taxi (2023-01) — 전처리 + train/test split + 표준화
# =============================================================================

library(arrow)
library(dplyr)
library(lubridate)

# =============================================================================
# 0. 데이터 로드
# =============================================================================
if (!file.exists("yellow_tripdata_2023-01.parquet")) {
  download.file(
    "https://d37ci6vzurychx.cloudfront.net/trip-data/yellow_tripdata_2023-01.parquet",
    destfile = "yellow_tripdata_2023-01.parquet",
    mode = "wb"
  )
}

df_raw <- read_parquet("yellow_tripdata_2023-01.parquet")

# =============================================================================
# 1. 전처리
# =============================================================================
df_clean <- df_raw |>
  filter(RatecodeID == 1) |>                         # Standard rate only
  filter(fare_amount > 2.5, fare_amount < 200) |>    # 요금 범위
  filter(trip_distance > 0, trip_distance < 50) |>   # 거리 범위
  filter(passenger_count >= 1) |>                    # 입력 오류(0) 제거
  collect() |>
  mutate(
    trip_duration = as.numeric(difftime(
      tpep_dropoff_datetime, tpep_pickup_datetime, units = "mins")),
    hour    = hour(tpep_pickup_datetime),
    weekday = wday(tpep_pickup_datetime, week_start = 1)  # 1=Mon, 7=Sun
  ) |>
  filter(trip_duration > 1, trip_duration < 120)

cat("전처리 후 N =", nrow(df_clean), "\n")

# =============================================================================
# 2. 변수 구성
#    y: fare_amount
#    X: trip_distance, trip_duration, hour, weekday  (p = 4)
# =============================================================================
y <- df_clean$fare_amount
X <- df_clean[, c("trip_distance", "trip_duration", "hour", "weekday")]
X <- as.matrix(X)

cat("y dim:", length(y), "\n")
cat("X dim:", dim(X), "\n")

# =============================================================================
# 3. Train / Test split  (8:2)
# =============================================================================
set.seed(42)
n      <- nrow(X)
n_train <- round(n * 0.8)
idx_train <- sample(n, n_train)
idx_test  <- setdiff(seq_len(n), idx_train)

X_train <- X[idx_train, ]
X_test  <- X[idx_test,  ]
y_train <- y[idx_train]
y_test  <- y[idx_test]

cat("Train N:", length(y_train), "\n")
cat("Test  N:", length(y_test),  "\n")

# =============================================================================
# 4. 표준화  (train mean/sd 기준으로 train, test 동일 적용)
# =============================================================================
mu  <- colMeans(X_train)
sig <- apply(X_train, 2, sd)

X_train_std <- scale(X_train, center = mu, scale = sig)
X_test_std  <- scale(X_test,  center = mu, scale = sig)

# intercept 포함 설계행렬
X_train_mat <- cbind(1, X_train_std)
X_test_mat  <- cbind(1, X_test_std)

cat("\n표준화 파라미터 (train 기준):\n")
print(round(rbind(mean = mu, sd = sig), 4))

# =============================================================================
# 5. Full OLS (pseudo-truth)  — train 전체로 적합
# =============================================================================
beta_full <- as.vector(solve(crossprod(X_train_mat),
                             crossprod(X_train_mat, y_train)))
names(beta_full) <- c("intercept",
                      "trip_distance", "trip_duration", "hour", "weekday")

cat("\nbeta_full (pseudo-truth):\n")
print(round(beta_full, 4))

# test MSE 확인 (benchmark)
y_pred_full <- X_test_mat %*% beta_full
mse_full    <- mean((y_test - y_pred_full)^2)
cat(sprintf("\nFULL test MSE: %.4f\n", mse_full))

# =============================================================================
# 6. A_hat 계산  (excess risk 계산용)
#    A_hat = (1/N_train) * X_train^T X_train
# =============================================================================
N_train <- nrow(X_train_mat)
A_hat   <- crossprod(X_train_mat) / N_train

# =============================================================================
# 7. 저장
# =============================================================================
if (!dir.exists("real_data/data"))    dir.create("real_data/data",    recursive = TRUE)
if (!dir.exists("real_data/results")) dir.create("real_data/results", recursive = TRUE)

saveRDS(
  list(
    X_train     = X_train_mat,
    X_test      = X_test_mat,
    y_train     = y_train,
    y_test      = y_test,
    beta_full   = beta_full,
    A_hat       = A_hat,
    mu          = mu,
    sig         = sig,
    N_train     = N_train,
    mse_full    = mse_full
  ),
  "data/taxi_prepared.rds"
)

cat("\n✓ 저장 완료: real_data/data/taxi_prepared.rds\n")