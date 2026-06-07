# =============================================================================
# 00_eda.R
# NYC Yellow Taxi (2023-01) — EDA & Heteroscedasticity Diagnostics
# =============================================================================

library(arrow)
library(dplyr)
library(lubridate)
library(lmtest)

# =============================================================================
# 0. 데이터 로드 & 전처리
# =============================================================================

# 파일 없으면 다운로드
if (!file.exists("yellow_tripdata_2023-01.parquet")) {
  download.file(
    "https://d37ci6vzurychx.cloudfront.net/trip-data/yellow_tripdata_2023-01.parquet",
    destfile = "yellow_tripdata_2023-01.parquet",
    mode = "wb"
  )
}

df_raw <- read_parquet("yellow_tripdata_2023-01.parquet")
cat("Raw dim:", dim(df_raw), "\n")
cat("Columns:", paste(names(df_raw), collapse = ", "), "\n\n")

# 전처리
df_clean <- df_raw |>
  filter(RatecodeID == 1) |>                        # Standard rate only
  filter(fare_amount > 2.5, fare_amount < 200) |>   # 요금 범위
  filter(trip_distance > 0, trip_distance < 50) |>  # 거리 범위
  collect() |>
  mutate(
    trip_duration = as.numeric(difftime(
      tpep_dropoff_datetime, tpep_pickup_datetime, units = "mins")),
    hour = hour(tpep_pickup_datetime),
    weekday = wday(tpep_pickup_datetime, week_start = 1)  # 1=Mon, 7=Sun
  ) |>
  filter(trip_duration > 1, trip_duration < 120)

cat("Cleaned dim:", dim(df_clean), "\n")
cat("N =", nrow(df_clean), "\n\n")

# =============================================================================
# 1. 변수별 기초 통계
# =============================================================================
cat("===== fare_amount =====\n")
print(summary(df_clean$fare_amount))

cat("\n===== trip_distance =====\n")
print(summary(df_clean$trip_distance))

cat("\n===== trip_duration =====\n")
print(summary(df_clean$trip_duration))

cat("\n===== passenger_count =====\n")
print(table(df_clean$passenger_count, useNA = "ifany"))

cat("\n===== hour distribution =====\n")
print(table(df_clean$hour))

cat("\n===== weekday distribution =====\n")
print(table(df_clean$weekday))

# =============================================================================
# 2. 반응변수 & 핵심 공변량 분포
# =============================================================================
set.seed(42)
idx_plot <- sample(nrow(df_clean), 20000)

par(mfrow = c(2, 3))

# fare_amount
hist(df_clean$fare_amount[idx_plot], breaks = 60,
     main = "fare_amount", xlab = "", col = "steelblue", border = NA)

# log(fare_amount)
hist(log(df_clean$fare_amount[idx_plot]), breaks = 60,
     main = "log(fare_amount)", xlab = "", col = "steelblue", border = NA)

# trip_distance
hist(df_clean$trip_distance[idx_plot], breaks = 60,
     main = "trip_distance", xlab = "", col = "tomato", border = NA)

# log(trip_distance)
hist(log(df_clean$trip_distance[idx_plot]), breaks = 60,
     main = "log(trip_distance)", xlab = "", col = "tomato", border = NA)

# trip_duration
hist(df_clean$trip_duration[idx_plot], breaks = 60,
     main = "trip_duration", xlab = "", col = "seagreen", border = NA)

# fare vs distance scatter
plot(df_clean$trip_distance[idx_plot], df_clean$fare_amount[idx_plot],
     pch = 16, cex = 0.2, col = rgb(0, 0, 0, 0.1),
     xlab = "trip_distance", ylab = "fare_amount",
     main = "fare vs distance")

par(mfrow = c(1, 1))

# =============================================================================
# 3. 공변량 간 상관관계 (공선성 점검)
# =============================================================================
cat("\n===== Correlation matrix (distance / duration / passenger / hour) =====\n")
cor_vars <- df_clean[, c("trip_distance", "trip_duration", "passenger_count", "hour")]
print(round(cor(cor_vars, use = "complete.obs"), 3))

# =============================================================================
# 4. Full OLS 적합 (pseudo-truth)
# =============================================================================
# 공변량 구성 — hour는 선형항으로 일단 포함 (진단 후 결정)
fit_full <- lm(fare_amount ~ trip_distance + trip_duration +
                 passenger_count + hour + weekday,
               data = df_clean)

cat("\n===== Full OLS summary =====\n")
print(summary(fit_full))

# =============================================================================
# 5. 잔차 진단 — 이분산 존재 확인
# =============================================================================

fitted_vals <- fitted(fit_full)
resid_vals  <- residuals(fit_full)

# 5-1. Residuals vs Fitted (10,000개 샘플)
set.seed(42)
idx_diag <- sample(length(fitted_vals), 10000)

plot(fitted_vals[idx_diag], resid_vals[idx_diag],
     pch = 16, cex = 0.3, col = rgb(0, 0, 0, 0.2),
     xlab = "Fitted values", ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red", lwd = 1.5)

# 5-2. Breusch-Pagan test
cat("\n===== Breusch-Pagan test =====\n")
# N이 크면 p-value는 항상 유의 → 검정통계량 크기(BP) 자체가 중요
bp_result <- bptest(fit_full)
print(bp_result)

# 5-3. Fitted decile별 잔차 분산 (이분산 magnitude)
cat("\n===== Fitted decile별 잔차 분산 =====\n")
df_diag <- data.frame(
  fitted = fitted_vals,
  resid  = resid_vals
)
df_diag$decile <- cut(df_diag$fitted,
                      breaks = quantile(df_diag$fitted, probs = seq(0, 1, 0.1)),
                      include.lowest = TRUE, labels = 1:10)

decile_var <- df_diag |>
  group_by(decile) |>
  summarise(
    n       = n(),
    mean_fitted = mean(fitted),
    var_resid   = var(resid),
    sd_resid    = sd(resid)
  )
print(decile_var)

var_ratio <- max(decile_var$var_resid) / min(decile_var$var_resid)
cat(sprintf("\nMax/Min 분산 비율 (decile): %.1f\n", var_ratio))

# =============================================================================
# 6. 이분산 구조 진단 — variance vs leverage alignment
#    핵심: cor(sigma2_hat, h_ii)
# =============================================================================
cat("\n===== Variance-Leverage alignment (핵심 진단) =====\n")

# 서브샘플로 계산 (메모리·속도)
set.seed(42)
n_sub <- min(50000, nrow(df_clean))
idx_sub <- sample(nrow(df_clean), n_sub)
df_sub <- df_clean[idx_sub, ]

fit_sub <- lm(fare_amount ~ trip_distance + trip_duration +
                passenger_count + hour + weekday,
              data = df_sub)

# leverage (hat values)
h_sub <- hatvalues(fit_sub)

# sigma^2 추정: residual^2
sigma2_hat <- residuals(fit_sub)^2

# log(sigma^2) ~ X 회귀 (variance model R^2)
log_sigma2 <- log(sigma2_hat + 1e-6)
fit_var_model <- lm(log_sigma2 ~ trip_distance + trip_duration +
                      passenger_count + hour + weekday,
                    data = df_sub)
r2_var <- summary(fit_var_model)$r.squared

cat(sprintf("Pearson  cor(sigma2_hat, h_ii)      : %.4f\n", cor(sigma2_hat, h_sub)))
cat(sprintf("Spearman cor(sigma2_hat, h_ii)      : %.4f\n",
            cor(sigma2_hat, h_sub, method = "spearman")))
cat(sprintf("R² of log(sigma2) ~ X (var model)  : %.4f\n", r2_var))

# leverage 분포 요약
cat("\n===== Leverage (h_ii) 분포 =====\n")
cat(sprintf("mean: %.6f  (= p/n = %.6f)\n", mean(h_sub), ncol(model.matrix(fit_sub)) / n_sub))
cat(sprintf("max : %.6f\n", max(h_sub)))
cat(sprintf("top 1%% mass: %.4f\n",
            sum(h_sub[h_sub > quantile(h_sub, 0.99)]) / sum(h_sub)))

# =============================================================================
# 7. passenger_count 정보량 점검
# =============================================================================
cat("\n===== passenger_count vs fare (분산 작은지 확인) =====\n")
df_clean |>
  group_by(passenger_count) |>
  summarise(n = n(), mean_fare = mean(fare_amount), sd_fare = sd(fare_amount)) |>
  print()

# =============================================================================
# 8. 결과 저장
# =============================================================================
diag_summary <- list(
  N             = nrow(df_clean),
  bp_stat       = as.numeric(bp_result$statistic),
  var_ratio     = var_ratio,
  cor_sigma2_h  = cor(sigma2_hat, h_sub),
  spearman_sigma2_h = cor(sigma2_hat, h_sub, method = "spearman"),
  r2_var_model  = r2_var,
  max_leverage  = max(h_sub),
  top1pct_mass  = sum(h_sub[h_sub > quantile(h_sub, 0.99)]) / sum(h_sub)
)

cat("\n===== 진단 요약 =====\n")
print(diag_summary)

# 저장
if (!dir.exists("real_data/results")) dir.create("real_data/results", recursive = TRUE)
saveRDS(diag_summary, "real_data/results/eda_diagnostics.rds")
saveRDS(df_clean,     "real_data/data/taxi_clean.rds")

cat("\n✓ 저장 완료: real_data/results/eda_diagnostics.rds\n")
cat("✓ 저장 완료: real_data/data/taxi_clean.rds\n")
