# 00a_tune_rf_hyperparams.R (재작성)
# 목적함수를 OOB MSE → true sigma2와의 상관계수로 교체.
# n0=200 기준, 여러 pilot draw 평균.

source("config/config_final.R")
source("R/00_utils.R")
source("R/01_dgp.R")

n0_tune   <- n0_default
n_draws   <- 50
tune_grid <- expand.grid(mtry = c(1L, 2L, 3L, 4L),
                          nodesize = c(1L, 2L, 3L, 5L))  # 50 근처는 제외, 더 작은 값 위주로 재탐색

dat_tune <- generate_data_final(N, beta_true, heteroscedastic = TRUE, seed = SEED_DATA)

cor_mat <- matrix(NA, nrow = n_draws, ncol = nrow(tune_grid))

set.seed(202607)
for (d in seq_len(n_draws)) {
  idx_pilot <- sample.int(nrow(dat_tune$X), n0_tune, replace = FALSE)
  X_pilot   <- dat_tune$X[idx_pilot, , drop = FALSE]

  fit_pilot <- lm(dat_tune$y[idx_pilot] ~ X_pilot - 1)
  res_sq    <- abs(residuals(fit_pilot))
  pilot_df  <- data.frame(res_sq = res_sq, X_pilot[, -1, drop = FALSE])

  Z_full <- as.data.frame(dat_tune$X[, -1, drop = FALSE])

  for (g in seq_len(nrow(tune_grid))) {
    rf <- randomForest(res_sq ~ ., data = pilot_df,
                        mtry     = tune_grid$mtry[g],
                        nodesize = tune_grid$nodesize[g])
    pred_full <- predict(rf, newdata = Z_full)
    cor_mat[d, g] <- cor(pred_full, dat_tune$sigma2_vec)  # ← 핵심: 목적함수를 상관계수로
  }
}

avg_cor <- colMeans(cor_mat)
best    <- tune_grid[which.max(avg_cor), ]   # 이제는 max (상관을 최대화)

cat("Best mtry:", best$mtry, " nodesize:", best$nodesize, "\n")
print(cbind(tune_grid, avg_cor = avg_cor))