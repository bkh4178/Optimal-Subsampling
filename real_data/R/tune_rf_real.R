# =============================================================================
# tune_rf_real.R
# Real data pilot에서 OOB error 기준 RF mtry/nodesize 1회 튜닝
# 실행 위치: real_data/R/
# =============================================================================

source("config_real.R")
library(randomForest)

dat     <- readRDS("data/taxi_prepared.rds")
X_train <- dat$X_train
y_train <- dat$y_train
N_train <- dat$N_train

# 본 실험 seed(seed*10000+b, b=1..B)와 절대 겹치지 않는 별도 seed 사용
SEED_TUNE <- 20260721

mtry_grid     <- c(1, 2, 3, 4)      # covariate 4개(intercept 제외) 기준
nodesize_grid <- c(1, 2, 3, 5)
n_seeds       <- 3                  # pilot(n0)이 작아 OOB가 흔들릴 수 있어 반복 평균

grid <- expand.grid(mtry = mtry_grid, nodesize = nodesize_grid)
grid$oob_mse <- NA_real_

for (i in seq_len(nrow(grid))) {
  oob_vals <- numeric(n_seeds)
  for (s in seq_len(n_seeds)) {
    set.seed(SEED_TUNE + i * 100L + s)
    idx_pilot <- sample(N_train, n0)

    beta_pilot <- as.vector(solve(crossprod(X_train[idx_pilot, ]),
                                   crossprod(X_train[idx_pilot, ], y_train[idx_pilot])))
    abs_resid  <- abs(as.vector(y_train[idx_pilot] - X_train[idx_pilot, ] %*% beta_pilot))

    fit <- randomForest(
      x = X_train[idx_pilot, -1, drop = FALSE],   # intercept(열1) 제외
      y = abs_resid,
      mtry = grid$mtry[i], nodesize = grid$nodesize[i]
    )
    oob_vals[s] <- tail(fit$mse, 1)
  }
  grid$oob_mse[i] <- mean(oob_vals)
  cat(sprintf("mtry=%d nodesize=%d  OOB_mse=%.4f\n",
              grid$mtry[i], grid$nodesize[i], grid$oob_mse[i]))
}

best <- grid[which.min(grid$oob_mse), ]
cat(sprintf("\nBest: mtry=%d, nodesize=%d (OOB_mse=%.4f)\n",
            best$mtry, best$nodesize, best$oob_mse))

if (!dir.exists("../results")) dir.create("../results", recursive = TRUE)
saveRDS(grid, "../results/rf_tuning_grid.rds")