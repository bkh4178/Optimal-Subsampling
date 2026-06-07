# smoke_test.R
# run_one_rep() 1회 실행 후 excess_risk / pred_mse / param_mse 가
# 모든 method 에서 NA 없이 숫자로 반환되는지 확인
# 실행 위치: sim_final/ 폴더에서  source("scripts/smoke_test.R")

source("config/config_final.R")
source("R/00_utils.R")
source("R/01_dgp.R")
source("R/02_sampling_methods.R")
source("R/03_variance_estimation.R")
source("R/04_metrics.R")

dat      <- generate_data_final(1000L, beta_true, heteroscedastic = TRUE,  seed = 1L)
dat_test <- generate_data_final(1000L, beta_true, heteroscedastic = TRUE,  seed = 2L)

set.seed(42L)
res <- run_one_rep(dat, k = 200L,
                   X_test  = dat_test$X,
                   mu_test = dat_test$mu,
                   y_test  = dat_test$y,
                   n0      = 50L)

ok <- TRUE
for (m in names(res)) {
  rv <- res[[m]]
  if (isTRUE(rv$fail)) {
    cat(sprintf("FAIL  [%s] fail=TRUE\n", m)); ok <- FALSE; next
  }
  for (col in c("excess_risk", "pred_mse", "param_mse")) {
    if (!is.finite(rv[[col]])) {
      cat(sprintf("FAIL  [%s] %s = %s\n", m, col, rv[[col]])); ok <- FALSE
    }
  }
}

if (ok) {
  cat("PASS  run_one_rep: all metrics finite\n\n")
  cat(sprintf("  %-25s %12s %12s %12s\n", "method", "excess_risk", "pred_mse", "param_mse"))
  cat(strrep("-", 65), "\n")
  for (m in names(res)) {
    rv <- res[[m]]
    cat(sprintf("  %-25s %12.6f %12.6f %12.6f\n",
                m, rv$excess_risk, rv$pred_mse, rv$param_mse))
  }
} else {
  cat("FAIL  smoke_test: 위 오류를 확인하세요.\n")
}
