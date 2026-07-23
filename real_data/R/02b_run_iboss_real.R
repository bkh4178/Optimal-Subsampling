# =============================================================================
# 02b_run_iboss_real.R
# NYC Yellow Taxi — IBOSS (deterministic, no repetition needed)
# 실행 위치: real_data/R/
# =============================================================================

source("../simulation_2plugin/R/00_utils.R")   # sample_iboss 등
source('../simulation_2plugin/R/02_sampling_methods.R')
source("config_real.R")

dat <- readRDS("data/taxi_prepared.rds")

X_train   <- dat$X_train
y_train   <- dat$y_train
X_test    <- dat$X_test
y_test    <- dat$y_test
beta_full <- dat$beta_full
A_hat     <- dat$A_hat
N_train   <- dat$N_train

# Metric2 baseline: 02_run_real.R과 동일하게 계산 (FULL로 mse0 근사)
test_mse_full <- mean((y_test - X_test %*% beta_full)^2)

# ⚠ sample_iboss()가 list$beta를 실제로 사용하는지 확인 필요.
#   sim에서는 true beta였지만 real data엔 true beta가 없어 beta_full을
#   placeholder로 전달함 — 내부에서 beta를 잔차 계산 등에 실제로 쓰는
#   구조라면 로직 재검토 필요.

# ── IBOSS: k당 1회 계산 후 rep 1..B로 broadcast (02_run_real.R의 results_df와
#    동일한 컬럼 스키마로 맞춰서 나중에 그대로 rbind 가능하게 함) ─────────────
compute_iboss_once <- function(k) {
  samp <- sample_iboss(list(X = X_train, y = y_train, beta = beta_full), k)
  bhat <- samp$beta_hat

  # IBOSS는 IPW weight가 없는 deterministic selection이라 clip/ess/max_weight/
  # cond_XtWX(weighted) 개념이 다른 method와 스케일이 안 맞음 -> NA 처리
  # (sim의 compute_iboss_once와 동일한 결정)
  if (is.null(bhat)) {
    return(data.frame(
      method = "IBOSS", k = k, fail = TRUE,
      er_leading = NA_real_, c1_pi = NA_real_,
      pseudo_er = NA_real_, pred_mse_N = NA_real_, param_mse_N = NA_real_,
      test_mse = NA_real_,
      n_realized = samp$n_realized, n_clip = NA_real_, clip_rate = NA_real_,
      ess_pi = NA_real_, max_weight = NA_real_, cond_XtWX = NA_real_,
      diag_sigma_cor = NA_real_, diag_score_cor = NA_real_,
      diag_sigma_cor_spearman = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  diff        <- bhat - beta_full
  pseudo_er   <- as.numeric(t(diff) %*% A_hat %*% diff)
  param_mse_N <- N_train * mean(diff^2)
  test_mse    <- mean((y_test - X_test %*% bhat)^2)
  pred_mse_N  <- N_train * (test_mse - test_mse_full)

  data.frame(
    method = "IBOSS", k = k, fail = FALSE,
    er_leading = NA_real_, c1_pi = NA_real_,
    pseudo_er = pseudo_er, pred_mse_N = pred_mse_N, param_mse_N = param_mse_N,
    test_mse = test_mse,
    n_realized = samp$n_realized, n_clip = NA_real_, clip_rate = NA_real_,
    ess_pi = NA_real_, max_weight = NA_real_, cond_XtWX = NA_real_,
    diag_sigma_cor = NA_real_, diag_score_cor = NA_real_,
    diag_sigma_cor_spearman = NA_real_,
    stringsAsFactors = FALSE
  )
}

iboss_rows <- vector("list", length(k_grid))
for (i in seq_along(k_grid)) {
  k <- k_grid[i]
  cat(sprintf("IBOSS k=%d ...\n", k))
  iboss_rows[[i]] <- compute_iboss_once(k)
}
iboss_once_df <- do.call(rbind, iboss_rows)

# rep 1..B로 broadcast (값은 전부 동일 -- deterministic이므로)
iboss_df <- do.call(rbind, lapply(seq_len(B), function(b) {
  rows <- iboss_once_df
  rows$rep <- b
  rows
}))
iboss_df <- iboss_df[, c("rep", setdiff(names(iboss_df), "rep"))]  # rep을 맨 앞으로

saveRDS(iboss_df, "../results/real_results_iboss.rds")
cat("✓ 저장 완료: ../results/real_results_iboss.rds\n")

cat("\n===== IBOSS 요약 (k당 1회 계산값) =====\n")
print(iboss_once_df[, c("method", "k", "fail", "pred_mse_N", "param_mse_N",
                          "pseudo_er", "test_mse", "n_realized")])