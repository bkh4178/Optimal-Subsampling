# 02_exp1_k_comparison.R
# Exp1: 이분산, k 비교 실험 (uses config k_grid and n0_default)
# Metric 1 (er_leading, 이론적 목적함수) / Metric 2 (pred_mse_N, 실측 ER, primary) /
# Metric 3 (param_mse_N)
# 실행 위치: sim_final/ 폴더에서  source("scripts/02_exp1_k_comparison.R")

source("config/config_final.R")
source("R/00_utils.R")
source("R/01_dgp.R")
source("R/02_sampling_methods.R")
source("R/03_variance_estimation.R")
source("R/04_metrics.R")
source("R/05_mcse_summary.R")

ensure_dir(DIR_EXP1)

cat("====================================================\n")
cat(" Exp1: 이분산, k 비교 실험\n")
cat("====================================================\n")
cat(sprintf(" N=%d / N_test=%d / n0=%d / n_rep=%d\n",
            N, N_test, n0_default, n_rep))
cat(sprintf(" k: %s\n", paste(k_grid, collapse = ", ")))
cat(sprintf(" %s\n", dgp_sigma_label))
cat(sprintf(" methods: %s\n", paste(METHOD_ORDER, collapse = ", ")))
cat("====================================================\n\n")

dat <- generate_data_final(N, beta_true, heteroscedastic = TRUE,
                           seed = SEED_DATA)
dat_test <- generate_data_final(N_test, beta_true, heteroscedastic = TRUE,
                                seed = SEED_TEST)
X_test  <- dat_test$X
mu_test <- dat_test$mu
y_test  <- dat_test$y

cat(sprintf("데이터 생성 완료 (N=%d, sigma 범위 [%.3f, %.3f])\n\n",
            N, min(dat$sigma_vec), max(dat$sigma_vec)))

# ── IBOSS: 결정론적이므로 k당 1회만 계산 후 n_rep행으로 broadcast ────────────
# (results_to_rows 재사용 -- 컬럼 순서/이름을 다른 method와 완전히 일치시킴)
compute_iboss_once <- function(dat, k, X_test, y_test, mse0) {
  N    <- nrow(dat$X)
  beta <- dat$beta
  samp <- sample_iboss(dat, k)
  bhat <- samp$beta_hat

  if (is.null(bhat)) {
    return(list(
      method = "IBOSS", fail = TRUE,
      er_leading = NA_real_,
      pred_mse = NA_real_, pred_mse_N = NA_real_,
      param_mse = NA_real_, param_mse_N = NA_real_,
      c1_pi = NA_real_,
      n_realized = samp$n_realized, n_clip = NA_real_,
      ess_pi = NA_real_, max_weight = NA_real_, cond_XtWX = NA_real_,
      diag_sigma_cor = NA_real_, diag_score_cor = NA_real_, diag_sigma_cor_spearman = NA_real_
    ))
  }

  # IBOSS는 IPW 가중치가 없어 다른 method의 weighted cond_XtWX와 스케일이 다름
  # -> 비교 불가하므로 NA 처리 (Metric 2, 3만 적용 가능한 method)
  cond_XtWX <- NA_real_

  mse_hat       <- mean(as.numeric(y_test - X_test %*% bhat)^2)
  pred_mse_raw  <- mse_hat - mse0
  param_mse_raw <- mean((bhat - beta)^2)

  list(
    method      = "IBOSS", fail = FALSE,
    er_leading  = NA_real_,                    # IBOSS엔 확률적 π 없음 -> 정의 안 됨
    pred_mse    = pred_mse_raw,
    pred_mse_N  = N * pred_mse_raw,
    param_mse   = param_mse_raw,
    param_mse_N = N * param_mse_raw,
    c1_pi       = NA_real_,
    n_realized  = samp$n_realized, n_clip = NA_real_,
    ess_pi = NA_real_, max_weight = NA_real_, cond_XtWX = cond_XtWX,
    diag_sigma_cor = NA_real_, diag_score_cor = NA_real_, diag_sigma_cor_spearman = NA_real_
  )
}

# k당 1회 계산한 IBOSS 결과를 n_rep개 행으로 복제 (rep=1..n_rep, 값은 전부 동일)
build_iboss_rows <- function(dat, k, X_test, y_test, mse0, n_rep) {
  res_once <- compute_iboss_once(dat, k, X_test, y_test, mse0)
  row1     <- results_to_rows(list(IBOSS = res_once), rep_id = 1L, k = k)
  rows     <- row1[rep(1L, n_rep), ]
  rows$rep <- seq_len(n_rep)
  rownames(rows) <- NULL
  rows
}

# Metric 2 baseline: MSE0 = 실현된 test set의 irreducible noise (전 method/k 공통 상수)
mse0 <- mean(as.numeric(y_test - mu_test)^2)

all_rows <- list()
idx_row  <- 1L

for (k in k_grid) {
  cat(sprintf("[k = %d]\n", k))
  t_k <- proc.time()

  for (r in seq_len(n_rep)) {
    set.seed(SEED_REP + r)
    res  <- run_one_rep(dat, k, X_test, mu_test, y_test, n0 = n0_default)
    df_r <- results_to_rows(res, rep_id = r, k = k)
    all_rows[[idx_row]] <- df_r
    idx_row <- idx_row + 1L

    if (r %% 100L == 0L) cat(sprintf("  %d / %d\n", r, n_rep))
  }

  # IBOSS: rep loop 밖에서 1회만 계산 (결정론적이므로 반복 불필요)
  all_rows[[idx_row]] <- build_iboss_rows(dat, k, X_test, y_test, mse0, n_rep)
  idx_row <- idx_row + 1L

  cat(sprintf("  완료 (%.1f초)\n", (proc.time() - t_k)[["elapsed"]]))
}

results_df <- do.call(rbind, all_rows)
outfile    <- file.path(DIR_EXP1, "exp1_k_comparison.rds")
saveRDS(results_df, outfile)
cat(sprintf("\n저장 완료: %s (%d행)\n", outfile, nrow(results_df)))

# ── 결과 출력 (Metric 1/2/3) ─────────────────────────────────────────────────
print_metric_table(results_df, "er_leading",  k_grid,
                   "[1] 평균 er_leading (Metric 1, leading term / 목적함수)")
print_metric_table(results_df, "pred_mse_N",  k_grid,
                   "[2] 평균 pred_mse_N (Metric 2, 실측 ER, Primary)")
print_metric_table(results_df, "param_mse_N", k_grid,
                   "[3] 평균 param_mse_N (Metric 3, parameter MSE)")
print_metric_table(results_df, "c1_pi",       k_grid,
                   "[4] 평균 c1_pi (diagnostic)")
print_clip_table(results_df, k_grid)
print_instability_table(results_df, k_grid)
print_paired_diffs(results_df, k_grid, metric = "pred_mse_N")
print_paired_diffs(results_df, k_grid, metric = "er_leading")   # 추가

# ── Plugin 추정 품질 (param / nonparam) ──────────────────────────────────────
cat(sprintf("\n%s\n", strrep("=", 65)))
cat(" Plugin 추정 품질 (sigma_cor, score_cor)\n")
cat(sprintf("%s\n", strrep("=", 65)))
cat(sprintf("  %-30s %6s %12s %12s\n", "method", "k", "sigma_cor", "score_cor"))
cat(strrep("-", 62), "\n")
for (meth in c("OPT-hetero-plugin-param", "OPT-hetero-plugin-nonparam")) {
  plug <- results_df[results_df$method == meth & !results_df$fail, ]
  for (k_ in k_grid) {
    sub <- plug[plug$k == k_, ]
    cat(sprintf("  %-30s %6d %12.4f %12.4f\n",
                meth, k_,
                mean(sub$diag_sigma_cor, na.rm = TRUE),
                mean(sub$diag_score_cor, na.rm = TRUE)))
  }
}

# ── 자동 판정 (Metric 2: pred_mse_N 기준) ────────────────────────────────────
cat(sprintf("\n%s\n", strrep("=", 65)))
cat(" 자동 판정 (pred_mse_N 기준)\n")
cat(sprintf("%s\n", strrep("=", 65)))
valid <- results_df[!results_df$fail, ]
for (k_ in k_grid) {
  sub         <- valid[valid$k == k_, ]
  homo_ex     <- mean(sub$pred_mse_N[sub$method == "OPT-homo"],                   na.rm = TRUE)
  oracle_ex   <- mean(sub$pred_mse_N[sub$method == "OPT-hetero-oracle"],          na.rm = TRUE)
  pparam_ex   <- mean(sub$pred_mse_N[sub$method == "OPT-hetero-plugin-param"],    na.rm = TRUE)
  pnonpar_ex  <- mean(sub$pred_mse_N[sub$method == "OPT-hetero-plugin-nonparam"], na.rm = TRUE)

  cat(sprintf("k=%4d: oracle/homo=%.4f [%s]  plugin-param/homo=%.4f [%s]  plugin-nonparam/homo=%.4f [%s]\n",
              k_,
              oracle_ex / homo_ex,
              ifelse(oracle_ex < homo_ex, "gain✓", "no gain"),
              pparam_ex / homo_ex,
              ifelse(pparam_ex < homo_ex, "gain✓", "no gain"),
              pnonpar_ex / homo_ex,
              ifelse(pnonpar_ex < homo_ex, "gain✓", "no gain")))
  cat(sprintf("         plugin-param/oracle=%.4f [%s]  plugin-nonparam/oracle=%.4f [%s]\n",
              pparam_ex / oracle_ex,
              ifelse(pparam_ex < oracle_ex * 1.10, "ok✓", "large gap"),
              pnonpar_ex / oracle_ex,
              ifelse(pnonpar_ex < oracle_ex * 1.10, "ok✓", "large gap")))
}