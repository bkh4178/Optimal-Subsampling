# 01_exp0_homoscedastic_sanity.R
# Exp0: 등분산 sanity check
#
# 목적:
#   등분산에서 OPT-hetero-oracle 의 score = sqrt(1 * ell) = OPT-homo score
#   → score(=pi_raw)가 완전히 동일하므로 metric1(er_leading, beta_hat 불필요)은
#     rep마다 거의 정확히 일치해야 함 (부동소수점 오차 수준)
#   → metric2/3(pred_mse_N, param_mse_N)은 beta_hat 기반이라 표본추출 noise가
#     남아있으므로 500 rep 평균 내야 비율이 1에 가까워짐
#   → plugin(param/nonparam)은 sigma^2 ≈ 1 로 추정될 것이므로 homo와 유사할 것
#
# 실행 위치: sim_final/ 폴더에서  source("scripts/01_exp0_homoscedastic_sanity.R")

source("config/config_final.R")
source("R/00_utils.R")
source("R/01_dgp.R")
source("R/02_sampling_methods.R")
source("R/03_variance_estimation.R")
source("R/04_metrics.R")
source("R/05_mcse_summary.R")

ensure_dir(DIR_EXP0)

N_REP_EXP0 <- 500L
KS_EXP0    <- c(400L, 600L, 800L, 1000L, 2000L)

source("config/config_final.R")
source("R/00_utils.R")
source("R/01_dgp.R")
source("R/02_sampling_methods.R")

dat_h <- generate_data_final(N, beta_true, heteroscedastic = FALSE, seed = SEED_DATA)

scores_h <- get_scores(dat_h, dat_h$sigma2_vec, dat_h$sigma2_vec)  # dummy plugin, 여기선 안 씀
pi_homo   <- compute_pi(scores_h[["OPT-homo"]], 400)
pi_oracle <- compute_pi(scores_h[["OPT-hetero-oracle"]], 400)

cat("identical:", identical(pi_homo, pi_oracle),
    "  max diff:", max(abs(pi_homo - pi_oracle)), "\n")

cat("====================================================\n")
cat(" Exp0: 등분산 sanity check\n")
cat("====================================================\n")
cat(sprintf(" N=%d / N_test=%d / n0=%d / n_rep=%d\n",
            N, N_test, n0_default, N_REP_EXP0))
cat(sprintf(" k: %s\n", paste(KS_EXP0, collapse = ", ")))
cat(" 등분산: sigma(Z) = 1\n")
cat(" 예상: OPT-hetero-oracle ≈ OPT-homo (모든 metric), plugin ≈ homo\n")
cat("====================================================\n\n")

dat <- generate_data_final(N, beta_true, heteroscedastic = FALSE,
                           seed = SEED_DATA)
dat_test <- generate_data_final(N_test, beta_true, heteroscedastic = FALSE,
                                seed = SEED_TEST)
X_test  <- dat_test$X
mu_test <- dat_test$mu
y_test  <- dat_test$y

all_rows <- list()
idx_row  <- 1L

for (k in KS_EXP0) {
  cat(sprintf("[k = %d]\n", k))
  t_k <- proc.time()

  for (r in seq_len(N_REP_EXP0)) {
    set.seed(SEED_REP + r)
    res  <- run_one_rep(dat, k, X_test, mu_test, y_test, n0 = n0_default)
    df_r <- results_to_rows(res, rep_id = r, k = k)
    all_rows[[idx_row]] <- df_r
    idx_row <- idx_row + 1L

    if (r %% 50L == 0L) cat(sprintf("  %d / %d\n", r, N_REP_EXP0))
  }

  cat(sprintf("  완료 (%.1f초)\n", (proc.time() - t_k)[["elapsed"]]))
}

results_df <- do.call(rbind, all_rows)
outfile    <- file.path(DIR_EXP0, "exp0_homo_sanity.rds")
saveRDS(results_df, outfile)
cat(sprintf("\n저장 완료: %s (%d행)\n", outfile, nrow(results_df)))

# ── 결과 요약 (Metric 1/2/3) ──
print_metric_table(results_df, "er_leading",  KS_EXP0,
                   "[1] 평균 er_leading (Metric 1, leading term)")
print_metric_table(results_df, "pred_mse_N",  KS_EXP0,
                   "[2] 평균 pred_mse_N (Metric 2, 실측 ER)")
print_metric_table(results_df, "param_mse_N", KS_EXP0,
                   "[3] 평균 param_mse_N (Metric 3, parameter MSE)")
print_clip_table(results_df, KS_EXP0)

# ── Sanity check 1: oracle vs homo, metric1 (rep별 정밀 일치 여부) ───────────
cat(sprintf("\n%s\n", strrep("=", 70)))
cat(" Sanity check 1: oracle vs homo -- metric1(er_leading), rep별 최대 절대오차\n")
cat(" (score가 이론적으로 동일하므로 거의 0에 가까워야 함)\n")
cat(sprintf("%s\n", strrep("=", 70)))
valid <- results_df[!results_df$fail, ]
for (k_ in KS_EXP0) {
  sub    <- valid[valid$k == k_, ]
  homo_v <- sub$er_leading[sub$method == "OPT-homo"]
  orac_v <- sub$er_leading[sub$method == "OPT-hetero-oracle"]
  mg     <- merge(data.frame(rep = sub$rep[sub$method == "OPT-homo"], homo = homo_v),
                  data.frame(rep = sub$rep[sub$method == "OPT-hetero-oracle"], oracle = orac_v),
                  by = "rep")
  max_abs_diff <- max(abs(mg$homo - mg$oracle), na.rm = TRUE)
  cat(sprintf("  k=%4d: max|homo - oracle| = %.3e  [%s]\n",
              k_, max_abs_diff,
              ifelse(max_abs_diff < 1e-8, "ok✓ (사실상 동일)", "CHECK!")))
}

# ── Sanity check 2: oracle/homo/plugin 평균 비율 (metric2, 3) ───────────────
cat(sprintf("\n%s\n", strrep("=", 70)))
cat(" Sanity check 2: 방법론 간 평균 비율 (≈ 1.00 이어야 함, MC noise 존재)\n")
cat(sprintf("%s\n", strrep("=", 70)))
for (k_ in KS_EXP0) {
  sub <- valid[valid$k == k_, ]

  homo_pred   <- mean(sub$pred_mse_N[sub$method == "OPT-homo"],                   na.rm = TRUE)
  oracle_pred <- mean(sub$pred_mse_N[sub$method == "OPT-hetero-oracle"],          na.rm = TRUE)
  pparam_pred <- mean(sub$pred_mse_N[sub$method == "OPT-hetero-plugin-param"],    na.rm = TRUE)
  pnonp_pred  <- mean(sub$pred_mse_N[sub$method == "OPT-hetero-plugin-nonparam"], na.rm = TRUE)

  cat(sprintf("  k=%4d [pred_mse_N]  oracle/homo=%.4f  plugin-param/homo=%.4f  plugin-nonparam/homo=%.4f\n",
              k_, oracle_pred / homo_pred, pparam_pred / homo_pred, pnonp_pred / homo_pred))

  homo_par   <- mean(sub$param_mse_N[sub$method == "OPT-homo"],                   na.rm = TRUE)
  oracle_par <- mean(sub$param_mse_N[sub$method == "OPT-hetero-oracle"],          na.rm = TRUE)
  pparam_par <- mean(sub$param_mse_N[sub$method == "OPT-hetero-plugin-param"],    na.rm = TRUE)
  pnonp_par  <- mean(sub$param_mse_N[sub$method == "OPT-hetero-plugin-nonparam"], na.rm = TRUE)

  cat(sprintf("           [param_mse_N] oracle/homo=%.4f  plugin-param/homo=%.4f  plugin-nonparam/homo=%.4f\n",
              oracle_par / homo_par, pparam_par / homo_par, pnonp_par / homo_par))
}