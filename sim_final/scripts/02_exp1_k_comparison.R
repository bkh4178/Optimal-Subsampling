# 02_exp1_k_comparison.R
# Exp1: 이분산, k 비교 실험 (uses config k_grid and n0_default)
# Primary metric: excess risk
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

  cat(sprintf("  완료 (%.1f초)\n", (proc.time() - t_k)[["elapsed"]]))
}

results_df <- do.call(rbind, all_rows)
outfile    <- file.path(DIR_EXP1, "exp1_k_comparison.rds")
saveRDS(results_df, outfile)
cat(sprintf("\n저장 완료: %s (%d행)\n", outfile, nrow(results_df)))

# ── 결과 출력 ─────────────────────────────────────────────────────────────────
print_metric_table(results_df, "excess_risk", k_grid,
                   "[1] 평균 excess_risk (Primary)")
print_metric_table(results_df, "param_mse",   k_grid,
                   "[2] 평균 param_mse")
print_metric_table(results_df, "pred_mse",    k_grid,
                   "[3] 평균 pred_mse")
print_metric_table(results_df, "c1_pi",       k_grid,
                   "[4] 평균 c1_pi (diagnostic)")
print_clip_table(results_df, k_grid)
print_instability_table(results_df, k_grid)
print_paired_diffs(results_df, k_grid, metric = "excess_risk")

# ── Plugin 추정 품질 ──────────────────────────────────────────────────────────
cat(sprintf("\n%s\n", strrep("=", 65)))
cat(" Plugin 추정 품질 (sigma_cor, score_cor)\n")
cat(sprintf("%s\n", strrep("=", 65)))
cat(sprintf("  %-6s %12s %12s\n", "k", "sigma_cor", "score_cor"))
cat(strrep("-", 34), "\n")
plug <- results_df[results_df$method == "OPT-hetero-plugin" & !results_df$fail, ]
for (k_ in k_grid) {
  sub <- plug[plug$k == k_, ]
  cat(sprintf("  %-6d %12.4f %12.4f\n",
              k_,
              mean(sub$diag_sigma_cor, na.rm = TRUE),
              mean(sub$diag_score_cor, na.rm = TRUE)))
}

# ── 자동 판정 ─────────────────────────────────────────────────────────────────
cat(sprintf("\n%s\n", strrep("=", 65)))
cat(" 자동 판정 (excess_risk 기준)\n")
cat(sprintf("%s\n", strrep("=", 65)))
valid <- results_df[!results_df$fail, ]
for (k_ in k_grid) {
  sub       <- valid[valid$k == k_, ]
  homo_ex   <- mean(sub$excess_risk[sub$method == "OPT-homo"],          na.rm = TRUE)
  oracle_ex <- mean(sub$excess_risk[sub$method == "OPT-hetero-oracle"], na.rm = TRUE)
  plugin_ex <- mean(sub$excess_risk[sub$method == "OPT-hetero-plugin"], na.rm = TRUE)
  cat(sprintf("k=%4d: oracle/homo=%.4f [%s]  plugin/homo=%.4f [%s]  plugin/oracle=%.4f [%s]\n",
              k_,
              oracle_ex / homo_ex,
              ifelse(oracle_ex < homo_ex, "gain✓", "no gain"),
              plugin_ex / homo_ex,
              ifelse(plugin_ex < homo_ex, "gain✓", "no gain"),
              plugin_ex / oracle_ex,
              ifelse(plugin_ex < oracle_ex * 1.10, "ok✓", "large gap")))
}
