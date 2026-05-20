# 01_exp0_homoscedastic_sanity.R
# Exp0: 등분산 sanity check
#
# 목적:
#   등분산에서 OPT-hetero-oracle 의 score = sqrt(1 * ell) = OPT-homo score
#   → oracle 과 homo 결과가 거의 동일해야 함
#   → plugin 은 sigma^2 ≈ 1 추정 시 homo 와 유사
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
KS_EXP0    <- c(400L, 600L, 800L)

cat("====================================================\n")
cat(" Exp0: 등분산 sanity check\n")
cat("====================================================\n")
cat(sprintf(" N=%d / N_test=%d / n0=%d / n_rep=%d\n",
            N, N_test, n0_default, N_REP_EXP0))
cat(sprintf(" k: %s\n", paste(KS_EXP0, collapse = ", ")))
cat(" 등분산: sigma(Z) = 1\n")
cat(" 예상: OPT-hetero-oracle ≈ OPT-homo (excess_risk)\n")
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

# ── 결과 요약 ──
print_metric_table(results_df, "excess_risk", KS_EXP0,
                   "[1] 평균 excess_risk (Primary)")
print_metric_table(results_df, "param_mse",   KS_EXP0,
                   "[2] 평균 param_mse")
print_clip_table(results_df, KS_EXP0)

# ── Sanity check ──
cat(sprintf("\n%s\n", strrep("=", 65)))
cat(" Sanity check: oracle/homo 비율 (≈ 1.00 이어야 함)\n")
cat(sprintf("%s\n", strrep("=", 65)))
valid <- results_df[!results_df$fail, ]
for (k_ in KS_EXP0) {
  sub       <- valid[valid$k == k_, ]
  homo_ex   <- mean(sub$excess_risk[sub$method == "OPT-homo"],          na.rm = TRUE)
  oracle_ex <- mean(sub$excess_risk[sub$method == "OPT-hetero-oracle"], na.rm = TRUE)
  plugin_ex <- mean(sub$excess_risk[sub$method == "OPT-hetero-plugin"], na.rm = TRUE)
  cat(sprintf("  k=%d: oracle/homo=%.4f [%s]  plugin/homo=%.4f\n",
              k_, oracle_ex / homo_ex,
              ifelse(abs(oracle_ex / homo_ex - 1) < 0.05, "ok✓", "CHECK!"),
              plugin_ex / homo_ex))
}
