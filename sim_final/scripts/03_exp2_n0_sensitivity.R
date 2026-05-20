# 03_exp2_n0_sensitivity.R
# Exp2: n0 sensitivity (uses config k_exp2 and n0_grid)
# 목적: pilot 크기가 plugin 성능에 미치는 영향
# 실행 위치: sim_final/ 폴더에서  source("scripts/03_exp2_n0_sensitivity.R")

source("config/config_final.R")
source("R/00_utils.R")
source("R/01_dgp.R")
source("R/02_sampling_methods.R")
source("R/03_variance_estimation.R")
source("R/04_metrics.R")
source("R/05_mcse_summary.R")

ensure_dir(DIR_EXP2)

cat("====================================================\n")
cat(" Exp2: n0 sensitivity\n")
cat("====================================================\n")
cat(sprintf(" N=%d / N_test=%d / k=%d / n_rep=%d\n",
            N, N_test, k_exp2, n_rep))
cat(sprintf(" n0_grid: %s\n", paste(n0_grid, collapse = ", ")))
cat(sprintf(" %s\n", dgp_sigma_label))
cat("====================================================\n\n")

dat <- generate_data_final(N, beta_true, heteroscedastic = TRUE,
                           seed = SEED_DATA)
dat_test <- generate_data_final(N_test, beta_true, heteroscedastic = TRUE,
                                seed = SEED_TEST)
X_test  <- dat_test$X
mu_test <- dat_test$mu
y_test  <- dat_test$y

all_rows <- list()
idx_row  <- 1L

for (n0_val in n0_grid) {
  cat(sprintf("[n0 = %d]\n", n0_val))
  t_n0 <- proc.time()

  for (r in seq_len(n_rep)) {
    set.seed(SEED_REP + r)
    res  <- run_one_rep(dat, k_exp2, X_test, mu_test, y_test, n0 = n0_val)
    df_r <- results_to_rows(res, rep_id = r, k = k_exp2)
    df_r$n0 <- n0_val
    all_rows[[idx_row]] <- df_r
    idx_row <- idx_row + 1L

    if (r %% 100L == 0L) cat(sprintf("  %d / %d\n", r, n_rep))
  }

  cat(sprintf("  완료 (%.1f초)\n", (proc.time() - t_n0)[["elapsed"]]))
}

results_df <- do.call(rbind, all_rows)
outfile    <- file.path(DIR_EXP2, "exp2_n0_sensitivity.rds")
saveRDS(results_df, outfile)
cat(sprintf("\n저장 완료: %s (%d행)\n", outfile, nrow(results_df)))

# ── n0별 Plugin 성능 요약 ──────────────────────────────────────────────────────
cat(sprintf("\n%s\n", strrep("=", 65)))
cat(sprintf(" n0별 OPT-hetero-plugin (k=%d, 이분산)\n", k_exp2))
cat(sprintf("%s\n", strrep("=", 65)))
cat(sprintf("  %-8s %12s %12s %12s %12s\n",
            "n0", "excess_risk", "mcse", "sigma_cor", "score_cor"))
cat(strrep("-", 60), "\n")

valid  <- results_df[!results_df$fail, ]
plugin <- valid[valid$method == "OPT-hetero-plugin", ]
for (n0_val in n0_grid) {
  sub <- plugin[plugin$n0 == n0_val, ]
  cat(sprintf("  %-8d %12.6f %12.6f %12.4f %12.4f\n",
              n0_val,
              mean(sub$excess_risk,    na.rm = TRUE),
              compute_mcse(sub$excess_risk),
              mean(sub$diag_sigma_cor, na.rm = TRUE),
              mean(sub$diag_score_cor, na.rm = TRUE)))
}

# ── n0별 전 method excess_risk ─────────────────────────────────────────────────
cat(sprintf("\n%s\n", strrep("=", 65)))
cat(sprintf(" n0별 excess_risk (k=%d)\n", k_exp2))
cat(sprintf("%s\n", strrep("=", 65)))
cat(sprintf("%-25s", "method \\ n0"))
for (n0_val in n0_grid) cat(sprintf("%12d", n0_val))
cat("\n")
cat(strrep("-", 25 + 12 * length(n0_grid)), "\n")
for (m in METHOD_ORDER) {
  cat(sprintf("%-25s", m))
  for (n0_val in n0_grid) {
    sub <- valid[valid$method == m & valid$n0 == n0_val, ]
    if (nrow(sub) == 0L) cat(sprintf("%12s", "NA"))
    else cat(sprintf("%12.5f", mean(sub$excess_risk, na.rm = TRUE)))
  }
  cat("\n")
}

# ── Paired diff (3 contrasts) with lower/upper ────────────────────────────────
contrasts <- list(
  list(m1 = "SRS",               m2 = "OPT-hetero-plugin",  label = "SRS - plugin"),
  list(m1 = "OPT-homo",          m2 = "OPT-hetero-plugin",  label = "homo - plugin"),
  list(m1 = "OPT-hetero-plugin", m2 = "OPT-hetero-oracle",  label = "plugin - oracle")
)

for (ct in contrasts) {
  cat(sprintf("\n%s\n", strrep("=", 65)))
  cat(sprintf(" n0별 페어 차이: %s  (excess_risk 기준)\n", ct$label))
  cat(sprintf("%s\n", strrep("=", 65)))
  cat(sprintf("  %-8s %12s %12s %12s %12s\n",
              "n0", "mean_diff", "mcse_diff", "lower", "upper"))
  cat(strrep("-", 60), "\n")

  for (n0_val in n0_grid) {
    sub <- results_df[results_df$n0 == n0_val, ]
    pd  <- tryCatch(
      paired_diff_summary(sub, ct$m1, ct$m2, "excess_risk"),
      error = function(e) NULL
    )
    if (!is.null(pd) && nrow(pd) > 0L) {
      cat(sprintf("  %-8d %12.6f %12.6f %12.6f %12.6f\n",
                  n0_val, pd$mean_diff[1L], pd$mcse_diff[1L],
                  pd$lower[1L], pd$upper[1L]))
    }
  }
}
