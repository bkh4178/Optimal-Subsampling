# 04_exp3_runtime_comparison.R
# Exp3: 전체 pipeline 런타임 측정
# 실행 위치: sim_final/ 폴더에서  source("scripts/04_exp3_runtime_comparison.R")

source("config/config_final.R")
source("R/00_utils.R")
source("R/01_dgp.R")
source("R/02_sampling_methods.R")
source("R/03_variance_estimation.R")
source("R/04_metrics.R")
source("R/06_runtime_utils.R")

ensure_dir(DIR_EXP3)

N_BENCH  <- 50L
KS_BENCH <- c(200L, 600L, 1000L)

cat("====================================================\n")
cat(" Exp3: 런타임 측정\n")
cat("====================================================\n")
cat(sprintf(" N=%d, n0=%d, n_bench=%d\n", N, n0_default, N_BENCH))
cat(sprintf(" k: %s\n", paste(KS_BENCH, collapse = ", ")))
cat("====================================================\n\n")

dat <- generate_data_final(N, beta_true, heteroscedastic = TRUE,
                           seed = SEED_DATA)
dat_test <- generate_data_final(N_test, beta_true, heteroscedastic = TRUE,
                                seed = SEED_TEST)
X_test  <- dat_test$X
mu_test <- dat_test$mu
y_test  <- dat_test$y

runtime_rows <- list()

for (k in KS_BENCH) {
  cat(sprintf("[k = %d] %d회 벤치마크 ...\n", k, N_BENCH))
  bench <- benchmark_rep(dat, k, X_test, mu_test, y_test,
                         n0 = n0_default,
                         n_bench = N_BENCH,
                         seed_start = 3001L)
  runtime_rows[[length(runtime_rows) + 1L]] <- data.frame(
    k        = k,
    n0       = n0_default,
    n_bench  = N_BENCH,
    mean_s   = bench$mean,
    sd_s     = bench$sd,
    median_s = bench$median,
    stringsAsFactors = FALSE
  )
  cat(sprintf("  mean=%.4f초  sd=%.4f초  median=%.4f초\n\n",
              bench$mean, bench$sd, bench$median))
}

runtime_df <- do.call(rbind, runtime_rows)
outfile    <- file.path(DIR_EXP3, "exp3_runtime.rds")
saveRDS(runtime_df, outfile)
cat(sprintf("저장 완료: %s\n", outfile))

cat(sprintf("\n%s\n", strrep("=", 55)))
cat(" 런타임 요약 (초 / rep)\n")
cat(sprintf("%s\n", strrep("=", 55)))
cat(sprintf("  %-8s %10s %10s %10s\n", "k", "mean", "sd", "median"))
cat(strrep("-", 42), "\n")
for (i in seq_len(nrow(runtime_df))) {
  cat(sprintf("  %-8d %10.4f %10.4f %10.4f\n",
              runtime_df$k[i], runtime_df$mean_s[i],
              runtime_df$sd_s[i], runtime_df$median_s[i]))
}
