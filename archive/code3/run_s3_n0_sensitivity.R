# run_s3_n0_sensitivity.R
# S3 v2 — Pilot size (n0) sensitivity
# n0 ∈ {200, 500, 1000} / k ∈ {500, 1000} / n_rep = 300
#
# 목적: pilot size가 커질수록 plugin-stab이 oracle 방향으로 수렴하는지 확인
# "pilot이 충분하면 plugin이 유효하다"는 claim 지지 또는 반박

source("code3/dgp.R")
source("code3/methods.R")
source("code3/metrics.R")

k_vals  <- c(500, 1000)     # S3에서 oracle이 비교적 경쟁력 있는 구간
n0_vals <- c(200, 500, 1000)
n_rep   <- 300
N       <- 10000

results <- list()
idx     <- 1

for (n0_val in n0_vals) {
  for (k_val in k_vals) {
    for (b in 1:n_rep) {

      set.seed(b * 1000 + 30 + k_val + n0_val)

      dat      <- generate_data(N = N, case = 3)
      dat_test <- generate_data(N = N, case = 3)

      X_test <- dat_test$X
      y_test <- dat_test$y

      methods_list <- list(
        tryCatch(run_full(dat),                                               error = function(e) NULL),
        tryCatch(run_opt_homo(dat, k = k_val),                               error = function(e) NULL),
        tryCatch(run_opt_hetero_oracle(dat, k = k_val),                      error = function(e) NULL),
        tryCatch(run_opt_hetero_plugin(dat, k = k_val, n0 = n0_val),         error = function(e) NULL),
        tryCatch(run_opt_hetero_plugin_stab(dat, k = k_val, n0 = n0_val),    error = function(e) NULL)
      )

      for (res in methods_list) {
        if (is.null(res)) next

        m        <- compute_metrics(res, dat, X_test, y_test)
        m$case   <- 3
        m$k      <- k_val
        m$n0     <- n0_val
        m$rep    <- b

        results[[idx]] <- m
        idx <- idx + 1
      }
    }
    cat(sprintf("S3-n0=%4d | k = %4d done\n", n0_val, k_val))
  }
}

sim_n0 <- do.call(rbind, results)

saveRDS(sim_n0, "results/sim_s3_n0_sensitivity.rds")
write.csv(sim_n0, "results/sim_s3_n0_sensitivity.csv", row.names = FALSE)

cat("Done. Total rows:", nrow(sim_n0), "\n")
# 예상 행 수: 3 n0 × 2 k × 300 reps × 5 methods = 9,000

# ── 터미널 확인 ───────────────────────────────────────────────────────────────

method_order <- c("FULL", "OPT-homo", "OPT-hetero-oracle",
                  "OPT-hetero-plugin",
                  "OPT-plugin-stab(λ=0.5,γ=0.5)")
sim_n0$method <- factor(sim_n0$method, levels = method_order)

cat("\n\n========== [S3 n0 sensitivity] mean excess risk ==========\n")
cat(sprintf("  %-40s %6s %6s %10s\n", "method", "k", "n0", "mean_er"))

agg <- aggregate(excess_risk ~ k + n0 + method,
                 data = sim_n0[sim_n0$method != "FULL", ], FUN = mean)
agg <- agg[order(agg$k, agg$n0, agg$excess_risk), ]

for (kv in k_vals) {
  cat(sprintf("\n  --- k = %d ---\n", kv))
  for (n0v in n0_vals) {
    cat(sprintf("  n0 = %d\n", n0v))
    sub <- agg[agg$k == kv & agg$n0 == n0v, ]
    for (i in seq_len(nrow(sub))) {
      cat(sprintf("    %-40s %10.5f\n",
                  as.character(sub$method[i]), sub$excess_risk[i]))
    }
  }
}

# plugin-stab vs. oracle gap
cat("\n\n========== [S3 n0 sensitivity] plugin-stab vs. oracle gap ==========\n")
cat(sprintf("  %6s %6s %12s %12s %12s\n", "k", "n0", "oracle", "stab", "stab/oracle"))

for (kv in k_vals) {
  for (n0v in n0_vals) {
    er_oracle <- mean(sim_n0$excess_risk[
      sim_n0$method == "OPT-hetero-oracle" &
      sim_n0$k == kv & sim_n0$n0 == n0v])
    er_stab <- mean(sim_n0$excess_risk[
      sim_n0$method == "OPT-plugin-stab(λ=0.5,γ=0.5)" &
      sim_n0$k == kv & sim_n0$n0 == n0v])
    cat(sprintf("  %6d %6d %12.5f %12.5f %12.3f\n",
                kv, n0v, er_oracle, er_stab, er_stab / er_oracle))
  }
}