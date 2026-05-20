# run_s2_mild.R
# S2-mild: Moderate Variance-Leverage Alignment
#
# DGP: sigma2 = exp(0.3 * ell_score)  ← gamma=0.3 (기존 S2는 gamma=0.7)
# 목적: moderate alignment에서 hetero-aware rule이 OPT-homo와 경쟁하는지 확인
#
# 빠른 확인용: k ∈ {500, 1000} / n_rep = 100
# 결과 좋으면 k ∈ {200, 500, 1000, 2000} / n_rep = 300으로 확장

source("code3/dgp.R")
source("code3/methods.R")
source("code3/metrics.R")

# ── Sanity Diagnostic ────────────────────────────────────────────────────────
cat("========== [S2-mild] Sanity Diagnostic (gamma=0.3) ==========\n")
set.seed(42)
dat_check <- generate_data(N = 10000, case = 2, gamma = 0.3)

cat(sprintf("\ncor(sigma2, ell)       = %.4f  (S2-strong ref: ~0.7)\n",
            cor(dat_check$sigma2, dat_check$ell)))
cat(sprintf("cor(sigma2, ell_rank)  = %.4f\n",
            cor(dat_check$sigma2, rank(dat_check$ell))))
cat("\nsummary(sigma2):\n")
print(summary(dat_check$sigma2))

# ESS, concentration 비교
for (kv in c(500, 1000)) {
  pi_oracle <- compute_pi(sqrt(dat_check$sigma2 * dat_check$ell), k = kv)
  pi_homo   <- compute_pi(sqrt(dat_check$ell), k = kv)
  ess_o <- sum(pi_oracle)^2 / sum(pi_oracle^2)
  ess_h <- sum(pi_homo)^2   / sum(pi_homo^2)
  top_o <- sum(pi_oracle[pi_oracle > quantile(pi_oracle, 0.99)]) / sum(pi_oracle)
  mm_o  <- max(pi_oracle) / (min(pi_oracle[pi_oracle > 0]) + 1e-15)
  cat(sprintf("\nk=%d | oracle ESS=%.1f  homo ESS=%.1f  oracle top1%%=%.4f  oracle max/min=%.1f\n",
              kv, ess_o, ess_h, top_o, mm_o))
}

# ── 실험 설정 ─────────────────────────────────────────────────────────────────
k_vals  <- c(500, 1000)    # 빠른 확인 먼저
n_rep   <- 100
N       <- 10000
n0      <- 200
gamma   <- 0.3

results <- list()
idx     <- 1

for (k_val in k_vals) {
  for (b in 1:n_rep) {

    set.seed(b * 1000 + 20 + k_val + round(gamma * 10))

    dat      <- generate_data(N = N, case = 2, gamma = gamma)
    dat_test <- generate_data(N = N, case = 2, gamma = gamma)

    X_test <- dat_test$X
    y_test <- dat_test$y

    methods_list <- list(
      tryCatch(run_full(dat),                                             error = function(e) NULL),
      tryCatch(run_uni_ipw(dat, k = k_val),                              error = function(e) NULL),
      tryCatch(run_lev_ipw(dat, k = k_val),                              error = function(e) NULL),
      tryCatch(run_opt_homo(dat, k = k_val),                             error = function(e) NULL),
      tryCatch(run_opt_hetero_oracle(dat, k = k_val),                    error = function(e) NULL),
      tryCatch(run_opt_hetero_plugin(dat, k = k_val, n0 = n0),           error = function(e) NULL),
      tryCatch(run_opt_hetero_tempering(dat, k = k_val, gamma = 0.5),    error = function(e) NULL),
      tryCatch(run_opt_hetero_plugin_stab(dat, k = k_val, n0 = n0),      error = function(e) NULL)
    )

    for (res in methods_list) {
      if (is.null(res)) next

      m        <- compute_metrics(res, dat, X_test, y_test)
      m$case   <- "S2-mild"
      m$gamma  <- gamma
      m$k      <- k_val
      m$rep    <- b

      results[[idx]] <- m
      idx <- idx + 1
    }
  }
  cat(sprintf("S2-mild (gamma=%.1f) | k = %4d done\n", gamma, k_val))
}

sim_mild <- do.call(rbind, results)

saveRDS(sim_mild, "results/sim_s2_mild_results.rds")
write.csv(sim_mild, "results/sim_s2_mild_results.csv", row.names = FALSE)

cat("Done. Total rows:", nrow(sim_mild), "\n")
# 예상 행 수: 2 k × 100 reps × 8 methods = 1,600

# ── 터미널 확인 ───────────────────────────────────────────────────────────────

method_order <- c("FULL", "UNI-IPW", "LEV-IPW", "OPT-homo",
                  "OPT-hetero-oracle", "OPT-hetero-plugin",
                  "OPT-hetero-tempering-0.5",
                  "OPT-plugin-stab(λ=0.5,γ=0.5)")
sim_mild$method <- factor(sim_mild$method, levels = method_order)

# 1. k별 순위
cat("\n\n========== [S2-mild] 1. k별 순위 (mean excess risk, FULL 제외) ==========\n")
agg <- aggregate(excess_risk ~ k + method,
                 data = sim_mild[sim_mild$method != "FULL", ], FUN = mean)
agg <- agg[order(agg$k, agg$excess_risk), ]
for (kv in k_vals) {
  sub <- agg[agg$k == kv, ]
  cat(sprintf("\n  k = %4d\n", kv))
  for (i in seq_len(nrow(sub))) {
    cat(sprintf("    %d. %-40s %.5f\n", i, as.character(sub$method[i]), sub$excess_risk[i]))
  }
}

# 2. Tail diagnostic
cat("\n\n========== [S2-mild] 2. Tail diagnostic — mean vs. median ==========\n")
cat(sprintf("  %-40s %6s %8s %8s %8s\n", "method", "k", "mean", "median", "mean/med"))
for (kv in k_vals) {
  for (m in levels(sim_mild$method)[-1]) {
    er <- sim_mild$excess_risk[sim_mild$method == m & sim_mild$k == kv]
    if (length(er) == 0) next
    cat(sprintf("  %-40s %6d %8.5f %8.5f %8.3f\n",
                m, kv, mean(er), median(er), mean(er) / median(er)))
  }
}

# 3. π concentration diagnostic
cat("\n\n========== [S2-mild] 3. π concentration diagnostic ==========\n")
cat(sprintf("  %-40s %6s %10s %10s %10s\n",
            "method", "k", "ESS_pi", "top1%_mass", "max/min_pi"))

for (kv in k_vals) {
  set.seed(1 * 1000 + 20 + kv + round(gamma * 10))
  dat_diag <- generate_data(N = N, case = 2, gamma = gamma)

  pi_list <- list(
    "UNI-IPW"                  = rep(kv / N, N),
    "LEV-IPW"                  = compute_pi(dat_diag$ell, k = kv),
    "OPT-homo"                 = compute_pi(sqrt(dat_diag$ell), k = kv),
    "OPT-hetero-oracle"        = compute_pi(sqrt(dat_diag$sigma2 * dat_diag$ell), k = kv),
    "OPT-hetero-tempering-0.5" = compute_pi(sqrt(dat_diag$sigma2 * dat_diag$ell)^0.5, k = kv)
  )

  for (nm in names(pi_list)) {
    pv      <- pi_list[[nm]]
    ess     <- sum(pv)^2 / sum(pv^2)
    top_idx <- pv > quantile(pv, 0.99)
    top_mass <- if (sum(top_idx) == 0) 0.01 else sum(pv[top_idx]) / sum(pv)
    mm_ratio <- max(pv) / (min(pv[pv > 0]) + 1e-15)
    cat(sprintf("  %-40s %6d %10.1f %10.4f %10.1f\n",
                nm, kv, ess, top_mass, mm_ratio))
  }
}

# 4. 분포 진단
cat("\n\n========== [S2-mild] 4. 분포 진단 (excess risk) ==========\n")
for (kv in k_vals) {
  cat(sprintf("\n  k = %4d\n", kv))
  sub <- sim_mild[sim_mild$method != "FULL" & sim_mild$k == kv, ]
  cat(sprintf("  %-40s %8s %8s %8s %8s\n", "method", "mean", "median", "sd", "q90"))
  for (m in levels(sim_mild$method)[-1]) {
    er <- sub$excess_risk[sub$method == m]
    if (length(er) == 0) next
    cat(sprintf("  %-40s %8.5f %8.5f %8.5f %8.5f\n",
                m, mean(er), median(er), sd(er), quantile(er, 0.90)))
  }
}

# 5. S2-strong vs. S2-mild ESS 비교 요약
cat("\n\n========== [S2-mild vs. S2-strong] oracle ESS 비교 ==========\n")
cat("  (S2-strong ref: oracle ESS ≈ 7,973 at k=500)\n")
set.seed(42)
for (kv in c(500, 1000)) {
  dat_mild   <- generate_data(N = N, case = 2, gamma = 0.3)
  dat_strong <- generate_data(N = N, case = 2, gamma = 0.7)
  pi_mild    <- compute_pi(sqrt(dat_mild$sigma2 * dat_mild$ell), k = kv)
  pi_strong  <- compute_pi(sqrt(dat_strong$sigma2 * dat_strong$ell), k = kv)
  pi_homo    <- compute_pi(sqrt(dat_mild$ell), k = kv)
  ess_mild   <- sum(pi_mild)^2   / sum(pi_mild^2)
  ess_strong <- sum(pi_strong)^2 / sum(pi_strong^2)
  ess_homo   <- sum(pi_homo)^2   / sum(pi_homo^2)
  cat(sprintf("  k=%d | oracle-mild ESS=%.1f  oracle-strong ESS=%.1f  homo ESS=%.1f\n",
              kv, ess_mild, ess_strong, ess_homo))
}