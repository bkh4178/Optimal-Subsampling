# run_s1.R
# Simulation Section v3 — S1: Homoscedastic Baseline
#
# DGP : case = 1, sigma2 = 1 (constant)
# Goal: OPT-homo ≈ OPT-hetero(tempering) > LEV-IPW > UNI-IPW
#
# N = 10000 / k ∈ {200, 500, 1000, 2000} / n_rep = 300 / n0 = 200

source("code2/dgp.R")
source("code2/methods.R")
source("code2/metrics.R")

# 실험 설정
k_vals <- c(200, 500, 1000, 2000)
n_rep  <- 300
N      <- 10000
n0     <- 200

results <- list()
idx     <- 1

for (k_val in k_vals) {
  for (b in 1:n_rep) {

    set.seed(b * 1000 + 10 + k_val)   # case=1 고정

    dat      <- generate_data(N = N, case = 1)
    dat_test <- generate_data(N = N, case = 1)

    X_test <- dat_test$X
    y_test <- dat_test$y

    methods_list <- list(
      tryCatch(run_full(dat),                                          error = function(e) NULL),
      tryCatch(run_uni_ipw(dat, k = k_val),                           error = function(e) NULL),
      tryCatch(run_lev_ipw(dat, k = k_val),                           error = function(e) NULL),
      tryCatch(run_opt_homo(dat, k = k_val),                          error = function(e) NULL),
      tryCatch(run_opt_hetero_oracle(dat, k = k_val),                 error = function(e) NULL),
      tryCatch(run_opt_hetero_plugin(dat, k = k_val, n0 = n0),        error = function(e) NULL),
      tryCatch(run_opt_hetero_tempering(dat, k = k_val, gamma = 0.5), error = function(e) NULL)
    )

    for (res in methods_list) {
      if (is.null(res)) next

      m        <- compute_metrics(res, dat, X_test, y_test)
      m$case   <- 1
      m$k      <- k_val
      m$rep    <- b

      results[[idx]] <- m
      idx <- idx + 1
    }
  }
  cat(sprintf("S1 | k = %4d done\n", k_val))
}

sim_s1 <- do.call(rbind, results)

saveRDS(sim_s1, "results/sim_s1_results.rds")
write.csv(sim_s1, "results/sim_s1_results.csv", row.names = FALSE)

cat("Done. Total rows:", nrow(sim_s1), "\n")
# 예상 행 수: 4 k × 300 reps × 7 methods = 8,400

# ── 터미널 1차 확인 ───────────────────────────────────────────────────────────

method_order <- c("FULL", "UNI-IPW", "LEV-IPW", "OPT-homo",
                  "OPT-hetero-oracle", "OPT-hetero-plugin",
                  "OPT-hetero-tempering-0.5")
sim_s1$method <- factor(sim_s1$method, levels = method_order)

# 1. k별 method 순위 (mean excess risk 기준, FULL 제외)
cat("\n\n========== [S1] 1. k별 순위 (mean excess risk, FULL 제외) ==========\n")
agg <- aggregate(excess_risk ~ k + method, data = sim_s1[sim_s1$method != "FULL", ],
                 FUN = mean)
agg <- agg[order(agg$k, agg$excess_risk), ]
for (kv in k_vals) {
  sub <- agg[agg$k == kv, ]
  cat(sprintf("\n  k = %4d\n", kv))
  for (i in seq_len(nrow(sub))) {
    cat(sprintf("    %d. %-35s %.5f\n", i, as.character(sub$method[i]), sub$excess_risk[i]))
  }
}

# 3. π concentration diagnostic — ESS_π, top 1% mass, max/min ratio
# clip_rate=0은 feasibility clipping 미발생을 의미하며, concentration 부재를 의미하지 않음
cat("\n\n========== [S1] 3. π concentration diagnostic ==========\n")
cat(sprintf("  %-35s %6s %8s %10s %10s %10s\n",
            "method", "k", "clip_rt", "ESS_pi", "top1%_mass", "max/min_pi"))
 
for (kv in k_vals) {
  # 해당 k, 한 replicate(b=1)에서 π 분포 계산
  set.seed(1 * 1000 + 10 + kv)
  dat_diag <- generate_data(N = N, case = 1)
 
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
    if (sum(top_idx) == 0) {
      top_mass <- 0.01  # uniform: top 1% = 정확히 1%
    } else {
      top_mass <- sum(pv[top_idx]) / sum(pv)
    }
    mm_ratio <- max(pv) / (min(pv[pv > 0]) + 1e-15)
    clip_rt  <- mean(pv >= 1 - 1e-10)
    cat(sprintf("  %-35s %6d %8.4f %10.1f %10.4f %10.1f\n",
                nm, kv, clip_rt, ess, top_mass, mm_ratio))
  }
}

# 4. 분포 진단 — mean / median / sd / q90 (FULL 제외)
cat("\n\n========== [S1] 4. 분포 진단 (excess risk) ==========\n")
for (kv in k_vals) {
  cat(sprintf("\n  k = %4d\n", kv))
  sub <- sim_s1[sim_s1$method != "FULL" & sim_s1$k == kv, ]
  cat(sprintf("  %-35s %8s %8s %8s %8s\n", "method", "mean", "median", "sd", "q90"))
  for (m in levels(sim_s1$method)[-1]) {
    er <- sub$excess_risk[sub$method == m]
    if (length(er) == 0) next
    cat(sprintf("  %-35s %8.5f %8.5f %8.5f %8.5f\n",
                m, mean(er), median(er), sd(er), quantile(er, 0.90)))
  }
}