# run_s3_v2.R
# Simulation Section v4 — S3 (수정 버전): Group-wise Heteroskedasticity
#
# DGP: sigma2 = 1 + 2*G_s3  (G ~ Bernoulli(0.5), uncentered)
#      G=0 → sigma2_raw=1, G=1 → sigma2_raw=3
#      정규화 후: G=0 그룹 ≈ 0.5, G=1 그룹 ≈ 1.5 (약 3배 contrast)
#
# 기존 S3 버그: 1 + 2*X[,7]^2 에서 centered Bernoulli 제곱이 거의 상수 → 사실상 homo
# 수정: uncentered binary indicator G_s3 사용
#
# N = 10000 / k ∈ {200, 500, 1000, 2000} / n_rep = 300 / n0 = 200

source("code3/dgp.R")
source("code3/methods.R")
source("code3/metrics.R")

# ── Sanity Diagnostic ────────────────────────────────────────────────────────
cat("========== [S3 v2] Sanity Diagnostic ==========\n")
set.seed(42)
dat_check <- generate_data(N = 10000, case = 3)

cat("\nsummary(dat$sigma2):\n")
print(summary(dat_check$sigma2))

cat("\ntapply(sigma2, G_s3, mean):\n")
print(tapply(dat_check$sigma2, dat_check$G_s3, mean))

cat(sprintf("\ncor(sigma2, ell) = %.5f\n", cor(dat_check$sigma2, dat_check$ell)))
cat(sprintf("table(G_s3): G=0: %d, G=1: %d\n",
            sum(dat_check$G_s3 == 0), sum(dat_check$G_s3 == 1)))
cat(sprintf("length(unique(sigma2)) = %d\n\n", length(unique(round(dat_check$sigma2, 10)))))

# 정상 기준:
# - G=0 그룹 mean sigma2 ≈ 0.5, G=1 그룹 ≈ 1.5
# - cor(sigma2, ell) ≈ 0 (leverage와 독립)
# - unique sigma2 값이 2개가 아닌 다양해야 함 (N이 커서 비율 약간씩 다름)

# ── 실험 설정 ─────────────────────────────────────────────────────────────────
k_vals <- c(200, 500, 1000, 2000)
n_rep  <- 300
N      <- 10000
n0     <- 200

results <- list()
idx     <- 1

for (k_val in k_vals) {
  for (b in 1:n_rep) {

    set.seed(b * 1000 + 30 + k_val)   # case=3 seed 유지

    dat      <- generate_data(N = N, case = 3)
    dat_test <- generate_data(N = N, case = 3)

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
      m$case   <- 3
      m$k      <- k_val
      m$rep    <- b

      results[[idx]] <- m
      idx <- idx + 1
    }
  }
  cat(sprintf("S3v2 | k = %4d done\n", k_val))
}

sim_s3v2 <- do.call(rbind, results)

saveRDS(sim_s3v2, "results/sim_s3v2_results.rds")
write.csv(sim_s3v2, "results/sim_s3v2_results.csv", row.names = FALSE)

cat("Done. Total rows:", nrow(sim_s3v2), "\n")
# 예상 행 수: 4 k × 300 reps × 8 methods = 9,600

# ── 터미널 1차 확인 ───────────────────────────────────────────────────────────

method_order <- c("FULL", "UNI-IPW", "LEV-IPW", "OPT-homo",
                  "OPT-hetero-oracle", "OPT-hetero-plugin",
                  "OPT-hetero-tempering-0.5",
                  "OPT-plugin-stab(λ=0.5,γ=0.5)")
sim_s3v2$method <- factor(sim_s3v2$method, levels = method_order)

# 1. k별 순위 (mean excess risk, FULL 제외)
cat("\n\n========== [S3v2] 1. k별 순위 (mean excess risk, FULL 제외) ==========\n")
agg <- aggregate(excess_risk ~ k + method,
                 data = sim_s3v2[sim_s3v2$method != "FULL", ], FUN = mean)
agg <- agg[order(agg$k, agg$excess_risk), ]
for (kv in k_vals) {
  sub <- agg[agg$k == kv, ]
  cat(sprintf("\n  k = %4d\n", kv))
  for (i in seq_len(nrow(sub))) {
    cat(sprintf("    %d. %-40s %.5f\n", i, as.character(sub$method[i]), sub$excess_risk[i]))
  }
}

# 2. Tail diagnostic — mean vs. median
cat("\n\n========== [S3v2] 2. Tail diagnostic — mean vs. median ==========\n")
cat(sprintf("  %-40s %6s %8s %8s %8s\n", "method", "k", "mean", "median", "mean/med"))
for (kv in k_vals) {
  for (m in levels(sim_s3v2$method)[-1]) {
    er <- sim_s3v2$excess_risk[sim_s3v2$method == m & sim_s3v2$k == kv]
    if (length(er) == 0) next
    cat(sprintf("  %-40s %6d %8.5f %8.5f %8.3f\n",
                m, kv, mean(er), median(er), mean(er) / median(er)))
  }
}

# 3. π concentration diagnostic
cat("\n\n========== [S3v2] 3. π concentration diagnostic ==========\n")
cat(sprintf("  %-40s %6s %8s %10s %10s %10s\n",
            "method", "k", "clip_rt", "ESS_pi", "top1%_mass", "max/min_pi"))

for (kv in k_vals) {
  set.seed(1 * 1000 + 30 + kv)
  dat_diag <- generate_data(N = N, case = 3)

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
    clip_rt  <- mean(pv >= 1 - 1e-10)
    cat(sprintf("  %-40s %6d %8.4f %10.1f %10.4f %10.1f\n",
                nm, kv, clip_rt, ess, top_mass, mm_ratio))
  }
}

# 4. 분포 진단 — mean / median / sd / q90
cat("\n\n========== [S3v2] 4. 분포 진단 (excess risk) ==========\n")
for (kv in k_vals) {
  cat(sprintf("\n  k = %4d\n", kv))
  sub <- sim_s3v2[sim_s3v2$method != "FULL" & sim_s3v2$k == kv, ]
  cat(sprintf("  %-40s %8s %8s %8s %8s\n", "method", "mean", "median", "sd", "q90"))
  for (m in levels(sim_s3v2$method)[-1]) {
    er <- sub$excess_risk[sub$method == m]
    if (length(er) == 0) next
    cat(sprintf("  %-40s %8.5f %8.5f %8.5f %8.5f\n",
                m, mean(er), median(er), sd(er), quantile(er, 0.90)))
  }
}

# 5. Group별 sampling rate (S3 핵심 diagnostic)
cat("\n\n========== [S3v2] 5. G별 평균 π (k=500, b=1 기준) ==========\n")
set.seed(1 * 1000 + 30 + 500)
dat_g <- generate_data(N = N, case = 3)

pi_g_list <- list(
  "UNI-IPW"                  = rep(500 / N, N),
  "LEV-IPW"                  = compute_pi(dat_g$ell, k = 500),
  "OPT-homo"                 = compute_pi(sqrt(dat_g$ell), k = 500),
  "OPT-hetero-oracle"        = compute_pi(sqrt(dat_g$sigma2 * dat_g$ell), k = 500),
  "OPT-hetero-tempering-0.5" = compute_pi(sqrt(dat_g$sigma2 * dat_g$ell)^0.5, k = 500)
)

cat(sprintf("  %-40s %10s %10s %10s\n", "method", "G=0 mean_pi", "G=1 mean_pi", "ratio"))
for (nm in names(pi_g_list)) {
  pv   <- pi_g_list[[nm]]
  m0   <- mean(pv[dat_g$G_s3 == 0])
  m1   <- mean(pv[dat_g$G_s3 == 1])
  cat(sprintf("  %-40s %10.5f %10.5f %10.3f\n", nm, m0, m1, m1 / m0))
}