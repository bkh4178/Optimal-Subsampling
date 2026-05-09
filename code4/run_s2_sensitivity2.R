# run_s2_sensitivity.R
# S2 Gamma Sensitivity — sigma2 = 1 + gamma * (ell / mean(ell))
#
# gamma ∈ {0.3, 0.7, 2.0}
#   gamma=0.3: mild contrast
#   gamma=0.7: moderate-strong contrast
#   gamma=2.0: strong contrast (concentration stress-test)
#
# N = 10000 / k ∈ {200, 500, 1000, 2000} / n_rep = 300 / n0 = 200

source("code4/dgp.R")
source("code4/methods.R")
source("code4/metrics_a.R")

k_vals     <- c(200, 500, 1000, 2000)
n_rep      <- 300
N          <- 10000
n0         <- 200
gamma_vals <- c(0.3, 0.7, 2.0)

method_order <- c("FULL", "UNI-IPW", "LEV-IPW", "OPT-homo",
                  "OPT-hetero-oracle", "OPT-hetero-plugin",
                  "OPT-hetero-tempering-0.5",
                  "OPT-plugin-stab(λ=0.5,γ=0.5)")

# ── Sanity Diagnostic (gamma별 sigma2 분포 및 π concentration) ───────────────
cat("========== Sanity Diagnostic: gamma별 sigma2 및 π concentration ==========\n")
cat(sprintf("\n%-8s %-22s %8s %8s %8s %8s %10s %10s %10s\n",
            "gamma", "summary(sigma2)", "min", "median", "max",
            "cor(s2,ell)", "oracle_ESS", "top1%_mass", "max/min_pi"))

set.seed(42)
for (gv in gamma_vals) {
  dat_d <- generate_data(N = N, case = 2, gamma = gv)
  pi_o  <- compute_pi(sqrt(dat_d$sigma2 * dat_d$ell), k = 500)
  ess_o <- sum(pi_o)^2 / sum(pi_o^2)
  top_i <- pi_o > quantile(pi_o, 0.99)
  top_m <- sum(pi_o[top_i]) / sum(pi_o)
  mm    <- max(pi_o) / (min(pi_o[pi_o > 0]) + 1e-15)
  cat(sprintf("%-8.1f [%.3f, %.3f, %.3f] %8.4f %10.1f %10.4f %10.1f\n",
              gv,
              min(dat_d$sigma2), median(dat_d$sigma2), max(dat_d$sigma2),
              cor(dat_d$sigma2, dat_d$ell),
              ess_o, top_m, mm))
}

# ── 실험 실행 ─────────────────────────────────────────────────────────────────
all_results <- list()

for (gv in gamma_vals) {
  results <- list()
  idx     <- 1

  for (k_val in k_vals) {
    for (b in 1:n_rep) {

      set.seed(b * 1000 + 20 + k_val + round(gv * 10))

      dat      <- generate_data(N = N, case = 2, gamma = gv)
      dat_test <- generate_data(N = N, case = 2, gamma = gv)

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
        m       <- compute_metrics(res, dat, dat_test$X, dat_test$y)
        m$case  <- "S2"
        m$gamma <- gv
        m$k     <- k_val
        m$rep   <- b
        results[[idx]] <- m
        idx <- idx + 1
      }
    }
    cat(sprintf("S2 gamma=%.1f | k = %4d done\n", gv, k_val))
  }

  sim_g <- do.call(rbind, results)
  saveRDS(sim_g, sprintf("results/sim_s2_gamma%.0f_results.rds", gv * 10))
  write.csv(sim_g, sprintf("results/sim_s2_gamma%.0f_results.csv", gv * 10), row.names = FALSE)
  all_results[[as.character(gv)]] <- sim_g
}

cat("\nAll simulations done.\n")

# ── 터미널 비교 출력 ──────────────────────────────────────────────────────────

sim_all <- do.call(rbind, all_results)
sim_all$method <- factor(sim_all$method, levels = method_order)

# 1. gamma × k별 순위표 (FULL 제외)
cat("\n\n========== 1. gamma × k별 순위 (mean excess risk) ==========\n")
for (gv in gamma_vals) {
  cat(sprintf("\n  ── gamma = %.1f ──\n", gv))
  agg <- aggregate(excess_risk ~ k + method,
                   data = sim_all[sim_all$gamma == gv & sim_all$method != "FULL", ],
                   FUN = mean)
  agg <- agg[order(agg$k, agg$excess_risk), ]
  for (kv in k_vals) {
    sub <- agg[agg$k == kv, ]
    cat(sprintf("\n    k = %4d\n", kv))
    for (i in seq_len(nrow(sub))) {
      cat(sprintf("      %d. %-40s %.5f\n", i,
                  as.character(sub$method[i]), sub$excess_risk[i]))
    }
  }
}

# 2. 핵심 method × gamma × k 비교표 (method간 cross-gamma 비교용)
cat("\n\n========== 2. 핵심 method 비교 (mean excess risk) ==========\n")
key_methods <- c("UNI-IPW", "OPT-homo", "OPT-hetero-oracle",
                 "OPT-hetero-tempering-0.5", "OPT-plugin-stab(λ=0.5,γ=0.5)")

cat(sprintf("  %-40s %6s %8s %8s %8s\n", "method", "k",
            "g=0.3", "g=0.7", "g=2.0"))
cat(sprintf("  %s\n", paste(rep("-", 80), collapse = "")))

for (kv in k_vals) {
  for (m in key_methods) {
    vals <- sapply(gamma_vals, function(gv) {
      er <- sim_all$excess_risk[sim_all$method == m &
                                sim_all$gamma == gv &
                                sim_all$k == kv]
      if (length(er) == 0) NA else mean(er)
    })
    cat(sprintf("  %-40s %6d %8.5f %8.5f %8.5f\n", m, kv, vals[1], vals[2], vals[3]))
  }
  cat("\n")
}

# 3. Tail diagnostic (mean/median ratio) — gamma별
cat("\n\n========== 3. Tail diagnostic (mean/median ratio) ==========\n")
cat(sprintf("  %-40s %6s %6s %8s %8s %8s\n",
            "method", "k", "gamma", "mean", "median", "mean/med"))
for (gv in gamma_vals) {
  for (kv in c(500, 1000)) {
    for (m in key_methods) {
      er <- sim_all$excess_risk[sim_all$method == m &
                                sim_all$gamma == gv &
                                sim_all$k == kv]
      if (length(er) == 0) next
      cat(sprintf("  %-40s %6d %6.1f %8.5f %8.5f %8.3f\n",
                  m, kv, gv, mean(er), median(er), mean(er) / median(er)))
    }
    cat("\n")
  }
}

# 4. π concentration diagnostic — gamma별 (k=500 기준)
cat("\n\n========== 4. π concentration diagnostic (k=500, single replicate) ==========\n")
cat(sprintf("  %-40s %6s %10s %10s %10s\n",
            "method", "gamma", "ESS_pi", "top1%_mass", "max/min_pi"))

set.seed(42)
for (gv in gamma_vals) {
  dat_diag <- generate_data(N = N, case = 2, gamma = gv)
  pi_list <- list(
    "UNI-IPW"                  = rep(500 / N, N),
    "OPT-homo"                 = compute_pi(sqrt(dat_diag$ell), k = 500),
    "OPT-hetero-oracle"        = compute_pi(sqrt(dat_diag$sigma2 * dat_diag$ell), k = 500),
    "OPT-hetero-tempering-0.5" = compute_pi(sqrt(dat_diag$sigma2 * dat_diag$ell)^0.5, k = 500)
  )
  for (nm in names(pi_list)) {
    pv      <- pi_list[[nm]]
    ess     <- sum(pv)^2 / sum(pv^2)
    top_idx <- pv > quantile(pv, 0.99)
    top_mass <- if (sum(top_idx) == 0) 0.01 else sum(pv[top_idx]) / sum(pv)
    mm_ratio <- max(pv) / (min(pv[pv > 0]) + 1e-15)
    cat(sprintf("  %-40s %6.1f %10.1f %10.4f %10.1f\n",
                nm, gv, ess, top_mass, mm_ratio))
  }
  cat("\n")
}

# 5. Oracle vs. OPT-homo gap (gamma별 oracle이 homo 대비 얼마나 나쁜지)
cat("\n\n========== 5. oracle vs. OPT-homo excess risk gap ==========\n")
cat(sprintf("  %6s %6s %10s %10s %10s\n",
            "gamma", "k", "oracle", "homo", "oracle/homo"))
for (gv in gamma_vals) {
  for (kv in k_vals) {
    er_o <- mean(sim_all$excess_risk[sim_all$method == "OPT-hetero-oracle" &
                                     sim_all$gamma == gv & sim_all$k == kv])
    er_h <- mean(sim_all$excess_risk[sim_all$method == "OPT-homo" &
                                     sim_all$gamma == gv & sim_all$k == kv])
    cat(sprintf("  %6.1f %6d %10.5f %10.5f %10.3f\n",
                gv, kv, er_o, er_h, er_o / er_h))
  }
  cat("\n")
}

# 6. C1(pi) vs. excess risk 비교
# 핵심 확인: oracle이 C1(pi)에서는 가장 좋은데 excess risk에서 나쁘면
#            → leading term은 최소화했지만 finite-sample remainder가 커진 것
#            oracle이 C1(pi)에서도 나쁘면 → 더 심각한 문제 (구현 또는 이론 불일치)
cat("\n\n========== 6. C1(pi) vs. excess risk 비교 (FULL 제외) ==========\n")
cat("  C1(pi) = sum(sigma2_i * ell_i / pi_i)  — theoretical leading term\n")
cat("  이론: oracle이 C1에서 가장 작아야 함\n\n")

cat(sprintf("  %-40s %6s %6s %12s %12s %10s\n",
            "method", "k", "gamma", "mean_C1", "mean_ER", "ER/C1_ratio"))

for (gv in gamma_vals) {
  cat(sprintf("\n  --- gamma = %.1f ---\n", gv))
  for (kv in c(500, 1000)) {
    cat(sprintf("  k = %d\n", kv))
    # C1 순위 (오름차순)
    c1_vals <- sapply(key_methods, function(m) {
      v <- sim_all$c1_pi[sim_all$method == m &
                         sim_all$gamma == gv & sim_all$k == kv]
      if (length(v) == 0) NA else mean(v)
    })
    er_vals <- sapply(key_methods, function(m) {
      v <- sim_all$excess_risk[sim_all$method == m &
                               sim_all$gamma == gv & sim_all$k == kv]
      if (length(v) == 0) NA else mean(v)
    })
    ord <- order(c1_vals, na.last = TRUE)
    for (i in ord) {
      m <- key_methods[i]
      cat(sprintf("    %-40s %12.1f %12.5f\n", m, c1_vals[i], er_vals[i]))
    }
    cat("\n")
  }
}

# C1 순위와 excess risk 순위 비교 요약
cat("\n========== 6b. C1 rank vs. ER rank 비교 (k=500) ==========\n")
cat(sprintf("  %-40s %6s %8s %8s\n", "method", "gamma", "C1_rank", "ER_rank"))
for (gv in gamma_vals) {
  c1_vals <- sapply(key_methods, function(m) {
    v <- sim_all$c1_pi[sim_all$method == m &
                       sim_all$gamma == gv & sim_all$k == 500]
    if (length(v) == 0) NA else mean(v)
  })
  er_vals <- sapply(key_methods, function(m) {
    v <- sim_all$excess_risk[sim_all$method == m &
                             sim_all$gamma == gv & sim_all$k == 500]
    if (length(v) == 0) NA else mean(v)
  })
  c1_rank <- rank(c1_vals, na.last = TRUE)
  er_rank <- rank(er_vals, na.last = TRUE)
  for (i in seq_along(key_methods)) {
    cat(sprintf("  %-40s %6.1f %8.0f %8.0f\n",
                key_methods[i], gv, c1_rank[i], er_rank[i]))
  }
  cat("\n")
}