# run_sim_p3.R
# Phase 1 시뮬레이션 실행 — p = 3
#
# 구조:
#   case in {1, 2, 3}  ×  k in {200, 500, 1000, 2000}  ×  500 replicates
#
# 공변량: X1 ~ N(0,1), X2 ~ Bernoulli(0.5) centered, X3 ~ N(0,1)
# beta0 = (0.5, 1.0, 2.0)
#
# 저장: sim/results_p3.rds

source("sim/dgp.R")
source("sim/methods.R")
source("sim/metrics.R")

# ── 설정 ──────────────────────────────────────────────────────────────────────

N       <- 10000
N_test  <- 5000
p       <- 3
beta0   <- c(0.5, 1.0, 2.0)
n0      <- 200
n_rep   <- 500
cases   <- c(1, 2, 3)
ks      <- c(200, 500, 1000, 2000)
gamma   <- 0.5   # Case 2: sigma2 = 1 + gamma * |X1|

# ── 실행 ──────────────────────────────────────────────────────────────────────

results    <- list()
idx_result <- 1L
time_start <- proc.time()

cat("====================================================\n")
cat(sprintf(" Phase 1 시뮬레이션 시작 (p = %d)\n", p))
cat(sprintf(" N = %d / N_test = %d / n0 = %d / n_rep = %d\n", N, N_test, n0, n_rep))
cat(sprintf(" beta0: (%s)\n", paste(beta0, collapse = ", ")))
cat(sprintf(" cases: %s\n", paste(cases, collapse = ", ")))
cat(sprintf(" k:     %s\n", paste(ks, collapse = ", ")))
cat("====================================================\n")

for (case in cases) {

  cat(sprintf("\n[Case %d 시작]\n", case))
  time_case <- proc.time()

  for (k in ks) {

    time_k <- proc.time()

    for (rep in seq_len(n_rep)) {

      # train data
      dat <- generate_data(N = N, p = p, case = case, gamma = gamma, beta0 = beta0)

      # test data (같은 DGP, 독립 생성)
      dat_test <- generate_data(N = N_test, p = p, case = case, gamma = gamma, beta0 = beta0)
      X_test   <- dat_test$X
      y_test   <- dat_test$y

      # 6개 method 실행 및 지표 계산
      methods_out <- list(
        run_full(dat),
        run_uni_ipw(dat, k),
        run_lev_ipw(dat, k),
        run_opt_homo(dat, k),
        run_opt_hetero_oracle(dat, k),
        run_opt_hetero_plugin(dat, k, n0 = n0)
      )

      for (res in methods_out) {
        row      <- compute_metrics(res, dat, X_test, y_test)
        row$case <- case
        row$k    <- k
        row$rep  <- rep
        results[[idx_result]] <- row
        idx_result <- idx_result + 1L
      }
    }

    elapsed_k <- round((proc.time() - time_k)["elapsed"], 1)
    cat(sprintf("  k = %4d 완료  (%5.1f초)\n", k, elapsed_k))
  }

  elapsed_case <- round((proc.time() - time_case)["elapsed"], 1)
  cat(sprintf("[Case %d 완료 — %.1f초]\n", case, elapsed_case))
}

# ── 저장 ──────────────────────────────────────────────────────────────────────

results_df <- do.call(rbind, results)
saveRDS(results_df, "sim/results_p3.rds")

elapsed_total <- round((proc.time() - time_start)["elapsed"], 1)

cat("\n====================================================\n")
cat(sprintf(" 저장 완료: sim/results_p3.rds\n"))
cat(sprintf(" 총 행 수:  %d\n", nrow(results_df)))
cat(sprintf(" 총 소요:   %.1f초\n", elapsed_total))
cat("====================================================\n")

# ── 요약 테이블 (case × k × method 별 평균 excess_risk + 순위) ───────────────

cat("\n[요약] 평균 excess_risk — case × k × method  (괄호: 순위, FULL 제외)\n")

method_order <- c("FULL", "UNI-IPW", "LEV-IPW", "OPT-homo",
                  "OPT-hetero-oracle", "OPT-hetero-plugin")
rank_methods <- setdiff(method_order, "FULL")

agg <- aggregate(excess_risk ~ case + k + method, data = results_df, FUN = mean)

for (cs in cases) {
  cat(sprintf("\n======= Case %d =======\n", cs))

  cat(sprintf("%-22s", "method \\ k"))
  for (k in ks) cat(sprintf("%16d", k))
  cat("\n")
  cat(strrep("-", 22 + 16 * length(ks)), "\n")

  # k별 순위 계산 (FULL 제외)
  rank_table <- list()
  for (k in ks) {
    sub <- agg[agg$case == cs & agg$k == k & agg$method %in% rank_methods, ]
    sub <- sub[order(sub$excess_risk), ]
    rank_table[[as.character(k)]] <- setNames(seq_len(nrow(sub)), sub$method)
  }

  for (m in method_order) {
    cat(sprintf("%-22s", m))
    for (k in ks) {
      val <- agg$excess_risk[agg$case == cs & agg$k == k & agg$method == m]
      if (length(val) == 0) {
        cat(sprintf("%16s", "NA"))
      } else if (m == "FULL") {
        cat(sprintf("%16.5f", val))
      } else {
        rk <- rank_table[[as.character(k)]][m]
        cat(sprintf("  %8.5f(%d등)", val, rk))
      }
    }
    cat("\n")
  }
}

# ── SE 확인 (Monte Carlo noise 점검) ─────────────────────────────────────────

cat("\n[SE 확인] case × k 별 excess_risk 표준오차\n")

agg_se <- aggregate(excess_risk ~ case + k + method, data = results_df,
                    FUN = function(x) sd(x) / sqrt(length(x)))
names(agg_se)[names(agg_se) == "excess_risk"] <- "se"

for (cs in cases) {
  cat(sprintf("\n--- Case %d ---\n", cs))
  cat(sprintf("%-22s", "method \\ k"))
  for (k in ks) cat(sprintf("%12d", k))
  cat("\n")

  for (m in method_order) {
    cat(sprintf("%-22s", m))
    for (k in ks) {
      val <- agg_se$se[agg_se$case == cs & agg_se$k == k & agg_se$method == m]
      if (length(val) == 0) cat(sprintf("%12s", "NA")) else cat(sprintf("%12.6f", val))
    }
    cat("\n")
  }
}

cat("\n")