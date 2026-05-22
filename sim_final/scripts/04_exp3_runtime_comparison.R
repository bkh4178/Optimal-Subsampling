# 04_exp3_runtime_comparison.R
# Exp3: 단계별 런타임 분해 측정 + N별 Scaling 실험
#
# 설정 A (분해): N=100,000, n0=500, n_bench=50, k in {1000, 3000, 5000}
# 설정 B (scaling): N in {10000,50000,100000,500000}, k=1000, n_bench=30
# 단계: pilot_draw, pilot_fit, variance_est, moment_matrix,
#       score_compute, subsample_draw, final_fit
#
# 실행 위치: sim_final/ 폴더에서  source("scripts/04_exp3_runtime_comparison.R")

source("config/config_final.R")
source("R/00_utils.R")
source("R/01_dgp.R")
source("R/02_sampling_methods.R")
source("R/03_variance_estimation.R")
source("R/04_metrics.R")
source("R/06_runtime_utils.R")

ensure_dir(DIR_EXP3)

# ── 공통 파라미터 ─────────────────────────────────────────────────────────────
METHOD_ORDER_EXP3 <- c(
  "FULL", "SRS", "LEV-IPW", "OPT-homo",
  "OPT-hetero-oracle", "OPT-hetero-plugin"
)
STAGE_ORDER <- c(
  "pilot_draw", "pilot_fit", "variance_est",
  "moment_matrix", "score_compute", "subsample_draw", "final_fit"
)
PILOT_STAGES <- c("pilot_draw", "pilot_fit", "variance_est")

# ═════════════════════════════════════════════════════════════════════════════
# Exp3-a: 단계별 분해 (N=100,000, k 변화)
# ═════════════════════════════════════════════════════════════════════════════
N_EXP3   <- 100000L
N0_EXP3  <- 500L
N_BENCH  <- 50L
KS_BENCH <- c(1000L, 3000L, 5000L)

cat("====================================================\n")
cat(" Exp3-a: 단계별 런타임 분해\n")
cat("====================================================\n")
cat(sprintf(" N=%d, n0=%d, n_bench=%d\n", N_EXP3, N0_EXP3, N_BENCH))
cat(sprintf(" k: %s\n", paste(KS_BENCH, collapse = ", ")))
cat(sprintf(" %s\n", dgp_sigma_label))
cat("====================================================\n\n")

cat(sprintf("데이터 생성 중 (N = %d) ...\n", N_EXP3))
dat <- generate_data_final(N_EXP3, beta_true,
                           heteroscedastic = TRUE, seed = SEED_DATA)
cat("  완료.\n\n")

all_summary <- list()
all_raw     <- list()

for (kv in KS_BENCH) {
  cat(sprintf("[k = %d] 워밍업 3회 후 %d회 벤치마크 ...\n", kv, N_BENCH))
  bench           <- benchmark_decomposed(dat, kv, n0 = N0_EXP3,
                                          n_bench = N_BENCH,
                                          seed_start = 3001L)
  bench$raw$k     <- kv
  all_summary[[length(all_summary) + 1L]] <- bench$summary
  all_raw[[length(all_raw) + 1L]]         <- bench$raw
  cat(sprintf("  k = %d 완료.\n\n", kv))
}

summary_df <- do.call(rbind, all_summary)
raw_df     <- do.call(rbind, all_raw)

rds_path <- file.path(DIR_EXP3, "exp3_runtime_decomposed.rds")
saveRDS(summary_df, rds_path)
cat(sprintf("저장 완료: %s\n\n", rds_path))

# ── 출력 1: 방법별 Total Runtime ─────────────────────────────────────────────
total_rows <- list()
for (m in METHOD_ORDER_EXP3) {
  for (kv in KS_BENCH) {
    sub <- raw_df[raw_df$method == m & raw_df$k == kv, ]
    if (nrow(sub) == 0L) next
    rep_totals <- tapply(sub$elapsed, sub$bench_i, sum)
    total_rows[[length(total_rows) + 1L]] <- data.frame(
      method     = m,
      k          = kv,
      mean_total = mean(rep_totals),
      sd_total   = sd(rep_totals),
      stringsAsFactors = FALSE
    )
  }
}
total_df <- do.call(rbind, total_rows)

cat(sprintf("%s\n", strrep("=", 72)))
cat(" [표 1] 방법별 Total Runtime (초)\n")
cat(sprintf("%s\n", strrep("=", 72)))
header <- sprintf("  %-22s", "method")
for (kv in KS_BENCH) header <- paste0(header, sprintf("  k=%-8d", kv))
cat(header, "\n")
cat(strrep("-", nchar(header)), "\n")
for (m in METHOD_ORDER_EXP3) {
  line <- sprintf("  %-22s", m)
  for (kv in KS_BENCH) {
    sub <- total_df[total_df$method == m & total_df$k == kv, ]
    if (nrow(sub) == 0L) {
      line <- paste0(line, sprintf("  %10s", "-"))
    } else {
      line <- paste0(line,
                     sprintf("  %6.4f+%5.4f",
                             sub$mean_total, sub$sd_total))
    }
  }
  cat(line, "\n")
}

# ── 출력 2: 단계별 분해 표 (k별) ─────────────────────────────────────────────
for (kv in KS_BENCH) {
  cat(sprintf("\n%s\n", strrep("=", 72)))
  cat(sprintf(" [표 2] 단계별 분해 (k = %d, 단위: 초)\n", kv))
  cat(sprintf("%s\n", strrep("=", 72)))
  cat(sprintf("  %-22s %-16s %9s %9s %9s\n",
              "method", "stage", "mean", "sd", "median"))
  cat(strrep("-", 68), "\n")
  for (m in METHOD_ORDER_EXP3) {
    sub <- summary_df[summary_df$method == m & summary_df$k == kv, ]
    if (nrow(sub) == 0L) next
    sub$stage <- factor(sub$stage, levels = STAGE_ORDER)
    sub <- sub[order(sub$stage), ]
    for (j in seq_len(nrow(sub))) {
      cat(sprintf("  %-22s %-16s %9.5f %9.5f %9.5f\n",
                  if (j == 1L) m else "",
                  as.character(sub$stage[j]),
                  sub$mean_time[j],
                  sub$sd_time[j],
                  sub$median_time[j]))
    }
    cat(strrep("-", 68), "\n")
  }
}

# ── 출력 3: OPT-hetero-plugin 단계 비중 ──────────────────────────────────────
cat(sprintf("\n%s\n", strrep("=", 72)))
cat(" [표 3] OPT-hetero-plugin 단계 비중 (mean_time 기준)\n")
cat(sprintf("%s\n", strrep("=", 72)))
cat(sprintf("  %-8s %-16s %10s %8s\n",
            "k", "stage", "mean(s)", "pct(%)"))
cat(strrep("-", 46), "\n")

for (kv in KS_BENCH) {
  sub <- summary_df[
    summary_df$method == "OPT-hetero-plugin" & summary_df$k == kv,
  ]
  if (nrow(sub) == 0L) next
  sub$stage  <- factor(sub$stage, levels = STAGE_ORDER)
  sub        <- sub[order(sub$stage), ]
  total_mean <- sum(sub$mean_time)
  pilot_var  <- sum(
    sub$mean_time[as.character(sub$stage) %in% PILOT_STAGES]
  )
  first_row <- TRUE
  for (j in seq_len(nrow(sub))) {
    pct <- sub$mean_time[j] / total_mean * 100
    cat(sprintf("  %-8s %-16s %10.5f %7.1f%%\n",
                if (first_row) as.character(kv) else "",
                as.character(sub$stage[j]),
                sub$mean_time[j], pct))
    first_row <- FALSE
  }
  cat(sprintf("  %-8s %-16s %10.5f %7.1f%%  <- pilot+variance\n",
              "", "[pilot+var]",
              pilot_var, pilot_var / total_mean * 100))
  cat(sprintf("  %-8s %-16s %10.5f  100.0%%\n",
              "", "[total]", total_mean))
  cat(strrep("-", 46), "\n")
}

# ═════════════════════════════════════════════════════════════════════════════
# Exp3-b: N별 Scaling 실험
# ═════════════════════════════════════════════════════════════════════════════
N_GRID_SCALE  <- c(10000L, 50000L, 100000L, 500000L)
K_SCALE       <- 1000L
N0_SCALE      <- 500L
N_BENCH_SCALE <- 30L

cat("\n====================================================\n")
cat(" Exp3-b: N별 Scaling 실험\n")
cat("====================================================\n")
cat(sprintf(" N_grid: %s\n", paste(N_GRID_SCALE, collapse = ", ")))
cat(sprintf(" k=%d, n0=%d, n_bench=%d\n",
            K_SCALE, N0_SCALE, N_BENCH_SCALE))
cat("====================================================\n\n")

scaling_summary <- list()
scaling_raw     <- list()

for (nv in N_GRID_SCALE) {
  cat(sprintf("[N = %d] 데이터 생성 중 ...\n", nv))
  dat_n <- generate_data_final(nv, beta_true,
                               heteroscedastic = TRUE, seed = SEED_DATA)
  cat(sprintf("[N = %d] 워밍업 3회 후 %d회 벤치마크 ...\n",
              nv, N_BENCH_SCALE))
  bench_n             <- benchmark_decomposed(dat_n, K_SCALE,
                                              n0 = N0_SCALE,
                                              n_bench = N_BENCH_SCALE,
                                              seed_start = 4001L)
  bench_n$raw$k       <- K_SCALE
  bench_n$raw$N       <- nv
  bench_n$summary$N   <- nv
  scaling_summary[[length(scaling_summary) + 1L]] <- bench_n$summary
  scaling_raw[[length(scaling_raw) + 1L]]         <- bench_n$raw
  cat(sprintf("  N = %d 완료.\n\n", nv))
}

scaling_summary_df <- do.call(rbind, scaling_summary)
scaling_raw_df     <- do.call(rbind, scaling_raw)

scaling_rds <- file.path(DIR_EXP3, "exp3_runtime_scaling.rds")
saveRDS(scaling_summary_df, scaling_rds)
cat(sprintf("저장 완료: %s\n\n", scaling_rds))

# ── 출력 4: N별 방법별 Total Runtime ─────────────────────────────────────────
scaling_total_rows <- list()
for (m in METHOD_ORDER_EXP3) {
  for (nv in N_GRID_SCALE) {
    sub <- scaling_raw_df[
      scaling_raw_df$method == m & scaling_raw_df$N == nv,
    ]
    if (nrow(sub) == 0L) next
    rep_totals <- tapply(sub$elapsed, sub$bench_i, sum)
    scaling_total_rows[[length(scaling_total_rows) + 1L]] <- data.frame(
      method     = m,
      N          = nv,
      mean_total = mean(rep_totals),
      sd_total   = sd(rep_totals),
      stringsAsFactors = FALSE
    )
  }
}
scaling_total_df <- do.call(rbind, scaling_total_rows)

cat(sprintf("%s\n", strrep("=", 72)))
cat(" [표 4] N별 방법별 Total Runtime (초)\n")
cat(sprintf("%s\n", strrep("=", 72)))
header4 <- sprintf("  %-22s", "method")
for (nv in N_GRID_SCALE) {
  header4 <- paste0(header4, sprintf("  N=%-10d", nv))
}
cat(header4, "\n")
cat(strrep("-", nchar(header4)), "\n")
for (m in METHOD_ORDER_EXP3) {
  line <- sprintf("  %-22s", m)
  for (nv in N_GRID_SCALE) {
    sub <- scaling_total_df[
      scaling_total_df$method == m & scaling_total_df$N == nv,
    ]
    if (nrow(sub) == 0L) {
      line <- paste0(line, sprintf("  %12s", "-"))
    } else {
      line <- paste0(line, sprintf("  %12.4f", sub$mean_total))
    }
  }
  cat(line, "\n")
}

# ── 출력 5: N별 OPT-hetero-plugin 단계 분해 ──────────────────────────────────
cat(sprintf("\n%s\n", strrep("=", 72)))
cat(" [표 5] OPT-hetero-plugin 단계 분해 (N별, k=1000)\n")
cat(sprintf("%s\n", strrep("=", 72)))
cat(sprintf("  %-12s %-16s %10s %8s\n",
            "N", "stage", "mean(s)", "pct(%)"))
cat(strrep("-", 50), "\n")
for (nv in N_GRID_SCALE) {
  sub <- scaling_summary_df[
    scaling_summary_df$method == "OPT-hetero-plugin" &
      scaling_summary_df$N == nv,
  ]
  if (nrow(sub) == 0L) next
  sub$stage  <- factor(sub$stage, levels = STAGE_ORDER)
  sub        <- sub[order(sub$stage), ]
  total_mean <- sum(sub$mean_time)
  first_row  <- TRUE
  for (j in seq_len(nrow(sub))) {
    pct <- sub$mean_time[j] / total_mean * 100
    cat(sprintf("  %-12s %-16s %10.5f %7.1f%%\n",
                if (first_row) as.character(nv) else "",
                as.character(sub$stage[j]),
                sub$mean_time[j], pct))
    first_row <- FALSE
  }
  cat(strrep("-", 50), "\n")
}

# ── 출력 6: variance_est 비중 N별 변화 ───────────────────────────────────────
cat(sprintf("\n%s\n", strrep("=", 60)))
cat(" [표 6] variance_est 비중의 N별 변화 (OPT-hetero-plugin)\n")
cat(sprintf("%s\n", strrep("=", 60)))
cat(sprintf("  %-12s %12s %12s %10s\n",
            "N", "var_est(s)", "total(s)", "pct(%)"))
cat(strrep("-", 50), "\n")
for (nv in N_GRID_SCALE) {
  sub <- scaling_summary_df[
    scaling_summary_df$method == "OPT-hetero-plugin" &
      scaling_summary_df$N == nv,
  ]
  if (nrow(sub) == 0L) next
  total_mean <- sum(sub$mean_time)
  var_time   <- sub$mean_time[sub$stage == "variance_est"]
  if (length(var_time) == 0L) var_time <- 0
  cat(sprintf("  %-12d %12.5f %12.5f %9.1f%%\n",
              nv, var_time, total_mean,
              var_time / total_mean * 100))
}

cat("\nExp3 완료.\n")
