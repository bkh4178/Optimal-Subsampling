# 07_exp5_misspec_n0_sensitivity.R
# Exp5: Misspecification × n0 Sensitivity
#
# 목적: 분산 모형 오명세 상황에서 n0가 커질수록
#       score alignment가 회복되는지, oracle gap이 줄어드는지 확인.
#   mild_omitted  : n0를 늘리면 어느 정도 커버됨
#   severe_omitted: n0를 늘려도 구조적으로 회복되지 않음
#
# 고정: heteroscedastic DGP, k=1000, n_rep=500
# 변화: n0 in c(100, 200, 400, 800)
# Plugin 3종: correct_scale, mild_omitted, severe_omitted
# 베이스라인 5종: FULL, SRS, LEV-IPW, OPT-homo, OPT-hetero-oracle
#
# 실행 위치: sim_final/ 폴더에서  source("scripts/07_exp5_misspec_n0_sensitivity.R")

source("config/config_final.R")
source("R/00_utils.R")
source("R/01_dgp.R")
source("R/02_sampling_methods.R")
source("R/03_variance_estimation.R")
source("R/04_metrics.R")
source("R/05_mcse_summary.R")

DIR_EXP5 <- file.path("results", "raw", "exp5")
ensure_dir(DIR_EXP5)
ensure_dir(file.path("results", "summary"))

# ── 실험 파라미터 ─────────────────────────────────────────────────────────────

K_EXP5       <- 1000L
N0_GRID_EXP5 <- c(50L, 100L, 200L, 300L, 400L, 500L)

PLUGIN_SETTINGS_EXP5 <- c("correct_scale", "mild_omitted", "severe_omitted")
PLUGIN_METHODS_EXP5  <- c("plugin-correct", "plugin-mild", "plugin-severe")
BASE_METHODS_EXP5 <- c(
  "FULL", "SRS", "LEV-IPW", "OPT-homo", "OPT-hetero-oracle"
)
ALL_METHODS_EXP5 <- c(BASE_METHODS_EXP5, PLUGIN_METHODS_EXP5)

# ── abs_s 보조 함수 ──────────────────────────────────────────────────────────
# X 열 순서: (intercept, Z1, Z2, Z3, Z4, Z5, Z6, Z7, Z8, Z9)
#   col:          1       2   3   4   5   6   7   8   9   10

make_abs_s_mild <- function(X) {
  # Z8 누락: |0.5*Z1 + 0.5*Z2 + Z6|
  abs(0.5 * X[, 2L] + 0.5 * X[, 3L] + X[, 7L])
}

make_abs_s_severe <- function(X) {
  # Z6, Z8 누락: |0.5*Z1 + 0.5*Z2|
  abs(0.5 * X[, 2L] + 0.5 * X[, 3L])
}

# ── Plugin 추정 함수 ─────────────────────────────────────────────────────────

.clip_positive <- function(v) {
  bad <- !is.finite(v) | v <= 0
  if (any(bad)) {
    good_vals <- v[!bad]
    v[bad]    <- if (length(good_vals) > 0L) median(good_vals) else 1
  }
  pmax(v, 1e-10)
}

# make_abs_s (정명세용) 는 R/03_variance_estimation.R 에서 정의됨
est_sigma2_correct_scale <- function(dat, idx_pilot, beta_pilot) {
  resid_pilot <- dat$y[idx_pilot] -
    as.numeric(dat$X[idx_pilot, , drop = FALSE] %*% beta_pilot)
  abs_s_pilot <- make_abs_s(dat$X[idx_pilot, , drop = FALSE])
  pilot_df    <- data.frame(abs_resid = abs(resid_pilot), abs_s = abs_s_pilot)
  sig_model   <- lm(abs_resid ~ abs_s, data = pilot_df)
  sigma_hat   <- as.numeric(predict(sig_model,
    newdata = data.frame(abs_s = make_abs_s(dat$X))))
  .clip_positive(sigma_hat)^2
}

est_sigma2_mild_omitted <- function(dat, idx_pilot, beta_pilot) {
  resid_pilot <- dat$y[idx_pilot] -
                 as.numeric(dat$X[idx_pilot, , drop = FALSE] %*% beta_pilot)
  abs_s_pilot <- make_abs_s_mild(dat$X[idx_pilot, , drop = FALSE])
  pilot_df    <- data.frame(abs_resid = abs(resid_pilot), abs_s = abs_s_pilot)
  sig_model   <- lm(abs_resid ~ abs_s, data = pilot_df)
  sigma_hat   <- as.numeric(predict(sig_model,
                   newdata = data.frame(abs_s = make_abs_s_mild(dat$X))))
  .clip_positive(sigma_hat)^2
}

est_sigma2_severe_omitted <- function(dat, idx_pilot, beta_pilot) {
  resid_pilot <- dat$y[idx_pilot] -
                 as.numeric(dat$X[idx_pilot, , drop = FALSE] %*% beta_pilot)
  abs_s_pilot <- make_abs_s_severe(dat$X[idx_pilot, , drop = FALSE])
  pilot_df    <- data.frame(abs_resid = abs(resid_pilot), abs_s = abs_s_pilot)
  sig_model   <- lm(abs_resid ~ abs_s, data = pilot_df)
  sigma_hat   <- as.numeric(predict(sig_model,
                   newdata = data.frame(abs_s = make_abs_s_severe(dat$X))))
  .clip_positive(sigma_hat)^2
}

# ── 단일 replicate: 8 methods (5 baseline + 3 plugin) ────────────────────────

run_one_rep_exp5 <- function(dat, k, X_test, mu_test, y_test, n0) {
  N    <- nrow(dat$X)
  beta <- dat$beta

  idx_pilot  <- sample.int(N, n0, replace = FALSE)
  beta_pilot <- fit_ipw(dat$X, dat$y, idx_pilot, pi = rep(n0 / N, N))
  if (is.null(beta_pilot)) beta_pilot <- beta

  sigma2_correct <- est_sigma2_correct_scale(dat, idx_pilot, beta_pilot)
  sigma2_mild    <- est_sigma2_mild_omitted(dat, idx_pilot, beta_pilot)
  sigma2_severe  <- est_sigma2_severe_omitted(dat, idx_pilot, beta_pilot)

  oracle_score <- sqrt(dat$sigma2_vec * dat$ell)

  make_diag <- function(sigma2_hat) {
    list(
      sigma_cor = safe_cor(sigma2_hat, dat$sigma2_vec),
      score_cor = safe_cor(sqrt(sigma2_hat * dat$ell), oracle_score)
    )
  }
  diag_correct <- make_diag(sigma2_correct)
  diag_mild    <- make_diag(sigma2_mild)
  diag_severe  <- make_diag(sigma2_severe)

  run_method <- function(method_name, misspec_setting, score,
                         sigma_cor = NA_real_, score_cor = NA_real_) {
    samp       <- sample_bernoulli(score, k)
    idx_final  <- samp$idx
    pi_raw     <- samp$pi
    n_clip     <- sum(pi_raw >= 1 - 1e-12)

    nz         <- pi_raw[pi_raw > 0]
    ess_pi     <- if (length(nz) > 0L) sum(pi_raw)^2 / sum(pi_raw^2) else NA_real_
    max_weight <- if (length(nz) > 0L) max(1 / pmax(nz, 1e-10))     else NA_real_

    beta_hat <- fit_ipw(dat$X, dat$y, idx_final, pi = pi_raw)

    if (is.null(beta_hat)) {
      return(data.frame(
        method          = method_name,
        misspec_setting = misspec_setting,
        fail            = TRUE,
        excess_risk     = NA_real_,
        param_mse       = NA_real_,
        sigma_cor       = sigma_cor,
        score_cor       = score_cor,
        n_realized      = length(idx_final),
        n_clip          = n_clip,
        ESS_pi          = ess_pi,
        max_weight      = max_weight,
        stringsAsFactors = FALSE
      ))
    }

    db <- beta_hat - beta
    data.frame(
      method          = method_name,
      misspec_setting = misspec_setting,
      fail            = FALSE,
      excess_risk     = as.numeric(t(db) %*% dat$A_hat %*% db),
      param_mse       = mean(db^2),
      sigma_cor       = sigma_cor,
      score_cor       = score_cor,
      n_realized      = length(idx_final),
      n_clip          = n_clip,
      ESS_pi          = ess_pi,
      max_weight      = max_weight,
      stringsAsFactors = FALSE
    )
  }

  # FULL (전체 OLS)
  full_samp <- sample_full(dat)
  bhat_full <- full_samp$beta_hat
  full_result <- if (is.null(bhat_full)) {
    data.frame(
      method = "FULL", misspec_setting = "baseline", fail = TRUE,
      excess_risk = NA_real_, param_mse = NA_real_,
      sigma_cor = NA_real_, score_cor = NA_real_,
      n_realized = full_samp$n_realized, n_clip = 0L,
      ESS_pi = NA_real_, max_weight = NA_real_,
      stringsAsFactors = FALSE
    )
  } else {
    db <- bhat_full - beta
    data.frame(
      method          = "FULL",
      misspec_setting = "baseline",
      fail            = FALSE,
      excess_risk     = as.numeric(t(db) %*% dat$A_hat %*% db),
      param_mse       = mean(db^2),
      sigma_cor       = NA_real_,
      score_cor       = NA_real_,
      n_realized      = full_samp$n_realized,
      n_clip          = 0L,
      ESS_pi          = NA_real_,
      max_weight      = NA_real_,
      stringsAsFactors = FALSE
    )
  }

  rbind(
    full_result,
    run_method("SRS",               "baseline",
               rep(1, N)),
    run_method("LEV-IPW",           "baseline",
               dat$ell),
    run_method("OPT-homo",          "baseline",
               sqrt(dat$ell)),
    run_method("OPT-hetero-oracle", "baseline",
               sqrt(dat$sigma2_vec * dat$ell)),
    run_method("plugin-correct",    "correct_scale",
               sqrt(sigma2_correct * dat$ell),
               sigma_cor = diag_correct$sigma_cor,
               score_cor = diag_correct$score_cor),
    run_method("plugin-mild",       "mild_omitted",
               sqrt(sigma2_mild * dat$ell),
               sigma_cor = diag_mild$sigma_cor,
               score_cor = diag_mild$score_cor),
    run_method("plugin-severe",     "severe_omitted",
               sqrt(sigma2_severe * dat$ell),
               sigma_cor = diag_severe$sigma_cor,
               score_cor = diag_severe$score_cor)
  )
}

# ── 데이터 생성 (고정) ────────────────────────────────────────────────────────

cat("====================================================\n")
cat(" Exp5: Misspecification × n0 Sensitivity\n")
cat("====================================================\n")
cat(sprintf(" N=%d, k=%d, n_rep=%d\n", N, K_EXP5, n_rep))
cat(sprintf(" n0_grid: %s\n", paste(N0_GRID_EXP5, collapse = ", ")))
cat(sprintf(" %s\n", dgp_sigma_label))
cat(" Plugin 3종: correct_scale / mild_omitted / severe_omitted\n")
cat("====================================================\n\n")

dat      <- generate_data_final(N,      beta_true, heteroscedastic = TRUE, seed = SEED_DATA)
dat_test <- generate_data_final(N_test, beta_true, heteroscedastic = TRUE, seed = SEED_TEST)
X_test   <- dat_test$X
mu_test  <- dat_test$mu
y_test   <- dat_test$y

# ── 메인 루프: n0 × rep ───────────────────────────────────────────────────────

all_rows <- list()
t_start  <- proc.time()

for (n0v in N0_GRID_EXP5) {
  cat(sprintf("── n0 = %d ─────────────────────────────\n", n0v))
  for (r in seq_len(n_rep)) {
    set.seed(SEED_REP + r)
    df_r      <- run_one_rep_exp5(dat, K_EXP5, X_test, mu_test, y_test, n0 = n0v)
    df_r$rep  <- r
    df_r$k    <- K_EXP5
    df_r$n0   <- n0v
    all_rows[[length(all_rows) + 1L]] <- df_r
    if (r %% 100L == 0L) cat(sprintf("  %d / %d\n", r, n_rep))
  }
}

results_df <- do.call(rbind, all_rows)

# 열 순서 정리 (요구 스펙 순)
results_df <- results_df[, c("method", "misspec_setting", "n0", "rep", "k",
                              "fail", "excess_risk", "param_mse",
                              "sigma_cor", "score_cor",
                              "n_realized", "n_clip", "ESS_pi", "max_weight")]

rds_path <- file.path(DIR_EXP5, "exp5_misspec_n0.rds")
saveRDS(results_df, rds_path)
cat(sprintf("\n저장 완료: %s (%d행, %.1f초)\n",
            rds_path, nrow(results_df), (proc.time() - t_start)[["elapsed"]]))

# ── 요약 출력 ─────────────────────────────────────────────────────────────────

valid <- results_df[!results_df$fail, ]

for (n0v in N0_GRID_EXP5) {
  cat(sprintf("\n%s\n", strrep("=", 70)))
  cat(sprintf(" n0 = %d\n", n0v))
  cat(sprintf("%s\n", strrep("=", 70)))

  vn <- valid[valid$n0 == n0v, ]

  # excess_risk 요약
  cat(sprintf("  %-25s %12s %12s\n", "method", "mean_excess", "mcse"))
  cat(strrep("-", 52), "\n")
  for (m in ALL_METHODS_EXP5) {
    sub <- vn[vn$method == m, ]
    if (nrow(sub) == 0L) next
    cat(sprintf("  %-25s %12.6f %12.6f\n",
                m,
                mean(sub$excess_risk, na.rm = TRUE),
                compute_mcse(sub$excess_risk)))
  }

  # sigma_cor / score_cor (plugin만)
  cat(sprintf("\n  %-20s %12s %12s\n", "misspec_setting", "sigma_cor", "score_cor"))
  cat(strrep("-", 46), "\n")
  for (i in seq_along(PLUGIN_SETTINGS_EXP5)) {
    sub <- vn[vn$method == PLUGIN_METHODS_EXP5[i], ]
    if (nrow(sub) == 0L) next
    cat(sprintf("  %-20s %12.4f %12.4f\n",
                PLUGIN_SETTINGS_EXP5[i],
                mean(sub$sigma_cor, na.rm = TRUE),
                mean(sub$score_cor, na.rm = TRUE)))
  }

  # Paired diff: plugin - oracle, homo - plugin
  cat(sprintf("\n  %-16s %-22s %10s %10s %10s %10s\n",
              "setting", "contrast",
              "mean_diff", "mcse_diff", "lower", "upper"))
  cat(strrep("-", 72), "\n")
  for (i in seq_along(PLUGIN_SETTINGS_EXP5)) {
    pm <- PLUGIN_METHODS_EXP5[i]
    ps <- PLUGIN_SETTINGS_EXP5[i]
    vn_n0 <- results_df[results_df$n0 == n0v, ]

    for (ct in list(
      list(m1 = pm,         m2 = "OPT-hetero-oracle", label = "plugin - oracle"),
      list(m1 = "OPT-homo", m2 = pm,                  label = "homo - plugin")
    )) {
      pd <- tryCatch(
        paired_diff_summary(vn_n0, ct$m1, ct$m2, "excess_risk"),
        error = function(e) NULL
      )
      if (!is.null(pd) && nrow(pd) > 0L) {
        cat(sprintf("  %-16s %-22s %10.6f %10.6f %10.6f %10.6f\n",
                    ps, ct$label,
                    pd$mean_diff[1L], pd$mcse_diff[1L],
                    pd$lower[1L],     pd$upper[1L]))
      }
    }
  }
}

# ── CSV 저장 ─────────────────────────────────────────────────────────────────

summary_rows <- list()
for (n0v in N0_GRID_EXP5) {
  vn <- valid[valid$n0 == n0v, ]
  for (m in ALL_METHODS_EXP5) {
    sub <- vn[vn$method == m, ]
    if (nrow(sub) == 0L) next
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      method           = m,
      misspec_setting  = sub$misspec_setting[1L],
      n0               = n0v,
      mean_excess_risk = mean(sub$excess_risk, na.rm = TRUE),
      mcse_excess_risk = compute_mcse(sub$excess_risk),
      mean_param_mse   = mean(sub$param_mse,   na.rm = TRUE),
      mean_sigma_cor   = mean(sub$sigma_cor,   na.rm = TRUE),
      mean_score_cor   = mean(sub$score_cor,   na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
}
summary_df  <- do.call(rbind, summary_rows)
csv_summary <- file.path("results", "summary", "exp5_misspec_n0_summary.csv")
write.csv(summary_df, csv_summary, row.names = FALSE)
cat(sprintf("\n저장 완료: %s\n", csv_summary))

# Paired diff CSV: plugin - oracle, homo - plugin × n0 × misspec_setting
paired_rows <- list()
for (n0v in N0_GRID_EXP5) {
  vn_n0 <- results_df[results_df$n0 == n0v, ]
  vn    <- valid[valid$n0 == n0v, ]
  for (i in seq_along(PLUGIN_SETTINGS_EXP5)) {
    pm <- PLUGIN_METHODS_EXP5[i]
    ps <- PLUGIN_SETTINGS_EXP5[i]
    mean_sc <- mean(vn$sigma_cor[vn$method == pm], na.rm = TRUE)
    mean_rc <- mean(vn$score_cor[vn$method == pm], na.rm = TRUE)
    for (ct in list(
      list(m1 = pm,         m2 = "OPT-hetero-oracle", label = "plugin - oracle"),
      list(m1 = "OPT-homo", m2 = pm,                  label = "homo - plugin")
    )) {
      pd <- tryCatch(
        paired_diff_summary(vn_n0, ct$m1, ct$m2, "excess_risk"),
        error = function(e) NULL
      )
      if (!is.null(pd) && nrow(pd) > 0L) {
        paired_rows[[length(paired_rows) + 1L]] <- data.frame(
          misspec_setting = ps,
          n0              = n0v,
          contrast        = ct$label,
          mean_diff       = pd$mean_diff[1L],
          mcse_diff       = pd$mcse_diff[1L],
          lower           = pd$lower[1L],
          upper           = pd$upper[1L],
          mean_sigma_cor  = mean_sc,
          mean_score_cor  = mean_rc,
          stringsAsFactors = FALSE
        )
      }
    }
  }
}
if (length(paired_rows) > 0L) {
  paired_df  <- do.call(rbind, paired_rows)
  csv_paired <- file.path("results", "summary", "exp5_misspec_n0_paired.csv")
  write.csv(paired_df, csv_paired, row.names = FALSE)
  cat(sprintf("저장 완료: %s\n", csv_paired))
}

cat("\nExp5 완료.\n")
