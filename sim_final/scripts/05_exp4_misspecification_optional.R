# 05_exp4_misspecification_optional.R
# Exp4 (선택): Plugin 분산 모형 오명세 실험
#
# Plugin 모형 4종:
#   correct_scale   : abs_resid ~ abs_s_full (sigma-scale, Z1,Z2,Z6,Z8)
#   variance_linear : resid^2   ~ abs_s_full (variance-scale, linear only)
#   mild_omitted    : abs_resid ~ abs_s_mild (sigma-scale, Z8 누락)
#   severe_omitted  : abs_resid ~ abs_s_severe (sigma-scale, Z6+Z8 누락)
# 베이스라인 5종: FULL, SRS, LEV-IPW, OPT-homo, OPT-hetero-oracle
# 실행 위치: sim_final/ 폴더에서  source("scripts/05_exp4_misspecification_optional.R")

source("config/config_final.R")
source("R/00_utils.R")
source("R/01_dgp.R")
source("R/02_sampling_methods.R")
source("R/03_variance_estimation.R")
source("R/04_metrics.R")
source("R/05_mcse_summary.R")

ensure_dir(DIR_EXP4)
ensure_dir(file.path("results", "summary"))
ensure_dir(file.path("results", "tables"))

# ── abs_s 보조 함수 (오명세 드라이버 변수) ───────────────────────────────────
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

# ── Plugin 추정 함수 4종 ─────────────────────────────────────────────────────

# 공통 fallback 처리 (양수 벡터 반환)
.clip_positive <- function(v) {
  bad <- !is.finite(v) | v <= 0
  if (any(bad)) {
    good_vals <- v[!bad]
    v[bad]    <- if (length(good_vals) > 0L) median(good_vals) else 1
  }
  pmax(v, 1e-10)
}

# (1) Sigma-scale 정명세: abs_resid ~ abs_s_full → sigma_hat^2
est_sigma2_correct_scale <- function(dat, idx_pilot, beta_pilot) {
  resid_pilot     <- dat$y[idx_pilot] -
                     as.numeric(dat$X[idx_pilot, , drop = FALSE] %*% beta_pilot)
  abs_s_pilot     <- make_abs_s(dat$X[idx_pilot, , drop = FALSE])
  pilot_df        <- data.frame(abs_resid = abs(resid_pilot), abs_s = abs_s_pilot)
  sig_model       <- lm(abs_resid ~ abs_s, data = pilot_df)
  sigma_hat       <- as.numeric(predict(sig_model,
                       newdata = data.frame(abs_s = make_abs_s(dat$X))))
  .clip_positive(sigma_hat)^2
}

# (2) Variance-scale 선형: resid^2 ~ abs_s_full (quadratic 없음)
est_sigma2_variance_linear <- function(dat, idx_pilot, beta_pilot) {
  resid_pilot  <- dat$y[idx_pilot] -
                  as.numeric(dat$X[idx_pilot, , drop = FALSE] %*% beta_pilot)
  abs_s_pilot  <- make_abs_s(dat$X[idx_pilot, , drop = FALSE])
  pilot_df     <- data.frame(resid2 = resid_pilot^2, abs_s = abs_s_pilot)
  var_model    <- lm(resid2 ~ abs_s, data = pilot_df)
  sigma2_hat   <- as.numeric(predict(var_model,
                    newdata = data.frame(abs_s = make_abs_s(dat$X))))
  .clip_positive(sigma2_hat)
}

# (3) Mild 오명세: abs_resid ~ abs_s_mild (Z8 누락)
est_sigma2_mild_omitted <- function(dat, idx_pilot, beta_pilot) {
  resid_pilot     <- dat$y[idx_pilot] -
                     as.numeric(dat$X[idx_pilot, , drop = FALSE] %*% beta_pilot)
  abs_s_pilot     <- make_abs_s_mild(dat$X[idx_pilot, , drop = FALSE])
  pilot_df        <- data.frame(abs_resid = abs(resid_pilot), abs_s = abs_s_pilot)
  sig_model       <- lm(abs_resid ~ abs_s, data = pilot_df)
  sigma_hat       <- as.numeric(predict(sig_model,
                       newdata = data.frame(abs_s = make_abs_s_mild(dat$X))))
  .clip_positive(sigma_hat)^2
}

# (4) Severe 오명세: abs_resid ~ abs_s_severe (Z6, Z8 누락)
est_sigma2_severe_omitted <- function(dat, idx_pilot, beta_pilot) {
  resid_pilot     <- dat$y[idx_pilot] -
                     as.numeric(dat$X[idx_pilot, , drop = FALSE] %*% beta_pilot)
  abs_s_pilot     <- make_abs_s_severe(dat$X[idx_pilot, , drop = FALSE])
  pilot_df        <- data.frame(abs_resid = abs(resid_pilot), abs_s = abs_s_pilot)
  sig_model       <- lm(abs_resid ~ abs_s, data = pilot_df)
  sigma_hat       <- as.numeric(predict(sig_model,
                       newdata = data.frame(abs_s = make_abs_s_severe(dat$X))))
  .clip_positive(sigma_hat)^2
}

# ── 단일 replicate: 9 methods (4 plugin + 5 baseline) ────────────────────────

run_one_rep_exp4 <- function(dat, k, X_test, mu_test, y_test, n0 = n0_default) {
  N    <- nrow(dat$X)
  beta <- dat$beta

  # 공유 Pilot
  idx_pilot  <- sample.int(N, n0, replace = FALSE)
  beta_pilot <- fit_ipw(dat$X, dat$y, idx_pilot, pi = rep(n0 / N, N))
  if (is.null(beta_pilot)) beta_pilot <- beta

  # 4종 plugin sigma^2 추정
  sigma2_correct   <- est_sigma2_correct_scale(dat, idx_pilot, beta_pilot)
  sigma2_varlinear <- est_sigma2_variance_linear(dat, idx_pilot, beta_pilot)
  sigma2_mild      <- est_sigma2_mild_omitted(dat, idx_pilot, beta_pilot)
  sigma2_severe    <- est_sigma2_severe_omitted(dat, idx_pilot, beta_pilot)

  oracle_score <- sqrt(dat$sigma2_vec * dat$ell)

  make_diag <- function(sigma2_hat) {
    list(
      sigma_cor = safe_cor(sigma2_hat, dat$sigma2_vec),
      score_cor = safe_cor(sqrt(sigma2_hat * dat$ell), oracle_score)
    )
  }
  diag_correct   <- make_diag(sigma2_correct)
  diag_varlinear <- make_diag(sigma2_varlinear)
  diag_mild      <- make_diag(sigma2_mild)
  diag_severe    <- make_diag(sigma2_severe)

  # 서브샘플 method 실행 (완전 컬럼 스키마)
  run_method <- function(method_name, misspec_setting, score,
                         sigma_cor = NA_real_, score_cor = NA_real_) {
    samp      <- sample_bernoulli(score, k)
    idx_final <- samp$idx
    pi_raw    <- samp$pi
    n_clip    <- sum(pi_raw >= 1 - 1e-12)

    nz         <- pi_raw[pi_raw > 0]
    ess_pi     <- if (length(nz) > 0L) sum(pi_raw)^2 / sum(pi_raw^2) else NA_real_
    max_weight <- if (length(nz) > 0L) max(1 / pmax(nz, 1e-10)) else NA_real_

    beta_hat <- fit_ipw(dat$X, dat$y, idx_final, pi = pi_raw)

    cond_XtWX <- tryCatch({
      if (length(idx_final) < ncol(dat$X)) return(NA_real_)
      w    <- 1 / pmax(pi_raw[idx_final], 1e-10)
      XtWX <- crossprod(dat$X[idx_final, , drop = FALSE] * sqrt(w))
      ev   <- eigen(XtWX, symmetric = TRUE, only.values = TRUE)$values
      ev   <- ev[is.finite(ev) & ev > 0]
      if (length(ev) == 0L) return(NA_real_)
      max(ev) / max(min(ev), 1e-12)
    }, error = function(e) NA_real_)

    if (is.null(beta_hat)) {
      return(data.frame(
        method = method_name, misspec_setting = misspec_setting, fail = TRUE,
        excess_risk = NA_real_, pred_mse = NA_real_, param_mse = NA_real_,
        c1_pi = NA_real_, ess_pi = ess_pi, max_weight = max_weight,
        cond_XtWX = cond_XtWX, n_realized = length(idx_final), n_clip = n_clip,
        sigma_cor = sigma_cor, score_cor = score_cor,
        stringsAsFactors = FALSE
      ))
    }

    db <- beta_hat - beta
    data.frame(
      method          = method_name,
      misspec_setting = misspec_setting,
      fail            = FALSE,
      excess_risk     = as.numeric(t(db) %*% dat$A_hat %*% db),
      pred_mse        = mean(as.numeric(mu_test - X_test %*% beta_hat)^2),
      param_mse       = mean(db^2),
      c1_pi           = mean(dat$sigma2_vec * dat$ell / pmax(pi_raw, 1e-10)),
      ess_pi          = ess_pi,
      max_weight      = max_weight,
      cond_XtWX       = cond_XtWX,
      n_realized      = length(idx_final),
      n_clip          = n_clip,
      sigma_cor       = sigma_cor,
      score_cor       = score_cor,
      stringsAsFactors = FALSE
    )
  }

  # FULL (전체 OLS)
  full_samp <- sample_full(dat)
  bhat_full <- full_samp$beta_hat
  full_result <- if (is.null(bhat_full)) {
    data.frame(
      method = "FULL", misspec_setting = "baseline", fail = TRUE,
      excess_risk = NA_real_, pred_mse = NA_real_, param_mse = NA_real_,
      c1_pi = NA_real_, ess_pi = NA_real_, max_weight = NA_real_,
      cond_XtWX = NA_real_, n_realized = full_samp$n_realized, n_clip = 0L,
      sigma_cor = NA_real_, score_cor = NA_real_,
      stringsAsFactors = FALSE
    )
  } else {
    db <- bhat_full - beta
    data.frame(
      method          = "FULL",
      misspec_setting = "baseline",
      fail            = FALSE,
      excess_risk     = as.numeric(t(db) %*% dat$A_hat %*% db),
      pred_mse        = mean(as.numeric(mu_test - X_test %*% bhat_full)^2),
      param_mse       = mean(db^2),
      c1_pi           = NA_real_,
      ess_pi          = NA_real_,
      max_weight      = NA_real_,
      cond_XtWX       = NA_real_,
      n_realized      = full_samp$n_realized,
      n_clip          = 0L,
      sigma_cor       = NA_real_,
      score_cor       = NA_real_,
      stringsAsFactors = FALSE
    )
  }

  rbind(
    full_result,
    run_method("SRS",               "baseline",        rep(1, N)),
    run_method("LEV-IPW",           "baseline",        dat$ell),
    run_method("OPT-homo",          "baseline",        sqrt(dat$ell)),
    run_method("OPT-hetero-oracle", "baseline",        sqrt(dat$sigma2_vec * dat$ell)),
    run_method("plugin-correct",    "correct_scale",   sqrt(sigma2_correct   * dat$ell),
               sigma_cor = diag_correct$sigma_cor,   score_cor = diag_correct$score_cor),
    run_method("plugin-varlinear",  "variance_linear", sqrt(sigma2_varlinear * dat$ell),
               sigma_cor = diag_varlinear$sigma_cor, score_cor = diag_varlinear$score_cor),
    run_method("plugin-mild",       "mild_omitted",    sqrt(sigma2_mild      * dat$ell),
               sigma_cor = diag_mild$sigma_cor,      score_cor = diag_mild$score_cor),
    run_method("plugin-severe",     "severe_omitted",  sqrt(sigma2_severe    * dat$ell),
               sigma_cor = diag_severe$sigma_cor,    score_cor = diag_severe$score_cor)
  )
}

# ── 실험 파라미터 ─────────────────────────────────────────────────────────────

K_EXP4 <- 1000L

PLUGIN_SETTINGS <- c("correct_scale", "variance_linear", "mild_omitted", "severe_omitted")
PLUGIN_METHODS  <- c("plugin-correct", "plugin-varlinear", "plugin-mild", "plugin-severe")
BASE_METHODS    <- c("FULL", "SRS", "LEV-IPW", "OPT-homo", "OPT-hetero-oracle")
ALL_METHODS_EXP4 <- c(BASE_METHODS, PLUGIN_METHODS)

cat("====================================================\n")
cat(" Exp4: Plugin 오명세 실험\n")
cat("====================================================\n")
cat(sprintf(" N=%d, k=%d, n0=%d, n_rep=%d\n", N, K_EXP4, n0_default, n_rep))
cat(sprintf(" %s\n", dgp_sigma_label))
cat(" Plugin 4종: correct_scale / variance_linear / mild_omitted / severe_omitted\n")
cat("====================================================\n\n")

dat     <- generate_data_final(N,      beta_true, heteroscedastic = TRUE, seed = SEED_DATA)
dat_test <- generate_data_final(N_test, beta_true, heteroscedastic = TRUE, seed = SEED_TEST)
X_test  <- dat_test$X
mu_test <- dat_test$mu
y_test  <- dat_test$y

all_rows <- list()
t_start  <- proc.time()

for (r in seq_len(n_rep)) {
  set.seed(SEED_REP + r)
  df_r      <- run_one_rep_exp4(dat, K_EXP4, X_test, mu_test, y_test, n0 = n0_default)
  df_r$rep  <- r
  df_r$k    <- K_EXP4
  all_rows[[r]] <- df_r
  if (r %% 100L == 0L) cat(sprintf("  %d / %d\n", r, n_rep))
}

results_df <- do.call(rbind, all_rows)
rds_path   <- file.path(DIR_EXP4, "exp4_misspecification.rds")
saveRDS(results_df, rds_path)
cat(sprintf("\n저장 완료: %s (%d행, %.1f초)\n",
            rds_path, nrow(results_df), (proc.time() - t_start)[["elapsed"]]))

# ── Method 요약: excess_risk ──────────────────────────────────────────────────
valid <- results_df[!results_df$fail, ]

cat(sprintf("\n%s\n", strrep("=", 70)))
cat(" excess_risk 비교\n")
cat(sprintf("%s\n", strrep("=", 70)))
cat(sprintf("  %-25s %12s %12s\n", "method", "mean", "mcse"))
cat(strrep("-", 52), "\n")
for (m in ALL_METHODS_EXP4) {
  sub <- valid[valid$method == m, ]
  if (nrow(sub) == 0L) next
  cat(sprintf("  %-25s %12.6f %12.6f\n",
              m, mean(sub$excess_risk, na.rm = TRUE),
              compute_mcse(sub$excess_risk)))
}

# ── Plugin 진단: sigma_cor / score_cor ───────────────────────────────────────
cat(sprintf("\n%s\n", strrep("=", 70)))
cat(" Plugin 진단: sigma_cor / score_cor\n")
cat(sprintf("%s\n", strrep("=", 70)))
cat(sprintf("  %-20s %-20s %12s %12s\n",
            "misspec_setting", "method", "sigma_cor", "score_cor"))
cat(strrep("-", 66), "\n")
for (i in seq_along(PLUGIN_SETTINGS)) {
  sub <- valid[valid$method == PLUGIN_METHODS[i], ]
  cat(sprintf("  %-20s %-20s %12.4f %12.4f\n",
              PLUGIN_SETTINGS[i], PLUGIN_METHODS[i],
              mean(sub$sigma_cor, na.rm = TRUE),
              mean(sub$score_cor, na.rm = TRUE)))
}

# ── Paired 차이: 3 contrasts × 4 plugin settings ─────────────────────────────
cat(sprintf("\n%s\n", strrep("=", 70)))
cat(" Paired 차이 (excess_risk 기준): 3 contrasts × 4 plugin 설정\n")
cat(sprintf("%s\n", strrep("=", 70)))
cat(sprintf("  %-16s %-22s %10s %10s %10s %10s %10s %10s\n",
            "setting", "contrast",
            "mean_diff", "mcse_diff", "lower", "upper",
            "sigma_cor", "score_cor"))
cat(strrep("-", 90), "\n")

paired_rows <- list()

for (i in seq_along(PLUGIN_SETTINGS)) {
  ps <- PLUGIN_SETTINGS[i]
  pm <- PLUGIN_METHODS[i]

  mean_sc <- mean(valid$sigma_cor[valid$method == pm], na.rm = TRUE)
  mean_rc <- mean(valid$score_cor[valid$method == pm], na.rm = TRUE)

  contrasts_list <- list(
    list(m1 = "SRS",               m2 = pm, label = "SRS - plugin"),
    list(m1 = "OPT-homo",          m2 = pm, label = "homo - plugin"),
    list(m1 = pm, m2 = "OPT-hetero-oracle", label = "plugin - oracle")
  )

  for (ct in contrasts_list) {
    pd <- tryCatch(
      paired_diff_summary(results_df, ct$m1, ct$m2, "excess_risk"),
      error = function(e) NULL
    )
    if (!is.null(pd) && nrow(pd) > 0L) {
      cat(sprintf("  %-16s %-22s %10.6f %10.6f %10.6f %10.6f %10.4f %10.4f\n",
                  ps, ct$label,
                  pd$mean_diff[1L], pd$mcse_diff[1L],
                  pd$lower[1L],     pd$upper[1L],
                  mean_sc,          mean_rc))
      paired_rows[[length(paired_rows) + 1L]] <- data.frame(
        misspec_setting = ps,
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

# ── CSV 저장 ─────────────────────────────────────────────────────────────────

# Method 요약 CSV
summary_rows <- lapply(ALL_METHODS_EXP4, function(m) {
  sub <- valid[valid$method == m, ]
  if (nrow(sub) == 0L) return(NULL)
  data.frame(
    method           = m,
    misspec_setting  = sub$misspec_setting[1L],
    mean_excess_risk = mean(sub$excess_risk, na.rm = TRUE),
    mcse_excess_risk = compute_mcse(sub$excess_risk),
    mean_pred_mse    = mean(sub$pred_mse,    na.rm = TRUE),
    mean_sigma_cor   = mean(sub$sigma_cor,   na.rm = TRUE),
    mean_score_cor   = mean(sub$score_cor,   na.rm = TRUE),
    stringsAsFactors = FALSE
  )
})
summary_df  <- do.call(rbind, Filter(Negate(is.null), summary_rows))
csv_summary <- file.path("results", "summary", "exp4_misspec_method_summary.csv")
write.csv(summary_df, csv_summary, row.names = FALSE)
cat(sprintf("\n저장 완료: %s\n", csv_summary))

# Paired MCSE CSV
if (length(paired_rows) > 0L) {
  paired_df  <- do.call(rbind, paired_rows)
  csv_paired <- file.path("results", "tables", "exp4_misspec_paired_mcse.csv")
  write.csv(paired_df, csv_paired, row.names = FALSE)
  cat(sprintf("저장 완료: %s\n", csv_paired))
}
