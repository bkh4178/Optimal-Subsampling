# 07b_exp5_naive_param_baseline.R
# Exp5b: "naive param" 베이스라인 -- rf-*와 공정하게 비교하기 위한 추가 실험
#
# 문제의식: param-correct/mild/severe는 정확한 결합계수(0.5,0.5,1,1)와 절댓값
# 변환까지 이미 주어진 abs_s를 입력받아 절편+기울기 2개 파라미터만 fitting함.
# 반면 rf-*는 raw Z만 받아 결합방식·비선형성을 전부 스스로 찾아야 함.
# 즉 "변수를 안다"는 말이 param과 rf에서 질적으로 다른 정보량을 의미했음.
#
# 이 스크립트는 rf-*와 "동일한 정보량"(raw Z만, 결합계수/변환 없음)을 param에도
# 줘서 -- 즉 단순 다중선형회귀 abs_resid ~ Z1+Z2+... -- 공정한 비교 지점을 만듦.
#
# param-naive-full   : Z1,Z2,Z6,Z8 (rf-full과 동일 feature)
# param-naive-mild   : Z1,Z2,Z6    (rf-mild와 동일 feature)
# param-naive-severe : Z1,Z2       (rf-severe와 동일 feature)
#
# 07_exp5(new)_misspec_n0_sensitivity.R과 완전히 동일한 dat/dat_test/seed 체계를
# 사용해 pilot(idx_pilot, beta_pilot)이 rep/n0 단위로 bit-identical하게 재현됨
# -- 이후 rbind해서 기존 저장 결과(exp5_misspec_n0.rds)와 paired 비교 가능.
#
# baseline(FULL/SRS/LEV-IPW/OPT-homo/OPT-hetero-oracle)과 param/rf 6종은
# 이미 exp5_misspec_n0.rds에 있으므로 재계산하지 않음 -- naive 3종만 신규 계산.
#
# 실행 위치: sim_final/ 폴더에서  source("scripts/07b_exp5_naive_param_baseline.R")

source("config/config_final.R")
source("R/00_utils.R")
source("R/01_dgp.R")
source("R/02_sampling_methods.R")
source("R/03_variance_estimation.R")
source("R/04_metrics.R")
source("R/05_mcse_summary.R")

DIR_EXP5 <- file.path("results", "raw", "exp5")
ensure_dir(DIR_EXP5)

# ── 07_exp5(new)와 반드시 동일해야 하는 상수 (pilot 재현성 보장) ─────────────
K_EXP5       <- 1000L
N0_GRID_EXP5 <- c(50L, 100L, 200L, 300L, 400L, 500L)

FEATURES_FULL   <- c("Z1", "Z2", "Z6", "Z8")
FEATURES_MILD   <- c("Z1", "Z2", "Z6")
FEATURES_SEVERE <- c("Z1", "Z2")

PLUGIN_SETTINGS_EXP5B <- c("correct_scale", "mild_omitted", "severe_omitted")
PLUGIN_METHODS_EXP5B  <- c("param-naive-full", "param-naive-mild", "param-naive-severe")

# ── naive param 추정 (raw feature 다중선형회귀, 결합계수/변환 없음) ──────────
# 주: sqrt(2/pi) 보정 상수를 안 곱해도 이후 π 정규화(sum pi=k) 과정에서
# 상수배는 자동 상쇄되고, sigma_cor(상관계수)도 양의 스칼라 배에 불변이라
# 결과에 영향 없음 -- param-correct/mild/severe와 동일한 컨벤션 유지.
est_sigma2_param_naive <- function(dat, idx_pilot, beta_pilot, feature_names) {
  X_pilot <- dat$X[idx_pilot, , drop = FALSE]
  resid_pilot <- dat$y[idx_pilot] - as.numeric(X_pilot %*% beta_pilot)

  pilot_df <- data.frame(abs_resid = abs(resid_pilot),
                         X_pilot[, feature_names, drop = FALSE])
  sig_model <- lm(abs_resid ~ ., data = pilot_df)

  newdata   <- as.data.frame(dat$X[, feature_names, drop = FALSE])
  sigma_hat <- as.numeric(predict(sig_model, newdata = newdata))

  bad <- !is.finite(sigma_hat) | sigma_hat <= 0
  if (any(bad)) {
    good_vals      <- sigma_hat[!bad]
    sigma_hat[bad] <- if (length(good_vals) > 0L) median(good_vals) else 1
  }
  pmax(sigma_hat, 1e-10)^2
}

# ── 단일 replicate: naive 3종만 ───────────────────────────────────────────────

run_one_rep_exp5b <- function(dat, k, X_test, mu_test, y_test, mse0, n0) {
  N    <- nrow(dat$X)
  beta <- dat$beta

  # 07_exp5(new)와 완전히 동일한 순서 -- pilot 재현성 보장 (RNG 소비 순서 일치)
  idx_pilot  <- sample.int(N, n0, replace = FALSE)
  beta_pilot <- fit_ipw(dat$X, dat$y, idx_pilot, pi = rep(n0 / N, N))
  if (is.null(beta_pilot)) beta_pilot <- beta

  sigma2_naive_full   <- est_sigma2_param_naive(dat, idx_pilot, beta_pilot, FEATURES_FULL)
  sigma2_naive_mild   <- est_sigma2_param_naive(dat, idx_pilot, beta_pilot, FEATURES_MILD)
  sigma2_naive_severe <- est_sigma2_param_naive(dat, idx_pilot, beta_pilot, FEATURES_SEVERE)

  oracle_score <- sqrt(dat$sigma2_vec * dat$ell)

  make_diag <- function(sigma2_hat) {
    list(
      sigma_cor = safe_cor(sigma2_hat, dat$sigma2_vec),
      score_cor = safe_cor(sqrt(sigma2_hat * dat$ell), oracle_score)
    )
  }
  diag_full   <- make_diag(sigma2_naive_full)
  diag_mild   <- make_diag(sigma2_naive_mild)
  diag_severe <- make_diag(sigma2_naive_severe)

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

    cond_XtWX <- tryCatch({
      if (length(idx_final) < ncol(dat$X)) return(NA_real_)
      w    <- 1 / pmax(pi_raw[idx_final], 1e-10)
      XtWX <- crossprod(dat$X[idx_final, , drop = FALSE] * sqrt(w))
      ev   <- eigen(XtWX, symmetric = TRUE, only.values = TRUE)$values
      ev   <- ev[is.finite(ev) & ev > 0]
      if (length(ev) == 0L) return(NA_real_)
      max(ev) / max(min(ev), 1e-12)
    }, error = function(e) NA_real_)

    er_leading <- mean(dat$sigma2_vec / pmax(pi_raw, 1e-10) * dat$ell)

    if (is.null(beta_hat)) {
      return(data.frame(
        method = method_name, misspec_setting = misspec_setting, fail = TRUE,
        er_leading = er_leading,
        pred_mse = NA_real_, pred_mse_N = NA_real_,
        param_mse = NA_real_, param_mse_N = NA_real_,
        sigma_cor = sigma_cor, score_cor = score_cor,
        n_realized = length(idx_final), n_clip = n_clip,
        ess_pi = ess_pi, max_weight = max_weight, cond_XtWX = cond_XtWX,
        stringsAsFactors = FALSE
      ))
    }

    db            <- beta_hat - beta
    mse_hat       <- mean(as.numeric(y_test - X_test %*% beta_hat)^2)
    pred_mse_raw  <- mse_hat - mse0
    param_mse_raw <- mean(db^2)

    data.frame(
      method = method_name, misspec_setting = misspec_setting, fail = FALSE,
      er_leading = er_leading,
      pred_mse = pred_mse_raw, pred_mse_N = N * pred_mse_raw,
      param_mse = param_mse_raw, param_mse_N = N * param_mse_raw,
      sigma_cor = sigma_cor, score_cor = score_cor,
      n_realized = length(idx_final), n_clip = n_clip,
      ess_pi = ess_pi, max_weight = max_weight, cond_XtWX = cond_XtWX,
      stringsAsFactors = FALSE
    )
  }

  rbind(
    run_method("param-naive-full",   "correct_scale",  sqrt(sigma2_naive_full   * dat$ell),
               sigma_cor = diag_full$sigma_cor,   score_cor = diag_full$score_cor),
    run_method("param-naive-mild",   "mild_omitted",   sqrt(sigma2_naive_mild   * dat$ell),
               sigma_cor = diag_mild$sigma_cor,   score_cor = diag_mild$score_cor),
    run_method("param-naive-severe", "severe_omitted", sqrt(sigma2_naive_severe * dat$ell),
               sigma_cor = diag_severe$sigma_cor, score_cor = diag_severe$score_cor)
  )
}

# ── 데이터 생성 (07_exp5(new)와 완전히 동일 -- 동일 seed) ────────────────────

cat("====================================================\n")
cat(" Exp5b: Naive param 베이스라인 (rf-*와 공정 비교용)\n")
cat("====================================================\n")
cat(sprintf(" N=%d, k=%d, n_rep=%d\n", N, K_EXP5, n_rep))
cat(sprintf(" n0_grid: %s\n", paste(N0_GRID_EXP5, collapse = ", ")))
cat(" param-naive-full/mild/severe: raw Z 다중선형회귀 (결합계수/변환 없음)\n")
cat("====================================================\n\n")

dat      <- generate_data_final(N,      beta_true, heteroscedastic = TRUE, seed = SEED_DATA)
dat_test <- generate_data_final(N_test, beta_true, heteroscedastic = TRUE, seed = SEED_TEST)
X_test   <- dat_test$X
mu_test  <- dat_test$mu
y_test   <- dat_test$y
mse0     <- mean(as.numeric(y_test - mu_test)^2)

all_rows <- list()
t_start  <- proc.time()

for (n0v in N0_GRID_EXP5) {
  cat(sprintf("-- n0 = %d ------------------------------\n", n0v))
  for (r in seq_len(n_rep)) {
    set.seed(SEED_REP + r)
    df_r      <- run_one_rep_exp5b(dat, K_EXP5, X_test, mu_test, y_test, mse0, n0 = n0v)
    df_r$rep  <- r
    df_r$k    <- K_EXP5
    df_r$n0   <- n0v
    all_rows[[length(all_rows) + 1L]] <- df_r
    if (r %% 100L == 0L) cat(sprintf("  %d / %d\n", r, n_rep))
  }
}

results_df_b <- do.call(rbind, all_rows)
results_df_b <- results_df_b[, c("method", "misspec_setting", "n0", "rep", "k",
                                  "fail", "er_leading",
                                  "pred_mse", "pred_mse_N",
                                  "param_mse", "param_mse_N",
                                  "sigma_cor", "score_cor",
                                  "n_realized", "n_clip", "ess_pi", "max_weight",
                                  "cond_XtWX")]

rds_path_b <- file.path(DIR_EXP5, "exp5b_naive_param.rds")
saveRDS(results_df_b, rds_path_b)
cat(sprintf("\n저장 완료: %s (%d행, %.1f초)\n",
            rds_path_b, nrow(results_df_b), (proc.time() - t_start)[["elapsed"]]))

# ── naive 단독 요약 (3-metric 전체) ──────────────────────────────────────────
valid_b <- results_df_b[!results_df_b$fail, ]

for (n0v in N0_GRID_EXP5) {
  cat(sprintf("\n=== n0 = %d ===\n", n0v))
  vn <- valid_b[valid_b$n0 == n0v, ]
  cat(sprintf("  %-20s %12s %12s %12s %12s %12s %12s\n",
              "method", "er_leading", "pred_mse_N", "mcse(pred)",
              "param_mse_N", "mcse(param)", ""))
  cat(sprintf("  %-20s %12s %12s %12s %12s %12s %12s\n",
              "", "", "sigma_cor", "score_cor", "", "", ""))
  for (m in PLUGIN_METHODS_EXP5B) {
    sub <- vn[vn$method == m, ]
    if (nrow(sub) == 0L) next
    cat(sprintf("  %-20s %12.4f %12.4f %12.4f %12.4f %12.4f\n",
                m,
                mean(sub$er_leading, na.rm = TRUE),
                mean(sub$pred_mse_N, na.rm = TRUE),
                compute_mcse(sub$pred_mse_N),
                mean(sub$param_mse_N, na.rm = TRUE),
                compute_mcse(sub$param_mse_N)))
    cat(sprintf("  %-20s %12s (sigma_cor=%.4f, score_cor=%.4f)\n",
                "", "",
                mean(sub$sigma_cor, na.rm = TRUE),
                mean(sub$score_cor, na.rm = TRUE)))
  }
}

# ── 핵심 비교: param-naive-* vs rf-* (기존 exp5 결과와 병합, 3-metric 전체) ──
rds_path_a <- file.path(DIR_EXP5, "exp5_misspec_n0.rds")
if (file.exists(rds_path_a)) {
  results_df_a <- readRDS(rds_path_a)
  combined     <- rbind(results_df_a, results_df_b)
  valid_a      <- results_df_a[!results_df_a$fail, ]

  cat(sprintf("\n%s\n", strrep("=", 90)))
  cat(" 핵심 비교: param-naive-* vs rf-* (동일 feature, 공정 비교, 3-metric 전체)\n")
  cat(sprintf("%s\n", strrep("=", 90)))

  pairs <- list(
    list(naive = "param-naive-full",   rf = "rf-full",   label = "full (Z1,Z2,Z6,Z8)"),
    list(naive = "param-naive-mild",   rf = "rf-mild",   label = "mild (Z1,Z2,Z6)"),
    list(naive = "param-naive-severe", rf = "rf-severe", label = "severe (Z1,Z2)")
  )
  metric_list <- c("er_leading", "pred_mse_N", "param_mse_N")

  for (n0v in N0_GRID_EXP5) {
    cat(sprintf("\n-- n0 = %d --\n", n0v))
    cn <- combined[combined$n0 == n0v, ]

    for (p in pairs) {
      sc_naive <- mean(valid_b$sigma_cor[valid_b$method == p$naive & valid_b$n0 == n0v], na.rm = TRUE)
      rc_naive <- mean(valid_b$score_cor[valid_b$method == p$naive & valid_b$n0 == n0v], na.rm = TRUE)
      sc_rf    <- mean(valid_a$sigma_cor[valid_a$method == p$rf    & valid_a$n0 == n0v], na.rm = TRUE)
      rc_rf    <- mean(valid_a$score_cor[valid_a$method == p$rf    & valid_a$n0 == n0v], na.rm = TRUE)

      cat(sprintf("  [%s]  sigma_cor: naive=%.4f rf=%.4f   score_cor: naive=%.4f rf=%.4f\n",
                  p$label, sc_naive, sc_rf, rc_naive, rc_rf))

      # 3-metric 각각의 naive - rf paired diff (rep 단위로 매칭)
      for (met in metric_list) {
        vec_naive <- valid_b[valid_b$method == p$naive & valid_b$n0 == n0v, c("rep", met)]
        vec_rf    <- valid_a[valid_a$method == p$rf    & valid_a$n0 == n0v, c("rep", met)]
        mrg <- merge(vec_naive, vec_rf, by = "rep", suffixes = c("_naive", "_rf"))
        if (nrow(mrg) == 0L) next
        met_naive_col <- paste0(met, "_naive")
        met_rf_col    <- paste0(met, "_rf")
        d      <- mrg[[met_naive_col]] - mrg[[met_rf_col]]
        d      <- d[is.finite(d)]
        if (length(d) < 3L) next
        m_d    <- mean(d)
        se_d   <- sd(d) / sqrt(length(d))
        cat(sprintf("      %-12s naive-rf diff: %10.4f  95%% CI [%10.4f, %10.4f]\n",
                    met, m_d, m_d - 1.96 * se_d, m_d + 1.96 * se_d))
      }
    }
  }
} else {
  cat("\n주의: exp5_misspec_n0.rds 없음 -- naive vs rf 비교는 별도로 수행 필요.\n")
}

cat("\nExp5b 완료.\n")