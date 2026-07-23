# 00_check_dgp_diagnostics.R
# DGP 진단: 공변량 분포, sigma 분포, leverage, plugin 추정 품질 확인
# 실행 위치: sim_final/ 폴더에서  source("scripts/00_check_dgp_diagnostics.R")

source("config/config_final.R")
source("R/00_utils.R")
source("R/01_dgp.R")
source("R/02_sampling_methods.R")
source("R/03_variance_estimation.R")

N_DIAG    <- 100000L
SEED_DIAG <- 42L

cat("====================================================\n")
cat(" DGP 진단 (sim_final 확정 DGP)\n")
cat("====================================================\n")
cat(sprintf(" N_diag=%d, seed=%d\n", N_DIAG, SEED_DIAG))
cat(sprintf(" beta = c(%s)\n", paste(round(beta_true, 3), collapse = ", ")))
cat(sprintf(" %s\n\n", dgp_sigma_label))

# ── 이분산 ────────────────────────────────────────────────────────────────────
dat <- generate_data_final(N_DIAG, beta_true, heteroscedastic = TRUE,
                           seed = SEED_DIAG)

cat("── 공변량 분포 ──\n")
for (nm in c("Z1", "Z4", "Z7")) {
  cat(sprintf("  %-4s: mean=%7.4f  var=%7.4f\n",
              nm, mean(dat$X[, nm]), var(dat$X[, nm])))
}

cat("\n── sigma 분포 (이분산) ──\n")
cat(sprintf("  sigma:   min=%.4f, mean=%.4f, max=%.4f, sd=%.4f\n",
            min(dat$sigma_vec), mean(dat$sigma_vec),
            max(dat$sigma_vec), sd(dat$sigma_vec)))
cat(sprintf("  sigma^2: min=%.4f, mean=%.4f, max=%.4f\n",
            min(dat$sigma2_vec), mean(dat$sigma2_vec), max(dat$sigma2_vec)))
cat(sprintf("  cor(sigma2, ell) = %.4f\n", safe_cor(dat$sigma2_vec, dat$ell)))
cat(sprintf("  cor(|s|, sigma2) = %.4f  (s=0.5*Z1+0.5*Z2+Z6+Z8)\n",
            safe_cor(abs(dat$s_var), dat$sigma2_vec)))

cat("\n── Leverage (ell) 분포 ──\n")
cat(sprintf("  min=%.4f, mean=%.4f, max=%.4f\n",
            min(dat$ell), mean(dat$ell), max(dat$ell)))
cat(sprintf("  sum(ell) = %.2f  (이론값 N*p=%.0f)\n", sum(dat$ell), N_DIAG * ncol(dat$X)))

# ── Sampling scores k=600 (확률 기반 method만 -- IBOSS는 별도 블록) ───────────
cat("\n── Sampling scores k=600 (SRS/LEV-IPW/OPT-*) ──\n")
k_diag <- 600L

# plugin score 포함하려면 pilot 한 번 뽑아 sigma2_hat 계산 필요
set.seed(SEED_DIAG)
idx_pilot_diag  <- sample.int(N_DIAG, n0_default, replace = FALSE)
pi_pilot_diag   <- rep(n0_default / N_DIAG, N_DIAG)
beta_pilot_diag <- fit_ipw(dat$X, dat$y, idx_pilot_diag, pi = pi_pilot_diag)
if (is.null(beta_pilot_diag)) beta_pilot_diag <- beta_true

s2hat_param_diag    <- estimate_sigma2_plugin_param(dat, idx_pilot_diag, beta_pilot_diag)
s2hat_nonparam_diag <- estimate_sigma2_plugin_nonparam(dat, idx_pilot_diag, beta_pilot_diag)

scores <- get_scores(dat, s2hat_param_diag, s2hat_nonparam_diag)
for (m in names(scores)) {
  pi_ <- compute_pi(scores[[m]], k_diag)
  cat(sprintf("  %-25s: n_clip=%4d, ESS=%.0f\n",
              m, sum(pi_ >= 1 - 1e-12),
              sum(pi_)^2 / sum(pi_^2)))
}

# ── IBOSS 진단 (결정론적 -- 코드 동작 확인용, 최소 버전) ─────────────────────
cat("\n── IBOSS 선택 진단 k=600 ──\n")
idx_iboss_diag <- select_iboss(dat$X, k_diag)
cat(sprintf("  목표 k=%d, 실제 선택 크기=%d (차이=%d, r 버림에 의한 것)\n",
            k_diag, length(idx_iboss_diag), k_diag - length(idx_iboss_diag)))

# ── Plugin 추정 품질 (N=10000, n0=100, 5회, param vs nonparam) ───────────────
cat("\n── Plugin 추정 품질 (N=10000, n0=200, 5회 반복) ──\n")
dat_s <- generate_data_final(N, beta_true, heteroscedastic = TRUE,
                             seed = SEED_DATA)
cors_sigma_param    <- cors_score_param    <- numeric(5L)
cors_sigma_nonparam <- cors_score_nonparam <- numeric(5L)
for (i in seq_len(5L)) {
  set.seed(SEED_REP + i)
  idx_p  <- sample.int(N, 200, replace = FALSE)
  bp     <- fit_ipw(dat_s$X, dat_s$y, idx_p, rep(200 / N, N))
  if (is.null(bp)) bp <- beta_true

  s2hat_param    <- estimate_sigma2_plugin_param(dat_s, idx_p, bp)
  s2hat_nonparam <- estimate_sigma2_plugin_nonparam(dat_s, idx_p, bp)

  # ── 추가: nonparam 예측이 얼마나 퍼져있는지 확인 ──
  cat(sprintf("  [rep %d] nonparam range=[%.4f, %.4f]  sd=%.4f  (true sigma2 range=[%.4f, %.4f])\n",
              i, min(s2hat_nonparam), max(s2hat_nonparam), sd(s2hat_nonparam),
              min(dat_s$sigma2_vec), max(dat_s$sigma2_vec)))

  cors_sigma_param[i]    <- safe_cor(s2hat_param, dat_s$sigma2_vec)
  cors_score_param[i]    <- safe_cor(sqrt(s2hat_param * dat_s$ell),
                                      sqrt(dat_s$sigma2_vec * dat_s$ell))
  cors_sigma_nonparam[i] <- safe_cor(s2hat_nonparam, dat_s$sigma2_vec)
  cors_score_nonparam[i] <- safe_cor(sqrt(s2hat_nonparam * dat_s$ell),
                                      sqrt(dat_s$sigma2_vec * dat_s$ell))
}
cat(sprintf("  [param]    sigma_cor: %.4f ± %.4f   score_cor: %.4f ± %.4f\n",
            mean(cors_sigma_param), sd(cors_sigma_param),
            mean(cors_score_param), sd(cors_score_param)))
cat(sprintf("  [nonparam] sigma_cor: %.4f ± %.4f   score_cor: %.4f ± %.4f \n",
            mean(cors_sigma_nonparam), sd(cors_sigma_nonparam),
            mean(cors_score_nonparam), sd(cors_score_nonparam)))

# ── 등분산 확인 ───────────────────────────────────────────────────────────────
cat("\n── 등분산 DGP 확인 ──\n")
dat_h <- generate_data_final(N_DIAG, beta_true, heteroscedastic = FALSE,
                             seed = SEED_DIAG)
cat(sprintf("  unique(sigma) 개수: %d  (1이어야 함)\n",
            length(unique(round(dat_h$sigma_vec, 10)))))
cat(sprintf("  var(y) = %.4f\n", var(dat_h$y)))

cat("\n진단 완료\n")