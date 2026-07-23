# 06_runtime_utils.R
# 런타임 측정 유틸리티 (total 측정 + 단계별 분해 측정)
# 의존: 00_utils.R, 01_dgp.R, 02_sampling_methods.R, 03_variance_estimation.R
#
# variance_est는 fit/predict 두 단계로 분리:
#   param   : lm() fit(pilot n0에 의존) vs predict(N에 선형이지만 극히 저렴)
#   nonparam: randomForest() fit(n0에 의존) vs predict(N에 비례 -- N-scaling 실험의 핵심 관찰 대상)
# OPT-hetero-plugin은 -param / -nonparam 두 method로 분리, pilot_draw/pilot_fit/moment_matrix는 공유(동일 절차이므로 1회만 측정)

# ── (구) 전체 rep 경과 시간 -- 필요시 유지, 미사용해도 무방 ──────────────────

time_full_rep <- function(dat, k, X_test, mu_test, y_test, n0) {
  t0 <- proc.time()
  run_one_rep(dat, k, X_test, mu_test, y_test, n0)
  (proc.time() - t0)[["elapsed"]]
}

benchmark_rep <- function(dat, k, X_test, mu_test, y_test, n0,
                          n_bench = 50L, seed_start = 2001L) {
  for (w in 1:3) time_full_rep(dat, k, X_test, mu_test, y_test, n0)
  times <- numeric(n_bench)
  for (i in seq_len(n_bench)) {
    set.seed(seed_start + i)
    times[i] <- time_full_rep(dat, k, X_test, mu_test, y_test, n0)
    if (i %% 10L == 0L)
      cat(sprintf("  bench %d/%d: %.3fs\n", i, n_bench, times[i]))
  }
  list(
    mean    = mean(times),
    sd      = sd(times),
    median  = median(times),
    times   = times,
    n_bench = n_bench
  )
}

# ── (신) 단계별 경과 시간 측정 ───────────────────────────────────────────────
#
# 측정 단계 (8단계):
#   pilot_draw            : 균등 pilot 표본 추출 (param/nonparam 공유)
#   pilot_fit             : pilot IPW-OLS fitting (param/nonparam 공유)
#   variance_est_fit      : sigma^2 추정 모형 학습 (param: lm / nonparam: randomForest)
#   variance_est_predict  : 전체 population에 대한 sigma^2 예측
#   moment_matrix         : X'X 및 leverage score ell 계산 (baseline 포함 전 method 공유)
#   score_compute         : pi_i 계산 및 정규화 (compute_pi)
#   subsample_draw        : Bernoulli sampling
#   final_fit             : weighted OLS fitting
#
# 방법별 해당 단계:
#   FULL                       : final_fit
#   SRS                        : subsample_draw, final_fit
#   LEV-IPW                    : moment_matrix, score_compute, subsample_draw, final_fit
#   OPT-homo                   : moment_matrix, score_compute, subsample_draw, final_fit
#   OPT-hetero-oracle           : moment_matrix, score_compute, subsample_draw, final_fit
#   OPT-hetero-plugin-param    : pilot_draw, pilot_fit, variance_est_fit, variance_est_predict,
#                                moment_matrix, score_compute, subsample_draw, final_fit
#   OPT-hetero-plugin-nonparam : (위와 동일 구조, variance_est만 RF)

# 단일 rep: 전 method 단계별 경과 시간
# 반환: data.frame(method, stage, elapsed)
time_stages_rep <- function(dat, k, n0) {
  N <- nrow(dat$X)
  X <- dat$X
  y <- dat$y

  rows <- list()
  add  <- function(method, stage, elapsed) {
    rows[[length(rows) + 1L]] <<- data.frame(
      method = method, stage = stage, elapsed = elapsed,
      stringsAsFactors = FALSE
    )
  }

  # ── 공유 사전 계산 ──────────────────────────────────────────────────────────

  # moment_matrix: X'X 역행렬 → leverage score (LEV-IPW / OPT-* 공통)
  t0      <- proc.time()
  XtX     <- crossprod(X)
  XtX_inv <- solve(XtX)
  ell_rep <- rowSums((X %*% XtX_inv) * X)
  t_moment <- (proc.time() - t0)[["elapsed"]]

  # Pilot (param/nonparam 공유 -- 동일 절차이므로 1회만 측정)
  t0           <- proc.time()
  idx_pilot    <- sample.int(N, n0, replace = FALSE)
  t_pilot_draw <- (proc.time() - t0)[["elapsed"]]

  pi_pilot <- rep(n0 / N, N)
  t0 <- proc.time()
  beta_pilot <- fit_ipw(X, y, idx_pilot, pi = pi_pilot)
  t_pilot_fit <- (proc.time() - t0)[["elapsed"]]
  if (is.null(beta_pilot)) beta_pilot <- dat$beta

  resid_pilot     <- y[idx_pilot] - as.numeric(X[idx_pilot, , drop = FALSE] %*% beta_pilot)
  abs_resid_pilot <- abs(resid_pilot)

  # ── variance_est: param (fit/predict 분리) ───────────────────────────────────
  t0 <- proc.time()
  abs_s_pilot   <- make_abs_s(X[idx_pilot, , drop = FALSE])
  pilot_df_p    <- data.frame(abs_resid = abs_resid_pilot, abs_s = abs_s_pilot)
  sig_model_p   <- lm(abs_resid ~ abs_s, data = pilot_df_p)
  t_var_fit_param <- (proc.time() - t0)[["elapsed"]]

  t0 <- proc.time()
  sigma_hat_p    <- as.numeric(predict(sig_model_p, newdata = data.frame(abs_s = make_abs_s(X))))
  sigma_hat_p    <- pmax(sigma_hat_p, 1e-10)
  sigma2_hat_param <- (sigma_hat_p / sqrt(2 / pi))^2
  t_var_predict_param <- (proc.time() - t0)[["elapsed"]]

  # ── variance_est: nonparam (fit/predict 분리, RF) ─────────────────────────────
  t0 <- proc.time()
  Z_pilot     <- as.data.frame(X[idx_pilot, -1, drop = FALSE])   # intercept 제외, Z1~Z9
  pilot_df_np <- data.frame(abs_resid = abs_resid_pilot, Z_pilot)
  rf_model    <- randomForest::randomForest(abs_resid ~ ., data = pilot_df_np,
                                            ntree = 100,
                                            mtry = RF_MTRY_FIXED,
                                            nodesize = RF_NODESIZE_FIXED)
  t_var_fit_nonparam <- (proc.time() - t0)[["elapsed"]]

  t0 <- proc.time()
  Z_full     <- as.data.frame(X[, -1, drop = FALSE])
  pred_abs   <- as.numeric(predict(rf_model, newdata = Z_full))
  pred_abs   <- pmax(pred_abs, 1e-10)
  sigma2_hat_nonparam <- (pred_abs / sqrt(2 / pi))^2
  t_var_predict_nonparam <- (proc.time() - t0)[["elapsed"]]

  # ── FULL ────────────────────────────────────────────────────────────────────
  t0 <- proc.time()
  .lm.fit(x = X, y = y)
  add("FULL", "final_fit", (proc.time() - t0)[["elapsed"]])

  # ── SRS ─────────────────────────────────────────────────────────────────────
  pi_srs <- pmin(rep(k / N, N), 1)
  t0 <- proc.time()
  idx_srs <- which(rbinom(N, 1L, pi_srs) == 1L)
  add("SRS", "subsample_draw", (proc.time() - t0)[["elapsed"]])
  t0 <- proc.time()
  fit_ipw(X, y, idx_srs, pi = pi_srs)
  add("SRS", "final_fit", (proc.time() - t0)[["elapsed"]])

  # ── LEV-IPW ─────────────────────────────────────────────────────────────────
  add("LEV-IPW", "moment_matrix", t_moment)
  t0 <- proc.time()
  pi_lev <- compute_pi(ell_rep, k)
  add("LEV-IPW", "score_compute", (proc.time() - t0)[["elapsed"]])
  t0 <- proc.time()
  idx_lev <- which(rbinom(N, 1L, pi_lev) == 1L)
  add("LEV-IPW", "subsample_draw", (proc.time() - t0)[["elapsed"]])
  t0 <- proc.time()
  fit_ipw(X, y, idx_lev, pi = pi_lev)
  add("LEV-IPW", "final_fit", (proc.time() - t0)[["elapsed"]])

  # ── OPT-homo ────────────────────────────────────────────────────────────────
  add("OPT-homo", "moment_matrix", t_moment)
  t0 <- proc.time()
  pi_homo <- compute_pi(sqrt(ell_rep), k)
  add("OPT-homo", "score_compute", (proc.time() - t0)[["elapsed"]])
  t0 <- proc.time()
  idx_homo <- which(rbinom(N, 1L, pi_homo) == 1L)
  add("OPT-homo", "subsample_draw", (proc.time() - t0)[["elapsed"]])
  t0 <- proc.time()
  fit_ipw(X, y, idx_homo, pi = pi_homo)
  add("OPT-homo", "final_fit", (proc.time() - t0)[["elapsed"]])

  # ── OPT-hetero-oracle ───────────────────────────────────────────────────────
  add("OPT-hetero-oracle", "moment_matrix", t_moment)
  t0 <- proc.time()
  oracle_score <- sqrt(dat$sigma2_vec * ell_rep)
  pi_oracle    <- compute_pi(oracle_score, k)
  add("OPT-hetero-oracle", "score_compute", (proc.time() - t0)[["elapsed"]])
  t0 <- proc.time()
  idx_oracle <- which(rbinom(N, 1L, pi_oracle) == 1L)
  add("OPT-hetero-oracle", "subsample_draw", (proc.time() - t0)[["elapsed"]])
  t0 <- proc.time()
  fit_ipw(X, y, idx_oracle, pi = pi_oracle)
  add("OPT-hetero-oracle", "final_fit", (proc.time() - t0)[["elapsed"]])

  # ── OPT-hetero-plugin-param ───────────────────────────────────────────────────
  add("OPT-hetero-plugin-param", "pilot_draw",           t_pilot_draw)
  add("OPT-hetero-plugin-param", "pilot_fit",             t_pilot_fit)
  add("OPT-hetero-plugin-param", "variance_est_fit",      t_var_fit_param)
  add("OPT-hetero-plugin-param", "variance_est_predict",  t_var_predict_param)
  add("OPT-hetero-plugin-param", "moment_matrix",         t_moment)
  t0 <- proc.time()
  score_param <- sqrt(sigma2_hat_param * ell_rep)
  pi_param    <- compute_pi(score_param, k)
  add("OPT-hetero-plugin-param", "score_compute", (proc.time() - t0)[["elapsed"]])
  t0 <- proc.time()
  idx_param <- which(rbinom(N, 1L, pi_param) == 1L)
  add("OPT-hetero-plugin-param", "subsample_draw", (proc.time() - t0)[["elapsed"]])
  t0 <- proc.time()
  fit_ipw(X, y, idx_param, pi = pi_param)
  add("OPT-hetero-plugin-param", "final_fit", (proc.time() - t0)[["elapsed"]])

  # ── OPT-hetero-plugin-nonparam ────────────────────────────────────────────────
  add("OPT-hetero-plugin-nonparam", "pilot_draw",          t_pilot_draw)
  add("OPT-hetero-plugin-nonparam", "pilot_fit",            t_pilot_fit)
  add("OPT-hetero-plugin-nonparam", "variance_est_fit",     t_var_fit_nonparam)
  add("OPT-hetero-plugin-nonparam", "variance_est_predict", t_var_predict_nonparam)
  add("OPT-hetero-plugin-nonparam", "moment_matrix",        t_moment)
  t0 <- proc.time()
  score_nonparam <- sqrt(sigma2_hat_nonparam * ell_rep)
  pi_nonparam    <- compute_pi(score_nonparam, k)
  add("OPT-hetero-plugin-nonparam", "score_compute", (proc.time() - t0)[["elapsed"]])
  t0 <- proc.time()
  idx_nonparam <- which(rbinom(N, 1L, pi_nonparam) == 1L)
  add("OPT-hetero-plugin-nonparam", "subsample_draw", (proc.time() - t0)[["elapsed"]])
  t0 <- proc.time()
  fit_ipw(X, y, idx_nonparam, pi = pi_nonparam)
  add("OPT-hetero-plugin-nonparam", "final_fit", (proc.time() - t0)[["elapsed"]])

  do.call(rbind, rows)
}

# n_bench 회 반복 → 단계별 요약 (워밍업 3회 포함)
benchmark_decomposed <- function(dat, k, n0,
                                 n_bench = 50L, seed_start = 2001L) {
  for (w in seq_len(3L)) {
    set.seed(seed_start - w)
    time_stages_rep(dat, k, n0)
  }

  rep_list <- vector("list", n_bench)
  for (i in seq_len(n_bench)) {
    set.seed(seed_start + i)
    df_i          <- time_stages_rep(dat, k, n0)
    df_i$bench_i  <- i
    rep_list[[i]] <- df_i
    if (i %% 10L == 0L) cat(sprintf("  bench %d/%d\n", i, n_bench))
  }
  raw_df <- do.call(rbind, rep_list)

  keys <- unique(raw_df[, c("method", "stage")])
  summary_rows <- lapply(seq_len(nrow(keys)), function(j) {
    m <- keys$method[j]
    s <- keys$stage[j]
    v <- raw_df$elapsed[raw_df$method == m & raw_df$stage == s]
    data.frame(
      method      = m,
      k           = k,
      stage       = s,
      mean_time   = mean(v),
      sd_time     = sd(v),
      median_time = median(v),
      stringsAsFactors = FALSE
    )
  })

  list(
    summary = do.call(rbind, summary_rows),
    raw     = raw_df
  )
}