# 06_runtime_utils.R
# 런타임 측정 유틸리티 (total 측정 + 단계별 분해 측정)
# 의존: 00_utils.R, 02_sampling_methods.R, 03_variance_estimation.R

# ── (구) 전체 rep 경과 시간 ───────────────────────────────────────────────────

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
# 측정 단계 (OPT-hetero-plugin 기준 7단계):
#   pilot_draw    : 균등 pilot 표본 추출
#   pilot_fit     : pilot IPW-OLS fitting
#   variance_est  : 잔차 기반 sigma^2 추정
#   moment_matrix : X'X 및 leverage score ell 계산
#   score_compute : pi_i 계산 및 정규화 (compute_pi)
#   subsample_draw: Bernoulli sampling
#   final_fit     : weighted OLS fitting
#
# 방법별 해당 단계:
#   FULL               : final_fit
#   SRS                : subsample_draw, final_fit
#   LEV-IPW            : moment_matrix, score_compute, subsample_draw, final_fit
#   OPT-homo           : moment_matrix, score_compute, subsample_draw, final_fit
#   OPT-hetero-oracle  : moment_matrix, score_compute, subsample_draw, final_fit
#   OPT-hetero-plugin  : pilot_draw, pilot_fit, variance_est,
#                        moment_matrix, score_compute, subsample_draw, final_fit

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

  # OPT-hetero-plugin pilot 단계 (sigma2_hat 도 공유)
  t0           <- proc.time()
  idx_pilot    <- sample.int(N, n0, replace = FALSE)
  t_pilot_draw <- (proc.time() - t0)[["elapsed"]]

  pi_pilot <- rep(n0 / N, N)
  t0 <- proc.time()
  beta_pilot <- fit_ipw(X, y, idx_pilot, pi = pi_pilot)
  t_pilot_fit <- (proc.time() - t0)[["elapsed"]]
  if (is.null(beta_pilot)) beta_pilot <- dat$beta

  t0         <- proc.time()
  sigma2_hat <- estimate_sigma2_plugin(dat, idx_pilot, beta_pilot)
  t_var_est  <- (proc.time() - t0)[["elapsed"]]

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
  add("OPT-hetero-oracle", "score_compute",
      (proc.time() - t0)[["elapsed"]])
  t0 <- proc.time()
  idx_oracle <- which(rbinom(N, 1L, pi_oracle) == 1L)
  add("OPT-hetero-oracle", "subsample_draw",
      (proc.time() - t0)[["elapsed"]])
  t0 <- proc.time()
  fit_ipw(X, y, idx_oracle, pi = pi_oracle)
  add("OPT-hetero-oracle", "final_fit",
      (proc.time() - t0)[["elapsed"]])

  # ── OPT-hetero-plugin ───────────────────────────────────────────────────────
  add("OPT-hetero-plugin", "pilot_draw",    t_pilot_draw)
  add("OPT-hetero-plugin", "pilot_fit",     t_pilot_fit)
  add("OPT-hetero-plugin", "variance_est",  t_var_est)
  add("OPT-hetero-plugin", "moment_matrix", t_moment)
  t0 <- proc.time()
  plugin_score <- sqrt(sigma2_hat * ell_rep)
  pi_plugin    <- compute_pi(plugin_score, k)
  add("OPT-hetero-plugin", "score_compute",
      (proc.time() - t0)[["elapsed"]])
  t0 <- proc.time()
  idx_plugin <- which(rbinom(N, 1L, pi_plugin) == 1L)
  add("OPT-hetero-plugin", "subsample_draw",
      (proc.time() - t0)[["elapsed"]])
  t0 <- proc.time()
  fit_ipw(X, y, idx_plugin, pi = pi_plugin)
  add("OPT-hetero-plugin", "final_fit",
      (proc.time() - t0)[["elapsed"]])

  do.call(rbind, rows)
}

# n_bench 회 반복 → 단계별 요약 (워밍업 3회 포함)
# 반환: list(
#   summary : data.frame(method, k, stage, mean_time, sd_time, median_time),
#   raw     : data.frame(method, stage, elapsed, bench_i)
# )
benchmark_decomposed <- function(dat, k, n0,
                                 n_bench = 50L, seed_start = 2001L) {
  # 워밍업 3회 (JIT / 캐시 효과 제거)
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

  # method × stage 집계
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
