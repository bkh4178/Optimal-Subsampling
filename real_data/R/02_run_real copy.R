# =============================================================================
# 02_run_real.R
# NYC Yellow Taxi — Repeated subsampling experiment (param + RF nonparam plugin)
# 실행 위치: real_data/R/
# =============================================================================

source("../simulation_2plugin/R/00_utils.R") # compute_pi, fit_ipw
source("config_real.R")               # k_grid, n0, B, seed, result_path, rf_mtry, rf_nodesize
library(randomForest)
library(parallel)

# =============================================================================
# 0. 준비 데이터 로드
# =============================================================================
dat <- readRDS("data/taxi_prepared.rds")

X_train   <- dat$X_train
y_train   <- dat$y_train
X_test    <- dat$X_test
y_test    <- dat$y_test
beta_full <- dat$beta_full
A_hat     <- dat$A_hat
N_train   <- dat$N_train

p1 <- ncol(X_train)   # intercept 포함
N_CORES <- parallel::detectCores()

cat(sprintf("N_train = %d / N_test = %d / p+1 = %d\n",
            N_train, length(y_test), p1))
cat(sprintf("k_grid : %s\n", paste(k_grid, collapse = ", ")))
cat(sprintf("n0 = %d / B = %d / RF mtry=%d nodesize=%d / cores=%d\n\n",
            n0, B, rf_mtry, rf_nodesize, N_CORES))

# =============================================================================
# 1. 전체 데이터 leverage (한 번만 계산)
# =============================================================================
A_hat_inv <- solve(A_hat)
ell       <- rowSums((X_train %*% A_hat_inv) * X_train)

# =============================================================================
# 2. Plugin sigma 추정 함수
# =============================================================================

# 병렬 predict (RF predict가 N_train 전체에서 병목 -> 청킹)
predict_rf_parallel <- function(fit, X, n_chunks = N_CORES) {
  N <- nrow(X)
  idx_split <- parallel::splitIndices(N, n_chunks)
  preds <- parallel::mclapply(idx_split,
                               function(idx) predict(fit, X[idx, , drop = FALSE]),
                               mc.cores = n_chunks)
  result <- numeric(N)
  for (i in seq_along(idx_split)) result[idx_split[[i]]] <- preds[[i]]
  unname(result)
}

# param: |resid| ~ trip_duration(열 3) 로그선형, 실제 함수형태 모르니 근사
estimate_sigma_plugin_param <- function(X_pilot, y_pilot, X_full) {
  beta_pilot <- tryCatch(
    as.vector(solve(crossprod(X_pilot), crossprod(X_pilot, y_pilot))),
    error = function(e) NULL
  )
  if (is.null(beta_pilot)) return(rep(1, nrow(X_full)))

  abs_resid <- abs(as.vector(y_pilot - X_pilot %*% beta_pilot))

  X_var   <- cbind(1, X_pilot[, 3])
  fit_var <- tryCatch(lm.fit(X_var, log(abs_resid + 1e-6)), error = function(e) NULL)
  if (is.null(fit_var)) return(rep(1, nrow(X_full)))

  X_var_full <- cbind(1, X_full[, 3])
  pred_abs   <- exp(as.vector(X_var_full %*% coef(fit_var)))

  sigma_hat  <- pmax(pred_abs, 1e-10) / sqrt(2 / pi)   # param과 동일한 정규분포 보정
  sigma_hat^2                                          # σ² 반환
}

# nonparam: RF regression of abs_resid on covariates
# ⚠ CORR_FACTOR는 simulation_2plugin/R/03_variance_estimation.R의
#   estimate_sigma2_plugin_nonparam()이 실제 쓰는 abs_resid->sigma 보정값과
#   반드시 맞춰야 함 — 지금은 정규분포 가정 하 이론값(sqrt(pi/2))으로 임시 세팅

estimate_sigma_plugin_nonparam <- function(X_pilot, y_pilot, X_full) {
  beta_pilot <- tryCatch(
    as.vector(solve(crossprod(X_pilot), crossprod(X_pilot, y_pilot))),
    error = function(e) NULL
  )
  if (is.null(beta_pilot)) return(rep(1, nrow(X_full)))

  abs_resid <- abs(as.vector(y_pilot - X_pilot %*% beta_pilot))

  fit <- tryCatch(
    randomForest(
      x = X_pilot[, -1, drop = FALSE],   # intercept(열1) 제외
      y = abs_resid,
      mtry = rf_mtry, nodesize = rf_nodesize
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) return(rep(1, nrow(X_full)))

  pred_abs_resid <- predict_rf_parallel(fit, X_full[, -1, drop = FALSE])
  sigma_hat <- pmax(pred_abs_resid, 1e-10) / sqrt(2 / pi)   # param과 동일한 정규분포 보정
  sigma_hat^2                                                # σ² 반환 (nonparam과 스케일 통일)
}

# =============================================================================
# 3. rep 단위: pilot 뽑고 sigma_hat 계산 (k와 무관 -> b당 1회만 계산)
# =============================================================================
compute_pilot_scores <- function(seed_rep) {
  set.seed(seed_rep)
  idx_pilot <- sample(N_train, n0)

  t0 <- proc.time()
  sigma_hat_param <- estimate_sigma_plugin_param(
    X_train[idx_pilot, ], y_train[idx_pilot], X_train
  )
  t_param <- (proc.time() - t0)[["elapsed"]]

  t0 <- proc.time()
  sigma_hat_nonparam <- estimate_sigma_plugin_nonparam(
    X_train[idx_pilot, ], y_train[idx_pilot], X_train
  )
  t_nonparam <- (proc.time() - t0)[["elapsed"]]

  list(
    score_uni      = rep(1, N_train),
    score_lev      = ell,
    score_homo     = sqrt(ell),
    score_param    = sqrt(sigma_hat_param    * ell),
    score_nonparam = sqrt(sigma_hat_nonparam * ell),
    diag = data.frame(
      sigma_param_mean    = mean(sigma_hat_param),    sigma_param_sd    = sd(sigma_hat_param),
      sigma_nonparam_mean = mean(sigma_hat_nonparam), sigma_nonparam_sd = sd(sigma_hat_nonparam),
      t_param = t_param, t_nonparam = t_nonparam
    )
  )
}

# =============================================================================
# 4. 주어진 score로 한 번 sampling + estimation
# =============================================================================
run_one_method <- function(name, score, k) {
  pi  <- compute_pi(score, k)
  idx <- which(rbinom(N_train, 1L, pi) == 1L)
  clip_rate <- mean(pi >= 1 - 1e-9)

  if (length(idx) < p1 + 1L) {
    return(data.frame(method = name, k = k, fail = TRUE,
                       pseudo_er = NA_real_, test_mse = NA_real_,
                       n_sub = length(idx), clip_rate = clip_rate))
  }

  beta_hat <- fit_ipw(X_train, y_train, idx, pi)
  if (is.null(beta_hat)) {
    return(data.frame(method = name, k = k, fail = TRUE,
                       pseudo_er = NA_real_, test_mse = NA_real_,
                       n_sub = length(idx), clip_rate = clip_rate))
  }

  diff      <- beta_hat - beta_full
  pseudo_er <- as.numeric(t(diff) %*% A_hat %*% diff)
  test_mse  <- mean((y_test - X_test %*% beta_hat)^2)

  data.frame(method = name, k = k, fail = FALSE,
             pseudo_er = pseudo_er, test_mse = test_mse,
             n_sub = length(idx), clip_rate = clip_rate)
}

# =============================================================================
# 5. 메인 루프: b 바깥 / k 안쪽 (pilot·sigma_hat을 k_grid 전체에서 재사용)
# =============================================================================
set.seed(seed)
all_rows  <- vector("list", B)
diag_rows <- vector("list", B)

cat(sprintf("총 %d reps × %d k = %d runs\n\n", B, length(k_grid), B * length(k_grid)))
t_start <- proc.time()

for (b in seq_len(B)) {
  seed_rep <- seed * 10000L + b
  scores   <- compute_pilot_scores(seed_rep)
  diag_rows[[b]] <- cbind(rep = b, scores$diag)

  methods_list <- list(
    list(name = "SRS",                       score = scores$score_uni),
    list(name = "LEV-IPW",                   score = scores$score_lev),
    list(name = "OPT-homo",                  score = scores$score_homo),
    list(name = "OPT-hetero-plugin-param",    score = scores$score_param),
    list(name = "OPT-hetero-plugin-nonparam", score = scores$score_nonparam)
  )

  rows_b <- vector("list", length(methods_list) * length(k_grid))
  idx_row <- 1L
  for (k in k_grid) {
    for (m in methods_list) {
      rows_b[[idx_row]] <- run_one_method(m$name, m$score, k)
      idx_row <- idx_row + 1L
    }
  }
  all_rows[[b]] <- cbind(rep = b, do.call(rbind, rows_b))

  if (b %% 10L == 0L) {
    elapsed <- (proc.time() - t_start)[["elapsed"]]
    cat(sprintf("  rep %d / %d (%.1f분 경과)\n", b, B, elapsed / 60))
  }
}

results_df <- do.call(rbind, all_rows)
diag_df    <- do.call(rbind, diag_rows)
rownames(results_df) <- NULL

# ── FULL benchmark: b/k와 무관 -> 1회 계산 후 broadcast ─────────────────────
test_mse_full <- mean((y_test - X_test %*% beta_full)^2)
full_rows <- expand.grid(rep = seq_len(B), k = k_grid)
full_rows$method    <- "FULL"
full_rows$fail      <- FALSE
full_rows$pseudo_er <- 0
full_rows$test_mse  <- test_mse_full
full_rows$n_sub     <- N_train
full_rows$clip_rate <- 0
full_rows <- full_rows[, names(results_df)]

results_df <- rbind(results_df, full_rows)

cat(sprintf("\n총 rows: %d\n", nrow(results_df)))
cat(sprintf("fail 비율: %.4f\n",
            mean(results_df$fail[results_df$method != "FULL"], na.rm = TRUE)))
cat(sprintf("RF 평균 fit+predict 시간: %.1f초/rep\n",
            mean(diag_df$t_param + diag_df$t_nonparam)))

# =============================================================================
# 6. 요약 출력
# =============================================================================
sub_df <- results_df[!results_df$fail & results_df$method != "FULL", ]
methods_order <- c("SRS", "LEV-IPW", "OPT-homo",
                    "OPT-hetero-plugin-param", "OPT-hetero-plugin-nonparam")

cat("\n===== Pseudo Excess Risk (mean ± MCSE) =====\n")
for (k in k_grid) {
  cat(sprintf("\n  k = %d\n", k))
  for (met in methods_order) {
    v <- sub_df$pseudo_er[sub_df$method == met & sub_df$k == k]
    cat(sprintf("    %-28s : %.6f ± %.6f\n", met, mean(v), sd(v) / sqrt(length(v))))
  }
}

cat("\n===== Test MSE (mean ± MCSE) =====\n")
for (k in k_grid) {
  cat(sprintf("\n  k = %d\n", k))
  for (met in c("FULL", methods_order)) {
    v <- results_df$test_mse[results_df$method == met & results_df$k == k & !results_df$fail]
    cat(sprintf("    %-28s : %.6f ± %.6f\n", met, mean(v), sd(v) / sqrt(length(v))))
  }
}

# =============================================================================
# 7. 저장
# =============================================================================
if (!dir.exists("../results")) dir.create("../results", recursive = TRUE)
saveRDS(results_df, result_path)
saveRDS(diag_df, sub("real_results", "real_diag", result_path))
cat(sprintf("\n✓ 저장 완료: %s\n", result_path))