# =============================================================================
# 02_run_real.R
# NYC Yellow Taxi — Repeated subsampling experiment
# 실행 위치: real_data/R/
# =============================================================================

source("../sim_final/R/00_utils.R")   # compute_pi, fit_ipw

source("config_real.R")

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

p1 <- ncol(X_train)   # intercept 포함 = 5

cat(sprintf("N_train = %d / N_test = %d / p+1 = %d\n",
            N_train, length(y_test), p1))
cat(sprintf("k_grid : %s\n", paste(k_grid, collapse = ", ")))
cat(sprintf("n0 = %d / B = %d\n\n", n0, B))

# =============================================================================
# 1. 전체 데이터 leverage (한 번만 계산)
# =============================================================================
A_hat_inv <- solve(A_hat)
ell       <- rowSums((X_train %*% A_hat_inv) * X_train)

# =============================================================================
# 2. Plugin sigma 추정
#    |ê_i| ~ X  전체 선형 모형 (DGP 구조 가정 없음)
# =============================================================================
#

estimate_sigma_plugin <- function(X_pilot, y_pilot, X_full) {
  beta_pilot <- tryCatch(
    as.vector(solve(crossprod(X_pilot), crossprod(X_pilot, y_pilot))),
    error = function(e) NULL
  )
  if (is.null(beta_pilot)) return(rep(1, nrow(X_full)))

  abs_resid <- abs(as.vector(y_pilot - X_pilot %*% beta_pilot))

  # trip_duration(열 3)만 분산 driver로 사용
  X_var <- cbind(1, X_pilot[, 3])
  fit_var <- tryCatch(lm.fit(X_var, log(abs_resid + 1e-6)), error = function(e) NULL)
  if (is.null(fit_var)) return(rep(1, nrow(X_full)))

  X_var_full <- cbind(1, X_full[, 3])
  sigma_hat <- exp(as.vector(X_var_full %*% coef(fit_var)))
  pmax(sigma_hat, 1e-6)
}

# =============================================================================
# 3. 단일 replicate
# =============================================================================
run_one_rep_real <- function(k, seed_rep) {
  set.seed(seed_rep)

  # pilot
  idx_pilot <- sample(N_train, n0)
  sigma_hat <- estimate_sigma_plugin(
    X_train[idx_pilot, ], y_train[idx_pilot], X_train
  )

  methods_list <- list(
    list(name = "UNI-IPW",          score = rep(1, N_train)),
    list(name = "LEV-IPW",          score = ell),
    list(name = "OPT-homo",         score = sqrt(ell)),
    list(name = "OPT-hetero-plugin",score = sqrt(sigma_hat * ell))
  )

  out <- vector("list", length(methods_list) + 1L)
  names(out) <- c(sapply(methods_list, `[[`, "name"), "FULL")

  for (m in methods_list) {
    pi  <- compute_pi(m$score, k)
    idx <- which(rbinom(N_train, 1L, pi) == 1L)

    clip_rate <- mean(pi >= 1 - 1e-9)

    if (length(idx) < p1 + 1L) {
      out[[m$name]] <- data.frame(
        method = m$name, k = k, fail = TRUE,
        pseudo_er = NA_real_, test_mse = NA_real_,
        n_sub = length(idx), clip_rate = clip_rate
      )
      next
    }

    beta_hat <- fit_ipw(X_train, y_train, idx, pi)

    if (is.null(beta_hat)) {
      out[[m$name]] <- data.frame(
        method = m$name, k = k, fail = TRUE,
        pseudo_er = NA_real_, test_mse = NA_real_,
        n_sub = length(idx), clip_rate = clip_rate
      )
      next
    }

    diff      <- beta_hat - beta_full
    pseudo_er <- as.numeric(t(diff) %*% A_hat %*% diff)
    test_mse  <- mean((y_test - X_test %*% beta_hat)^2)

    out[[m$name]] <- data.frame(
      method    = m$name, k = k, fail = FALSE,
      pseudo_er = pseudo_er, test_mse = test_mse,
      n_sub     = length(idx), clip_rate = clip_rate
    )
  }

  # FULL benchmark
  test_mse_full <- mean((y_test - X_test %*% beta_full)^2)
  out[["FULL"]] <- data.frame(
    method    = "FULL", k = k, fail = FALSE,
    pseudo_er = 0, test_mse = test_mse_full,
    n_sub     = N_train, clip_rate = 0
  )

  do.call(rbind, out)
}

# =============================================================================
# 4. 메인 루프
# =============================================================================
set.seed(seed)
all_results <- vector("list", length(k_grid))

cat(sprintf("총 %d runs\n\n", length(k_grid) * B))
t_start <- proc.time()

for (i in seq_along(k_grid)) {
  k <- k_grid[i]
  cat(sprintf("k = %d ...", k))

  reps <- vector("list", B)
  for (b in seq_len(B)) {
    reps[[b]] <- run_one_rep_real(k, seed_rep = seed * 10000L + b)
  }
  all_results[[i]] <- do.call(rbind, reps)

  elapsed <- (proc.time() - t_start)["elapsed"]
  cat(sprintf(" 완료 (%.1f초)\n", elapsed))
}

results_df <- do.call(rbind, all_results)
rownames(results_df) <- NULL

cat(sprintf("\n총 rows: %d\n", nrow(results_df)))
cat(sprintf("fail 비율: %.4f\n",
            mean(results_df$fail[results_df$method != "FULL"], na.rm = TRUE)))

# =============================================================================
# 5. 요약 출력
# =============================================================================
sub_df <- results_df[!results_df$fail & results_df$method != "FULL", ]

cat("\n===== Pseudo Excess Risk (mean ± MCSE) =====\n")
for (k in k_grid) {
  cat(sprintf("\n  k = %d\n", k))
  for (met in c("UNI-IPW", "LEV-IPW", "OPT-homo", "OPT-hetero-plugin")) {
    v <- sub_df$pseudo_er[sub_df$method == met & sub_df$k == k]
    cat(sprintf("    %-22s : %.6f ± %.6f\n",
                met, mean(v), sd(v) / sqrt(length(v))))
  }
}

cat("\n===== Test MSE (mean ± MCSE) =====\n")
for (k in k_grid) {
  cat(sprintf("\n  k = %d\n", k))
  for (met in c("FULL", "UNI-IPW", "LEV-IPW", "OPT-homo", "OPT-hetero-plugin")) {
    v <- results_df$test_mse[results_df$method == met & results_df$k == k &
                               !results_df$fail]
    cat(sprintf("    %-22s : %.6f ± %.6f\n",
                met, mean(v), sd(v) / sqrt(length(v))))
  }
}

# =============================================================================
# 6. 저장
# =============================================================================
if (!dir.exists("../results")) dir.create("../results", recursive = TRUE)
saveRDS(results_df, result_path)
cat(sprintf("\n✓ 저장 완료: %s\n", result_path))