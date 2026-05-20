# run_sim_lowdim5.R
# 저차원 시뮬레이션 (p=2, intercept 포함) — v5
#
# 변경사항:
#   - plugin second-stage: pilot 제외 free set 기준으로 compute_pi 재계산
#   - oracle-two-stage 추가 (pilot cost vs estimation error 분리용)
#   - sigma2_hat vs true sigma2 correlation diagnostic 추가
#
# Methods:
#   SRS                  : k개 uniform
#   OPT-homo             : k개, score = sqrt(ell)
#   OPT-hetero-oracle    : k개, score = sqrt(sigma2 * ell)
#   OPT-hetero-oracle-2s : pilot n0 + second-stage k-n0, true sigma2 사용
#   OPT-hetero-plugin    : pilot n0 + second-stage k-n0, sigma2_hat 사용

set.seed(42)
if (!dir.exists("sim")) dir.create("sim")

# ── 설정 ──────────────────────────────────────────────────────────────────────

N      <- 10000
N_test <- 10000
n0     <- 100
n_rep  <- 500
ks     <- c(500, 1000, 2000)
beta   <- c(1, 1, 1)

# ── 유틸리티 ──────────────────────────────────────────────────────────────────

compute_pi <- function(score, k, max_iter = 100) {
  N  <- length(score)
  pi <- score / sum(score) * k
  fixed <- rep(FALSE, N)
  for (iter in seq_len(max_iter)) {
    new_clip <- !fixed & (pi > 1)
    if (!any(new_clip)) break
    fixed[new_clip] <- TRUE
    pi[fixed] <- 1
    budget <- k - sum(fixed)
    free   <- !fixed
    if (budget <= 0 || !any(free)) { pi[free] <- 0; break }
    pi[free] <- score[free] / sum(score[free]) * budget
  }
  pmax(pmin(pi, 1), 0)
}

fit_ipw <- function(X_design, y, idx, pi) {
  w  <- 1 / pi[idx]
  df <- as.data.frame(X_design[idx, , drop = FALSE])
  colnames(df) <- paste0("V", seq_len(ncol(df)))
  df$y <- y[idx]
  tryCatch(
    coef(lm(y ~ . - 1, data = df, weights = w)),
    error = function(e) NULL
  )
}

# ── DGP ───────────────────────────────────────────────────────────────────────

generate_data <- function(N, beta, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  x1 <- rnorm(N)
  x2 <- rnorm(N)
  X  <- cbind(1, x1, x2)
  sigma_vec  <- abs(x1 + x2) + 0.5
  sigma2_vec <- sigma_vec^2
  mu <- as.vector(X %*% beta)
  y  <- mu + rnorm(N, 0, sigma_vec)
  A_hat     <- crossprod(X) / N
  A_hat_inv <- solve(A_hat)
  ell       <- rowSums((X %*% A_hat_inv) * X)
  list(X = X, y = y, mu = mu,
       sigma_vec = sigma_vec, sigma2_vec = sigma2_vec,
       ell = ell, A_hat = A_hat, beta = beta)
}

# ── 단일 replicate 실행 ────────────────────────────────────────────────────────

run_one_rep <- function(dat, k, X_test, mu_test, y_test) {
  N <- nrow(dat$X)
  p <- ncol(dat$X)
  results <- list()

  # ── Pilot (two-stage methods 공용) ─────────────────────────────────────────
  idx_pilot <- sample(N, n0, replace = FALSE)
  free      <- rep(TRUE, N)
  free[idx_pilot] <- FALSE

  pi_pilot   <- rep(n0 / N, N)
  beta_pilot <- fit_ipw(dat$X, dat$y, idx_pilot, pi_pilot)
  if (is.null(beta_pilot)) beta_pilot <- beta

  # sigma2_hat: pilot 잔차 → log-linear → 전체 복원
  resid_pilot <- dat$y[idx_pilot] - dat$X[idx_pilot, , drop = FALSE] %*% beta_pilot
  log_resid2  <- log(as.vector(resid_pilot)^2 + 1e-10)
  pilot_df    <- as.data.frame(dat$X[idx_pilot, , drop = FALSE])
  colnames(pilot_df) <- paste0("V", seq_len(p))
  var_model   <- lm(log_resid2 ~ . - 1, data = pilot_df)
  full_df     <- as.data.frame(dat$X)
  colnames(full_df) <- paste0("V", seq_len(p))
  sigma2_hat  <- exp(predict(var_model, newdata = full_df))
  bad         <- !is.finite(sigma2_hat)
  if (any(bad)) sigma2_hat[bad] <- median(sigma2_hat[!bad])

  # diagnostic: sigma2_hat vs true sigma2
  diag_sigma_cor <- cor(sigma2_hat, dat$sigma2_vec)
  plugin_score   <- sqrt(sigma2_hat * dat$ell)
  oracle_score   <- sqrt(dat$sigma2_vec * dat$ell)
  diag_score_cor <- cor(plugin_score, oracle_score)

  # ── Method 실행 함수 ────────────────────────────────────────────────────────

  # single-stage: pilot 없이 k개 전체
  run_single <- function(method, score) {
    pi_raw    <- compute_pi(score, k)
    idx_final <- which(rbinom(N, 1, pi_raw) == 1)
    n_clip    <- sum(pi_raw >= 1)
    beta_hat  <- fit_ipw(dat$X, dat$y, idx_final, pi_raw)
    if (is.null(beta_hat)) return(list(
      method=method, fail=TRUE, excess_risk=NA, pred_mse=NA,
      param_mse=NA, c1_pi=NA, n_realized=length(idx_final),
      n_clip=n_clip, diag_sigma_cor=NA, diag_score_cor=NA))
    list(
      method      = method, fail = FALSE,
      excess_risk = mean((X_test %*% (beta_hat - beta))^2),
      pred_mse    = mean((y_test - X_test %*% beta_hat)^2),
      param_mse   = sum((beta_hat - beta)^2),
      c1_pi       = mean(dat$sigma2_vec * dat$ell / pmax(pi_raw, 1e-10)),
      n_realized  = length(idx_final), n_clip = n_clip,
      diag_sigma_cor = NA, diag_score_cor = NA
    )
  }

  # two-stage: pilot n0 고정 + free set에서 k-n0개
  run_two_stage <- function(method, score_free) {
    budget <- k - n0

    # free set에서만 compute_pi → budget 손실 없음
    pi_free        <- rep(0, N)
    pi_free[free]  <- compute_pi(score_free[free], budget)

    idx_second <- which(rbinom(N, 1, pi_free) == 1)
    idx_final  <- c(idx_pilot, idx_second)

    # IPW pi: pilot은 n0/N, second-stage는 pi_free
    pi_ipw          <- pi_free
    pi_ipw[idx_pilot] <- n0 / N

    n_clip   <- sum(pi_free[free] >= 1)
    beta_hat <- fit_ipw(dat$X, dat$y, idx_final, pi_ipw)

    if (is.null(beta_hat)) return(list(
      method=method, fail=TRUE, excess_risk=NA, pred_mse=NA,
      param_mse=NA, c1_pi=NA, n_realized=length(idx_final),
      n_clip=n_clip, diag_sigma_cor=diag_sigma_cor, diag_score_cor=diag_score_cor))
    list(
      method      = method, fail = FALSE,
      excess_risk = mean((X_test %*% (beta_hat - beta))^2),
      pred_mse    = mean((y_test - X_test %*% beta_hat)^2),
      param_mse   = sum((beta_hat - beta)^2),
      c1_pi       = mean(dat$sigma2_vec * dat$ell / pmax(pi_ipw, 1e-10)),
      n_realized  = length(idx_final), n_clip = n_clip,
      diag_sigma_cor = diag_sigma_cor, diag_score_cor = diag_score_cor
    )
  }

  # single-stage methods (k개 온전히)
  results[["SRS"]]               <- run_single("SRS",               rep(1, N))
  results[["OPT-homo"]]          <- run_single("OPT-homo",          sqrt(dat$ell))
  results[["OPT-hetero-oracle"]] <- run_single("OPT-hetero-oracle", sqrt(dat$sigma2_vec * dat$ell))

  # two-stage methods (pilot n0 포함)
  results[["OPT-hetero-oracle-2s"]] <- run_two_stage("OPT-hetero-oracle-2s", sqrt(dat$sigma2_vec * dat$ell))
  results[["OPT-hetero-plugin"]]    <- run_two_stage("OPT-hetero-plugin",    sqrt(sigma2_hat * dat$ell))

  results
}

# ── 메인 루프 ─────────────────────────────────────────────────────────────────

cat("====================================================\n")
cat(" 저차원 시뮬레이션 v5 (p=2, intercept 포함)\n")
cat(sprintf(" N=%d / N_test=%d / n0=%d / n_rep=%d\n", N, N_test, n0, n_rep))
cat(sprintf(" k: %s\n", paste(ks, collapse = ", ")))
cat(" sigma(x) = |x1+x2| + 0.5\n")
cat(" single-stage: SRS/homo/oracle (k개)\n")
cat(" two-stage: oracle-2s/plugin (pilot n0 + second k-n0, free set 기준)\n")
cat("====================================================\n")

dat      <- generate_data(N, beta, seed = 123)
dat_test <- generate_data(N_test, beta, seed = 999)
X_test   <- dat_test$X
mu_test  <- dat_test$mu
y_test   <- dat_test$y

all_rows <- list()
idx_row  <- 1L

for (k in ks) {
  cat(sprintf("\n[k = %d]\n", k))
  time_k <- proc.time()

  for (r in seq_len(n_rep)) {
    res <- run_one_rep(dat, k, X_test, mu_test, y_test)

    for (m in names(res)) {
      rv <- res[[m]]
      all_rows[[idx_row]] <- data.frame(
        k              = k, rep = r,
        method         = rv$method,
        fail           = rv$fail,
        excess_risk    = rv$excess_risk,
        pred_mse       = rv$pred_mse,
        param_mse      = rv$param_mse,
        c1_pi          = rv$c1_pi,
        n_realized     = rv$n_realized,
        n_clip         = rv$n_clip,
        diag_sigma_cor = rv$diag_sigma_cor,
        diag_score_cor = rv$diag_score_cor,
        stringsAsFactors = FALSE
      )
      idx_row <- idx_row + 1L
    }

    if (r %% 100 == 0) cat(sprintf("  %d / %d\n", r, n_rep))
  }

  elapsed <- round((proc.time() - time_k)["elapsed"], 1)
  cat(sprintf("  완료 (%.1f초)\n", elapsed))
}

results_df <- do.call(rbind, all_rows)
saveRDS(results_df, "sim/results_lowdim5.rds")
cat(sprintf("\n저장 완료: sim/results_lowdim5.rds (%d행)\n", nrow(results_df)))

# ── 출력 ──────────────────────────────────────────────────────────────────────

method_order <- c("SRS", "OPT-homo", "OPT-hetero-oracle",
                  "OPT-hetero-oracle-2s", "OPT-hetero-plugin")

print_metric <- function(df, var, title) {
  cat(sprintf("\n====================================================\n"))
  cat(sprintf(" %s\n", title))
  cat(sprintf("====================================================\n"))
  agg <- aggregate(as.formula(paste(var, "~ method + k")),
                   data = df[!df$fail, ], FUN = mean, na.rm = TRUE)
  cat(sprintf("%-25s", "method \\ k"))
  for (k in ks) cat(sprintf("%12d", k))
  cat("\n")
  cat(strrep("-", 25 + 12 * length(ks)), "\n")

  rank_tbl <- lapply(ks, function(k_) {
    sub <- agg[agg$k == k_, ]
    sub <- sub[order(sub[[var]]), ]
    setNames(seq_len(nrow(sub)), sub$method)
  })

  for (m in method_order) {
    cat(sprintf("%-25s", m))
    for (i in seq_along(ks)) {
      k_  <- ks[i]
      val <- agg[[var]][agg$method == m & agg$k == k_]
      rk  <- rank_tbl[[i]][m]
      if (length(val) == 0) cat(sprintf("%12s", "NA"))
      else cat(sprintf(" %9.5f(%d)", val, rk))
    }
    cat("\n")
  }
}

print_metric(results_df, "excess_risk", "[1] 평균 excess_risk (Primary)")
print_metric(results_df, "c1_pi",       "[2] 평균 c1_pi (leading term diagnostic)")
print_metric(results_df, "pred_mse",    "[3] 평균 pred_mse (Secondary)")
print_metric(results_df, "param_mse",   "[4] 평균 param_mse")

# ── Clip 및 realized sample size ─────────────────────────────────────────────

cat("\n====================================================\n")
cat(" [5] 평균 clip 개수 및 실현 샘플 크기\n")
cat("====================================================\n")
cat(sprintf("%-25s %6s %10s %10s\n", "method \\ k", "k", "n_clip", "n_realized"))
cat(strrep("-", 55), "\n")
agg_clip <- aggregate(cbind(n_clip, n_realized) ~ method + k,
                      data = results_df[!results_df$fail, ], FUN = mean)
for (m in method_order) for (k_ in ks) {
  sub <- agg_clip[agg_clip$method == m & agg_clip$k == k_, ]
  if (nrow(sub) == 0) next
  cat(sprintf("%-25s %6d %10.2f %10.2f\n", m, k_, sub$n_clip, sub$n_realized))
}

# ── Diagnostic: sigma2_hat 품질 ───────────────────────────────────────────────

cat("\n====================================================\n")
cat(" [6] sigma2_hat 추정 품질 (two-stage methods)\n")
cat("====================================================\n")
cat(sprintf("%-25s %6s %14s %14s\n", "method", "k", "sigma_cor", "score_cor"))
cat(strrep("-", 65), "\n")

two_stage_methods <- c("OPT-hetero-oracle-2s", "OPT-hetero-plugin")
agg_diag <- aggregate(cbind(diag_sigma_cor, diag_score_cor) ~ method + k,
                      data = results_df[!results_df$fail & results_df$method %in% two_stage_methods, ],
                      FUN = mean, na.rm = TRUE)
for (m in two_stage_methods) for (k_ in ks) {
  sub <- agg_diag[agg_diag$method == m & agg_diag$k == k_, ]
  if (nrow(sub) == 0) next
  cat(sprintf("%-25s %6d %14.4f %14.4f\n",
              m, k_, sub$diag_sigma_cor, sub$diag_score_cor))
}

# ── Ratio 비교 (homo 대비) ────────────────────────────────────────────────────

cat("\n====================================================\n")
cat(" [7] homo 대비 비교 (ratio of means)\n")
cat("====================================================\n")
cat(sprintf("%-10s %-25s %12s %12s %10s %10s\n",
            "k", "method", "excess_diff", "excess_ratio", "c1_diff", "c1_ratio"))
cat(strrep("-", 85), "\n")

valid <- results_df[!results_df$fail, ]
for (k_ in ks) {
  sub_k   <- valid[valid$k == k_, ]
  homo_ex <- mean(sub_k$excess_risk[sub_k$method == "OPT-homo"], na.rm = TRUE)
  homo_c1 <- mean(sub_k$c1_pi[sub_k$method == "OPT-homo"], na.rm = TRUE)
  for (comp in c("OPT-hetero-oracle", "OPT-hetero-oracle-2s", "OPT-hetero-plugin")) {
    comp_ex <- mean(sub_k$excess_risk[sub_k$method == comp], na.rm = TRUE)
    comp_c1 <- mean(sub_k$c1_pi[sub_k$method == comp], na.rm = TRUE)
    cat(sprintf("%-10d %-25s %12.6f %12.4f %10.4f %10.4f\n",
                k_, comp, comp_ex - homo_ex, comp_ex / homo_ex,
                comp_c1 - homo_c1, comp_c1 / homo_c1))
  }
}

# ── 자동 판정 ─────────────────────────────────────────────────────────────────

cat("\n====================================================\n")
cat(" [자동 판정]\n")
cat("====================================================\n")

for (k_ in ks) {
  sub_k <- valid[valid$k == k_, ]
  homo_ex    <- mean(sub_k$excess_risk[sub_k$method == "OPT-homo"],           na.rm = TRUE)
  oracle_ex  <- mean(sub_k$excess_risk[sub_k$method == "OPT-hetero-oracle"],  na.rm = TRUE)
  oracle2_ex <- mean(sub_k$excess_risk[sub_k$method == "OPT-hetero-oracle-2s"], na.rm = TRUE)
  plugin_ex  <- mean(sub_k$excess_risk[sub_k$method == "OPT-hetero-plugin"],  na.rm = TRUE)

  cat(sprintf("\nk=%d:\n", k_))
  cat(sprintf("  oracle(single) vs homo: %.4f → %s\n",
              oracle_ex / homo_ex,
              ifelse(oracle_ex < homo_ex, "oracle wins ✓", "oracle loses")))
  cat(sprintf("  oracle-2s vs homo:      %.4f → %s\n",
              oracle2_ex / homo_ex,
              ifelse(oracle2_ex < homo_ex, "oracle-2s wins ✓", "pilot cost가 손해")))
  cat(sprintf("  plugin vs oracle-2s:    %.4f → %s\n",
              plugin_ex / oracle2_ex,
              ifelse(plugin_ex < oracle2_ex * 1.05,
                     "plugin estimation error 작음",
                     "variance model misspecification 의심")))
}

cat("\n")
