# simulation2_addon.R
# ============================================================
# 친구 코드(simulation2.R) 이분산 파트와 동일한 DGP 사용
# 기존 코드 구조 최대한 유지 + 추가 metric 계산
#
# 전제: simulation2.R 이분산 파트가 먼저 실행되어 아래가 존재해야 함:
#   X_pop, full_data_hetero, test_data_hetero,
#   leverage_part, beta_true, N, sigma_vec
#
# 추가 method: OPT-homo / OPT-hetero-oracle / OPT-hetero-plugin
# 추가 metric: excess_risk, param_mse, c1_pi, pi diagnostics,
#              weighted design diagnostics
# ============================================================

if (!dir.exists("sim")) dir.create("sim")

# ── Helper functions ──────────────────────────────────────────────────────────

compute_pi_addon <- function(score, k, max_iter = 100) {
  # iterative clipping: sum(pi) = k, pi <= 1
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

pi_diagnostics <- function(pi) {
  inv_pi <- 1 / pmax(pi, 1e-12)
  list(
    min_pi      = min(pi),
    q01_pi      = unname(quantile(pi, 0.01)),
    median_pi   = median(pi),
    q99_pi      = unname(quantile(pi, 0.99)),
    max_pi      = max(pi),
    q99_inv_pi  = unname(quantile(inv_pi, 0.99)),
    max_inv_pi  = max(inv_pi),
    clip_rate   = mean(pi >= 1 - 1e-12)
  )
}

weighted_design_diag <- function(X_mat, delta, pi) {
  # A_hat = N^{-1} sum_i (delta_i/pi_i) x_i x_i^T
  idx <- which(delta == 1)
  if (length(idx) < 2) return(list(cond_A_hat = NA, min_eig = NA, max_eig = NA))
  tryCatch({
    w     <- delta[idx] / pi[idx]
    XtW   <- t(X_mat[idx, , drop = FALSE]) * w
    A_hat <- XtW %*% X_mat[idx, , drop = FALSE] / nrow(X_mat)
    eigs  <- eigen(A_hat, symmetric = TRUE, only.values = TRUE)$values
    list(cond_A_hat = max(eigs) / max(min(eigs), 1e-12),
         min_eig    = min(eigs),
         max_eig    = max(eigs))
  }, error = function(e) list(cond_A_hat = NA, min_eig = NA, max_eig = NA))
}

compute_c1 <- function(sigma2_true, ell, pi) {
  mean(sigma2_true * ell / pmax(pi, 1e-12))
}

safe_metric_row <- function(method, n_total, rep,
                             beta_hat, pi, delta,
                             X_mat, X_test, y_test, beta_true,
                             sigma2_true, ell) {
  # beta_hat이 NULL이면 fail
  if (is.null(beta_hat)) {
    return(data.frame(
      method = method, n_total = n_total, rep = rep,
      fail = TRUE,
      excess_risk = NA, pred_mse = NA, param_mse = NA, c1_pi = NA,
      n_sub = sum(delta),
      min_pi = NA, q01_pi = NA, median_pi = NA, q99_pi = NA, max_pi = NA,
      q99_inv_pi = NA, max_inv_pi = NA, clip_rate = NA,
      cond_A_hat = NA, min_eig_A_hat = NA, max_eig_A_hat = NA,
      stringsAsFactors = FALSE
    ))
  }

  # metrics
  excess_risk <- mean((X_test %*% (beta_hat - beta_true))^2)
  pred_mse    <- mean((y_test - X_test %*% beta_hat)^2)
  param_mse   <- sum((beta_hat - beta_true)^2)
  c1_pi       <- compute_c1(sigma2_true, ell, pi)

  # diagnostics
  pd  <- pi_diagnostics(pi)
  wdd <- weighted_design_diag(X_mat, delta, pi)

  data.frame(
    method = method, n_total = n_total, rep = rep,
    fail = FALSE,
    excess_risk = excess_risk,
    pred_mse    = pred_mse,
    param_mse   = param_mse,
    c1_pi       = c1_pi,
    n_sub       = sum(delta),
    min_pi      = pd$min_pi,
    q01_pi      = pd$q01_pi,
    median_pi   = pd$median_pi,
    q99_pi      = pd$q99_pi,
    max_pi      = pd$max_pi,
    q99_inv_pi  = pd$q99_inv_pi,
    max_inv_pi  = pd$max_inv_pi,
    clip_rate   = pd$clip_rate,
    cond_A_hat     = wdd$cond_A_hat,
    min_eig_A_hat  = wdd$min_eig,
    max_eig_A_hat  = wdd$max_eig,
    stringsAsFactors = FALSE
  )
}

# ── 공통 설정 ─────────────────────────────────────────────────────────────────

X_mat_addon  <- as.matrix(full_data_hetero[, -1])   # N x (p+1), intercept 포함
y_addon      <- full_data_hetero$y
X_test_addon <- as.matrix(test_data_hetero[, -1])
y_test_addon <- test_data_hetero$y

# true sigma2: sigma_vec = |x1 + x2| + 0.5, sigma2 = sigma_vec^2
sigma2_true_addon <- sigma_vec^2
sigma2_norm_addon <- sigma2_true_addon / mean(sigma2_true_addon)  # 정규화 (oracle용)

ell_addon <- leverage_part   # simulation2.R에서 계산된 leverage score 재사용

n_reps_addon  <- 1000
n_pilot_addon <- 100
ns_addon      <- c(200, 400, 800)
methods_addon <- c("OPT-homo", "OPT-hetero-oracle", "OPT-hetero-plugin")

cat("\n")
cat("============================================================\n")
cat(" [추가 실험] OPT-homo / OPT-hetero-oracle / OPT-hetero-plugin\n")
cat(" DGP: simulation2.R 이분산 (sigma = |x1+x2| + 0.5)\n")
cat(sprintf(" N=%d / N_test=%d / n_reps=%d / n_pilot=%d\n",
            N, nrow(X_test_addon), n_reps_addon, n_pilot_addon))
cat("============================================================\n")

# ── 메인 루프 ─────────────────────────────────────────────────────────────────

all_rows <- list()
idx_row  <- 1L

for (n_total in ns_addon) {
  cat(sprintf("\n[n_total = %d]\n", n_total))

  for (r in seq_len(n_reps_addon)) {

    # ── OPT-homo ──────────────────────────────────────────────────────────────
    pi_homo   <- compute_pi_addon(sqrt(ell_addon), n_total)
    delta_homo <- rbinom(N, 1, pi_homo)
    idx_homo   <- which(delta_homo == 1)

    beta_homo <- tryCatch({
      w <- 1 / pi_homo[idx_homo]
      fit <- lm(y ~ . - 1, data = full_data_hetero[idx_homo, ],
                weights = w)
      coef(fit)
    }, error = function(e) NULL)

    all_rows[[idx_row]] <- safe_metric_row(
      "OPT-homo", n_total, r, beta_homo,
      pi_homo, delta_homo,
      X_mat_addon, X_test_addon, y_test_addon, beta_true,
      sigma2_norm_addon, ell_addon)
    idx_row <- idx_row + 1L

    # ── OPT-hetero-oracle ─────────────────────────────────────────────────────
    pi_oracle   <- compute_pi_addon(sqrt(sigma2_norm_addon * ell_addon), n_total)
    delta_oracle <- rbinom(N, 1, pi_oracle)
    idx_oracle   <- which(delta_oracle == 1)

    beta_oracle <- tryCatch({
      w <- 1 / pi_oracle[idx_oracle]
      fit <- lm(y ~ . - 1, data = full_data_hetero[idx_oracle, ],
                weights = w)
      coef(fit)
    }, error = function(e) NULL)

    all_rows[[idx_row]] <- safe_metric_row(
      "OPT-hetero-oracle", n_total, r, beta_oracle,
      pi_oracle, delta_oracle,
      X_mat_addon, X_test_addon, y_test_addon, beta_true,
      sigma2_norm_addon, ell_addon)
    idx_row <- idx_row + 1L

    # ── OPT-hetero-plugin ─────────────────────────────────────────────────────
    # Step 1: uniform pilot (친구 코드 방식 유지)
    idx_pilot <- sample(N, n_pilot_addon, replace = FALSE)
    d_pilot   <- full_data_hetero[idx_pilot, ]
    fit_pilot <- tryCatch(lm(y ~ . - 1, data = d_pilot), error = function(e) NULL)

    pi_plugin   <- rep(n_total / N, N)  # fallback
    beta_plugin <- NULL

    if (!is.null(fit_pilot)) {
      # Step 2: sigma2_hat (친구 코드 방식: 전체 N에 predict → raw squared residual)
      #y_hat_pilot  <- predict(fit_pilot, full_data_hetero)
      #sigma_sq_hat <- (full_data_hetero$y - y_hat_pilot)^2

      y_hat_pilot  <- predict(fit_pilot, d_pilot)            # pilot에만 predict
      sigma_sq_hat_pilot <- (d_pilot$y - y_hat_pilot)^2     # pilot y만 사용

      # 전체 N으로 확장은 predict로 (y 없이 X만 사용)
      sigma_sq_hat <- predict(
        lm(sigma_sq_hat_pilot ~ . - 1, data = d_pilot),     # pilot에서 모델 적합
        newdata = full_data_hetero                            # 전체 X로 predict
      )
      sigma_sq_hat <- pmax(sigma_sq_hat, 0)                  # 음수 방지
      
      # Step 3: pi 계산
      omega_i    <- sigma_sq_hat * ell_addon
      pi_plugin  <- compute_pi_addon(sqrt(omega_i), n_total)

      delta_plugin <- rbinom(N, 1, pi_plugin)
      idx_plugin   <- which(delta_plugin == 1)

      beta_plugin <- tryCatch({
        w <- 1 / pi_plugin[idx_plugin]
        fit <- lm(y ~ . - 1, data = full_data_hetero[idx_plugin, ],
                  weights = w)
        coef(fit)
      }, error = function(e) NULL)
    } else {
      delta_plugin <- rbinom(N, 1, pi_plugin)
    }

    all_rows[[idx_row]] <- safe_metric_row(
      "OPT-hetero-plugin", n_total, r, beta_plugin,
      pi_plugin, delta_plugin,
      X_mat_addon, X_test_addon, y_test_addon, beta_true,
      sigma2_norm_addon, ell_addon)
    idx_row <- idx_row + 1L

    if (r %% 200 == 0) cat(sprintf("  %d / %d 완료\n", r, n_reps_addon))
  }
}

results_df <- do.call(rbind, all_rows)
saveRDS(results_df, "sim/simulation2_addon_metrics.rds")
cat(sprintf("\nraw 결과 저장: sim/simulation2_addon_metrics.rds (%d행)\n",
            nrow(results_df)))

# ── 요약 계산 ─────────────────────────────────────────────────────────────────

agg_mean <- function(df, var) {
  aggregate(as.formula(paste(var, "~ method + n_total")),
            data = df[!df$fail, ], FUN = mean, na.rm = TRUE)
}

agg_sd <- function(df, var) {
  aggregate(as.formula(paste(var, "~ method + n_total")),
            data = df[!df$fail, ], FUN = sd, na.rm = TRUE)
}

summary_list <- list(
  excess_risk = agg_mean(results_df, "excess_risk"),
  pred_mse    = agg_mean(results_df, "pred_mse"),
  param_mse   = agg_mean(results_df, "param_mse"),
  c1_pi       = agg_mean(results_df, "c1_pi"),
  fail_rate   = aggregate(fail ~ method + n_total, data = results_df, FUN = mean)
)
saveRDS(summary_list, "sim/simulation2_addon_summary.rds")
cat("summary 저장: sim/simulation2_addon_summary.rds\n")

# ── 출력표 ────────────────────────────────────────────────────────────────────

print_table <- function(title, var, df, higher_is_worse = TRUE) {
  cat(sprintf("\n====================================================\n"))
  cat(sprintf(" %s\n", title))
  cat(sprintf("====================================================\n"))
  agg <- agg_mean(df, var)
  cat(sprintf("%-25s", paste0("method \\ n_total")))
  for (n in ns_addon) cat(sprintf("%10d", n))
  cat("\n")
  cat(strrep("-", 25 + 10 * length(ns_addon)), "\n")

  for (m in methods_addon) {
    cat(sprintf("%-25s", m))
    for (n in ns_addon) {
      val <- agg[[var]][agg$method == m & agg$n_total == n]
      if (length(val) == 0) cat(sprintf("%10s", "NA"))
      else cat(sprintf("%10.5f", val))
    }
    cat("\n")
  }

  # 순위
  cat("\n  [순위]\n")
  for (n in ns_addon) {
    vals <- sapply(methods_addon, function(m) {
      v <- agg[[var]][agg$method == m & agg$n_total == n]
      if (length(v) == 0) NA else v
    })
    ord <- order(vals, na.last = TRUE)
    if (!higher_is_worse) ord <- rev(ord)
    cat(sprintf("  n=%d: %s\n", n,
                paste(methods_addon[ord], collapse = " > ")))
  }
}

print_table("[출력표 1] 평균 pred_mse (y_test 기준, 친구 코드 비교용)",
            "pred_mse", results_df)
print_table("[출력표 2] 평균 excess_risk (true mean 기준, 이론 지표)",
            "excess_risk", results_df)
print_table("[출력표 3] 평균 param_mse",
            "param_mse", results_df)
print_table("[출력표 4] 평균 c1_pi (낮을수록 leading objective 좋음)",
            "c1_pi", results_df)

# 출력표 5: instability diagnostic
cat("\n====================================================\n")
cat(" [출력표 5] IPW/Sampling instability diagnostic\n")
cat("====================================================\n")
diag_vars <- c("n_sub", "q99_inv_pi", "max_inv_pi", "clip_rate",
                "cond_A_hat", "min_eig_A_hat")

for (n in ns_addon) {
  cat(sprintf("\n--- n_total = %d ---\n", n))
  cat(sprintf("%-25s %8s %10s %10s %8s %10s %10s %7s\n",
              "method", "n_sub", "q99_inv_pi", "max_inv_pi",
              "clip%", "cond_A", "min_eig", "fail%"))
  cat(strrep("-", 90), "\n")

  for (m in methods_addon) {
    sub <- results_df[results_df$method == m & results_df$n_total == n, ]
    valid <- !sub$fail
    fail_r <- mean(sub$fail)
    cat(sprintf("%-25s %8.1f %10.2f %10.2f %8.4f %10.2f %10.6f %7.4f\n",
                m,
                mean(sub$n_sub[valid], na.rm = TRUE),
                mean(sub$q99_inv_pi[valid], na.rm = TRUE),
                mean(sub$max_inv_pi[valid], na.rm = TRUE),
                mean(sub$clip_rate[valid], na.rm = TRUE),
                mean(sub$cond_A_hat[valid], na.rm = TRUE),
                mean(sub$min_eig_A_hat[valid], na.rm = TRUE),
                fail_r))
  }
}

# 출력표 6: homo 대비 비교 (ratio of means)
cat("\n====================================================\n")
cat(" [출력표 6] homo 대비 oracle/plugin 비교 (ratio of means)\n")
cat("====================================================\n")
cat(sprintf("%-10s %-20s %12s %12s %12s %12s\n",
            "n_total", "comparison",
            "excess_diff", "excess_ratio", "c1_diff", "c1_ratio"))
cat(strrep("-", 82), "\n")

for (n in ns_addon) {
  sub <- results_df[results_df$n_total == n & !results_df$fail, ]

  homo_ex <- mean(sub$excess_risk[sub$method == "OPT-homo"], na.rm = TRUE)
  homo_c1 <- mean(sub$c1_pi[sub$method == "OPT-homo"], na.rm = TRUE)

  for (comp in c("OPT-hetero-oracle", "OPT-hetero-plugin")) {
    comp_ex <- mean(sub$excess_risk[sub$method == comp], na.rm = TRUE)
    comp_c1 <- mean(sub$c1_pi[sub$method == comp], na.rm = TRUE)
    cat(sprintf("%-10d %-20s %12.6f %12.4f %12.4f %12.4f\n",
                n, comp,
                comp_ex - homo_ex,
                comp_ex / homo_ex,
                comp_c1 - homo_c1,
                comp_c1 / homo_c1))
  }
}

# ── 자동 판정 ─────────────────────────────────────────────────────────────────

cat("\n====================================================\n")
cat(" [자동 판정]\n")
cat("====================================================\n")

for (n in ns_addon) {
  sub <- results_df[results_df$n_total == n & !results_df$fail, ]
  homo_ex   <- mean(sub$excess_risk[sub$method == "OPT-homo"], na.rm = TRUE)
  homo_c1   <- mean(sub$c1_pi[sub$method == "OPT-homo"], na.rm = TRUE)
  oracle_ex <- mean(sub$excess_risk[sub$method == "OPT-hetero-oracle"], na.rm = TRUE)
  oracle_c1 <- mean(sub$c1_pi[sub$method == "OPT-hetero-oracle"], na.rm = TRUE)

  c1_ratio  <- oracle_c1 / homo_c1
  ex_ratio  <- oracle_ex / homo_ex

  if (c1_ratio < 1 & ex_ratio < 1) {
    verdict <- "oracle improves both leading objective and finite-sample excess risk"
  } else if (c1_ratio < 1 & ex_ratio >= 1) {
    verdict <- "oracle improves leading objective but worsens finite-sample excess risk; check IPW instability"
  } else {
    verdict <- "possible issue in pi/ell/objective calculation"
  }
  cat(sprintf("n=%d: c1_ratio=%.4f excess_ratio=%.4f → %s\n",
              n, c1_ratio, ex_ratio, verdict))
}

cat("\n====================================================\n")
cat(" 추가 실험 완료\n")
cat("====================================================\n")
