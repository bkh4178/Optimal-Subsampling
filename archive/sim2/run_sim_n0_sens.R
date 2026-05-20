# run_sim_n0_sens.R
# n0 sensitivity 실험 (p=2, 이분산)
#
# k=1000 고정, n0 in {50, 100, 200} 변화
# 기대 메시지:
#   - n0가 커질수록 plugin이 oracle에 가까워지는가
#   - 작은 n0에서도 plugin이 homo보다 이기는가

set.seed(42)
if (!dir.exists("sim")) dir.create("sim")

# ── 설정 ──────────────────────────────────────────────────────────────────────

N      <- 10000
N_test <- 10000
k      <- 1000          # 고정
n0s    <- c(50, 100, 200)
n_rep  <- 500
beta   <- c(1, 1, 1)

# ── 유틸리티 ──────────────────────────────────────────────────────────────────

compute_pi <- function(score, k, max_iter = 100) {
  N <- length(score)
  if (k <= 0) return(rep(0, N))
  if (all(!is.finite(score)) || sum(score, na.rm = TRUE) <= 0)
    return(rep(min(k / N, 1), N))
  score[!is.finite(score)] <- 0
  score <- pmax(score, 0)
  if (sum(score) <= 0) return(rep(min(k / N, 1), N))
  pi    <- score / sum(score) * k
  fixed <- rep(FALSE, N)
  for (iter in seq_len(max_iter)) {
    new_clip <- !fixed & (pi > 1)
    if (!any(new_clip)) break
    fixed[new_clip] <- TRUE
    pi[fixed] <- 1
    budget <- k - sum(fixed)
    free   <- !fixed
    if (budget <= 0 || !any(free)) { pi[free] <- 0; break }
    if (sum(score[free]) <= 0) pi[free] <- budget / sum(free)
    else pi[free] <- score[free] / sum(score[free]) * budget
  }
  pmax(pmin(pi, 1), 0)
}

fit_ipw <- function(X_design, y, idx, pi) {
  if (length(idx) == 0) return(NULL)
  w <- 1 / pmax(pi[idx], 1e-10)
  if (any(!is.finite(w)) || any(w <= 0)) return(NULL)
  df <- as.data.frame(X_design[idx, , drop = FALSE])
  colnames(df) <- paste0("V", seq_len(ncol(df)))
  df$y <- y[idx]
  tryCatch(coef(lm(y ~ . - 1, data = df, weights = w)), error = function(e) NULL)
}

safe_cor <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < 3 || sd(a[ok]) == 0 || sd(b[ok]) == 0) return(NA_real_)
  cor(a[ok], b[ok])
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
       ell = ell, beta = beta)
}

# ── Variance model: resid^2 ~ |x1+x2| ────────────────────
# E[resid^2] = sigma^2(x) 이용
# pilot residual 제곱을 response로, abs(x1+x2)를 predictor로 선형 적합
# → sigma2_hat = predict(var_model), 음수 방지 처리 포함

estimate_sigma2 <- function(dat, idx_pilot, beta_pilot) {
  x1_p <- dat$X[idx_pilot, 2]
  x2_p <- dat$X[idx_pilot, 3]
  resid_pilot  <- dat$y[idx_pilot] -
    dat$X[idx_pilot, , drop = FALSE] %*% beta_pilot
  resid2_pilot <- as.vector(resid_pilot)^2
  abs_s_p  <- abs(x1_p + x2_p)
  var_model <- lm(resid2_pilot ~ abs_s_p)
  x1_f <- dat$X[, 2]; x2_f <- dat$X[, 3]
  abs_s_f  <- abs(x1_f + x2_f)
  sigma2_hat <- predict(var_model, newdata = data.frame(abs_s_p = abs_s_f))
  bad <- !is.finite(sigma2_hat) | sigma2_hat <= 0
  if (any(bad)) sigma2_hat[bad] <- median(sigma2_hat[!bad])
  pmax(sigma2_hat, 1e-10)
}

# ── 단일 replicate ─────────────────────────────────────────────────────────────

run_one_rep <- function(dat, k, n0, X_test, mu_test, y_test) {
  N <- nrow(dat$X)

  # pilot
  idx_pilot  <- sample(N, n0, replace = FALSE)
  pi_pilot   <- rep(n0 / N, N)
  beta_pilot <- fit_ipw(dat$X, dat$y, idx_pilot, pi = pi_pilot)
  if (is.null(beta_pilot)) beta_pilot <- beta

  sigma2_hat <- estimate_sigma2(dat, idx_pilot, beta_pilot)
  plugin_score <- sqrt(sigma2_hat * dat$ell)
  oracle_score <- sqrt(dat$sigma2_vec * dat$ell)
  diag_sigma_cor <- safe_cor(sigma2_hat, dat$sigma2_vec)
  diag_score_cor <- safe_cor(plugin_score, oracle_score)

  run_method <- function(method, score, is_plugin = FALSE) {
    pi_raw    <- compute_pi(score, k)
    idx_final <- which(rbinom(N, 1, pi_raw) == 1)
    n_clip    <- sum(pi_raw >= 1)
    beta_hat  <- fit_ipw(dat$X, dat$y, idx_final, pi = pi_raw)
    if (is.null(beta_hat)) return(list(
      method=method, fail=TRUE, n0=n0,
      excess_risk=NA, pred_mse=NA, param_mse=NA, c1_pi=NA,
      n_realized=length(idx_final), n_clip=n_clip,
      diag_sigma_cor=ifelse(is_plugin, diag_sigma_cor, NA),
      diag_score_cor=ifelse(is_plugin, diag_score_cor, NA)))
    list(
      method      = method, fail = FALSE, n0 = n0,
      excess_risk = mean((X_test %*% (beta_hat - beta))^2),
      pred_mse    = mean((y_test - X_test %*% beta_hat)^2),
      param_mse   = sum((beta_hat - beta)^2),
      c1_pi       = mean(dat$sigma2_vec * dat$ell / pmax(pi_raw, 1e-10)),
      n_realized  = length(idx_final), n_clip = n_clip,
      diag_sigma_cor = ifelse(is_plugin, diag_sigma_cor, NA),
      diag_score_cor = ifelse(is_plugin, diag_score_cor, NA)
    )
  }

  list(
    SRS                  = run_method("SRS",               rep(1, N)),
    OPT_homo             = run_method("OPT-homo",          sqrt(dat$ell)),
    OPT_hetero_oracle    = run_method("OPT-hetero-oracle", sqrt(dat$sigma2_vec * dat$ell)),
    OPT_hetero_plugin    = run_method("OPT-hetero-plugin", sqrt(sigma2_hat * dat$ell), is_plugin = TRUE)
  )
}

# ── 메인 루프 ─────────────────────────────────────────────────────────────────

cat("====================================================\n")
cat(" n0 Sensitivity 실험 (p=2, 이분산)\n")
cat(sprintf(" N=%d / N_test=%d / k=%d / n_rep=%d\n", N, N_test, k, n_rep))
cat(sprintf(" n0: %s\n", paste(n0s, collapse = ", ")))
cat(" sigma(x) = |x1+x2| + 0.5\n")
cat("====================================================\n")

dat      <- generate_data(N, beta, seed = 123)
dat_test <- generate_data(N_test, beta, seed = 999)
X_test   <- dat_test$X
y_test   <- dat_test$y

all_rows <- list()
idx_row  <- 1L

for (n0 in n0s) {
  cat(sprintf("\n[n0 = %d]\n", n0))
  time_n0 <- proc.time()

  for (r in seq_len(n_rep)) {
    res <- run_one_rep(dat, k, n0, X_test, dat_test$mu, y_test)
    for (m in names(res)) {
      rv <- res[[m]]
      all_rows[[idx_row]] <- data.frame(
        n0=rv$n0, rep=r, method=rv$method, fail=rv$fail,
        excess_risk=rv$excess_risk, pred_mse=rv$pred_mse,
        param_mse=rv$param_mse, c1_pi=rv$c1_pi,
        n_realized=rv$n_realized, n_clip=rv$n_clip,
        diag_sigma_cor=rv$diag_sigma_cor,
        diag_score_cor=rv$diag_score_cor,
        stringsAsFactors = FALSE)
      idx_row <- idx_row + 1L
    }
    if (r %% 100 == 0) cat(sprintf("  %d / %d\n", r, n_rep))
  }

  elapsed <- round((proc.time() - time_n0)["elapsed"], 1)
  cat(sprintf("  완료 (%.1f초)\n", elapsed))
}

results_df <- do.call(rbind, all_rows)
saveRDS(results_df, "sim/results_n0_sens.rds")
cat(sprintf("\n저장 완료: sim/results_n0_sens.rds (%d행)\n", nrow(results_df)))

# ── 출력 ──────────────────────────────────────────────────────────────────────

method_order <- c("SRS", "OPT-homo", "OPT-hetero-oracle", "OPT-hetero-plugin")
valid <- results_df[!results_df$fail, ]

# [1] excess_risk
cat("\n====================================================\n")
cat(" [1] 평균 excess_risk (k=1000 고정, n0 변화)\n")
cat("====================================================\n")
cat(sprintf("%-25s", "method \\ n0"))
for (n0 in n0s) cat(sprintf("%10d", n0)); cat("\n")
cat(strrep("-", 25 + 10*length(n0s)), "\n")
agg_ex <- aggregate(excess_risk ~ method + n0, data=valid, FUN=mean)
for (m in method_order) {
  cat(sprintf("%-25s", m))
  for (n0 in n0s) {
    val <- agg_ex$excess_risk[agg_ex$method==m & agg_ex$n0==n0]
    if (length(val)==0) cat(sprintf("%10s","NA")) else cat(sprintf("%10.5f", val))
  }; cat("\n")
}

# [2] c1_pi
cat("\n====================================================\n")
cat(" [2] 평균 c1_pi\n")
cat("====================================================\n")
cat(sprintf("%-25s", "method \\ n0"))
for (n0 in n0s) cat(sprintf("%10d", n0)); cat("\n")
cat(strrep("-", 25 + 10*length(n0s)), "\n")
agg_c1 <- aggregate(c1_pi ~ method + n0, data=valid, FUN=mean)
for (m in method_order) {
  cat(sprintf("%-25s", m))
  for (n0 in n0s) {
    val <- agg_c1$c1_pi[agg_c1$method==m & agg_c1$n0==n0]
    if (length(val)==0) cat(sprintf("%10s","NA")) else cat(sprintf("%10.3f", val))
  }; cat("\n")
}

# [3] plugin ratio (homo 대비, oracle 대비)
cat("\n====================================================\n")
cat(" [3] plugin ratio (ratio of means)\n")
cat("====================================================\n")
cat(sprintf("%-12s %14s %14s %14s %14s\n",
            "n0", "plugin/homo(ex)", "plugin/oracle(ex)", "plugin/homo(c1)", "plugin/oracle(c1)"))
cat(strrep("-", 72), "\n")
for (n0 in n0s) {
  sub    <- valid[valid$n0==n0, ]
  homo_ex   <- mean(sub$excess_risk[sub$method=="OPT-homo"],          na.rm=TRUE)
  oracle_ex <- mean(sub$excess_risk[sub$method=="OPT-hetero-oracle"], na.rm=TRUE)
  plugin_ex <- mean(sub$excess_risk[sub$method=="OPT-hetero-plugin"], na.rm=TRUE)
  homo_c1   <- mean(sub$c1_pi[sub$method=="OPT-homo"],                na.rm=TRUE)
  oracle_c1 <- mean(sub$c1_pi[sub$method=="OPT-hetero-oracle"],       na.rm=TRUE)
  plugin_c1 <- mean(sub$c1_pi[sub$method=="OPT-hetero-plugin"],       na.rm=TRUE)
  cat(sprintf("%-12d %14.4f %14.4f %14.4f %14.4f\n",
              n0,
              plugin_ex/homo_ex, plugin_ex/oracle_ex,
              plugin_c1/homo_c1, plugin_c1/oracle_c1))
}

# [4] sigma2_hat 추정 품질
cat("\n====================================================\n")
cat(" [4] sigma2_hat 추정 품질 (plugin)\n")
cat("====================================================\n")
cat(sprintf("%-12s %14s %14s\n", "n0", "sigma_cor", "score_cor"))
cat(strrep("-", 44), "\n")
agg_diag <- aggregate(cbind(diag_sigma_cor, diag_score_cor) ~ n0,
                      data=valid[valid$method=="OPT-hetero-plugin",],
                      FUN=mean, na.rm=TRUE)
for (n0 in n0s) {
  sub <- agg_diag[agg_diag$n0==n0, ]
  if (nrow(sub)==0) next
  cat(sprintf("%-12d %14.4f %14.4f\n", n0, sub$diag_sigma_cor, sub$diag_score_cor))
}

# [5] realized sample size / clip
cat("\n====================================================\n")
cat(" [5] 평균 n_realized, n_clip\n")
cat("====================================================\n")
cat(sprintf("%-25s %6s %10s %10s\n", "method \\ n0", "n0", "n_realized", "n_clip"))
cat(strrep("-", 55), "\n")
agg_clip <- aggregate(cbind(n_realized, n_clip) ~ method+n0, data=valid, FUN=mean)
for (m in method_order) for (n0 in n0s) {
  sub <- agg_clip[agg_clip$method==m & agg_clip$n0==n0, ]
  if (nrow(sub)==0) next
  cat(sprintf("%-25s %6d %10.2f %10.2f\n", m, n0, sub$n_realized, sub$n_clip))
}

# [6] 자동 판정
cat("\n====================================================\n")
cat(" [자동 판정]\n")
cat("====================================================\n")
for (n0 in n0s) {
  sub       <- valid[valid$n0==n0, ]
  homo_ex   <- mean(sub$excess_risk[sub$method=="OPT-homo"],          na.rm=TRUE)
  oracle_ex <- mean(sub$excess_risk[sub$method=="OPT-hetero-oracle"], na.rm=TRUE)
  plugin_ex <- mean(sub$excess_risk[sub$method=="OPT-hetero-plugin"], na.rm=TRUE)
  cat(sprintf("n0=%3d: plugin/homo=%.4f [%s]  plugin/oracle=%.4f [%s]\n",
              n0,
              plugin_ex/homo_ex,
              ifelse(plugin_ex < homo_ex, "plugin wins ✓", "plugin loses"),
              plugin_ex/oracle_ex,
              ifelse(plugin_ex < oracle_ex*1.1, "near oracle ✓", "gap remains")))
}
cat("\n")
