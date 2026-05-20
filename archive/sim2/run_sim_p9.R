# run_sim_p9.R
# Step 4-1: p=9 Normal-only, aligned low-dimensional variance
#
# Purpose:
#   Check whether oracle/plugin gains persist as dimension increases.
#   Only dimension changes vs p=3; everything else identical.
#
# DGP:
#   X = cbind(1, x1, ..., x9)
#   x1,...,x9 ~ iid N(0,1)
#   beta = rep(1, 10)
#   sigma(x) = |x1 + x2| + 0.5
#
# Plug-in variance model:
#   Fit resid^2 ~ |x1 + x2| on the pilot sample.
#   sigma^2(x) estimated via linear model; score = sqrt(sigma2_hat * ell).
#
# Additional diagnostics: ESS_pi, max_weight, cond_XtWX

set.seed(42)
if (!dir.exists("sim")) dir.create("sim")

# ── 설정 ──────────────────────────────────────────────────────────────────────

N      <- 10000
N_test <- 10000
n0     <- 100
n_rep  <- 500
ks     <- c(1000)          # k 고정
beta   <- rep(1, 10)       # intercept + x1,...,x9

# ── 유틸리티 ──────────────────────────────────────────────────────────────────

compute_pi <- function(score, k, max_iter = 100) {
  N <- length(score)
  if (k <= 0) return(rep(0, N))
  if (all(!is.finite(score)) || sum(score, na.rm = TRUE) <= 0) {
    return(rep(min(k / N, 1), N))
  }
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
  x1 <- rnorm(N); x2 <- rnorm(N); x3 <- rnorm(N)
  x4 <- rnorm(N); x5 <- rnorm(N); x6 <- rnorm(N)
  x7 <- rnorm(N); x8 <- rnorm(N); x9 <- rnorm(N)
  X  <- cbind(1, x1, x2, x3, x4, x5, x6, x7, x8, x9)
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

# ── Plug-in variance model: estimate sigma^2(x) from resid^2 ────────────────
# Fit resid^2 ~ |x1 + x2| on the pilot sample.
# This estimates sigma^2(x), not sigma(x) directly.
# The sampling score uses sqrt(sigma2_hat * ell).

estimate_sigma2_abslinear <- function(dat, idx_pilot, beta_pilot) {
  x1_p <- dat$X[idx_pilot, 2]
  x2_p <- dat$X[idx_pilot, 3]

  resid_pilot <- dat$y[idx_pilot] -
    dat$X[idx_pilot, , drop = FALSE] %*% beta_pilot
  resid2_pilot <- as.vector(resid_pilot)^2

  # E[resid^2 | X] = sigma^2(X), so fit resid^2 ~ |x1+x2|.
  # x3 is included in the mean model but intentionally excluded from the variance model.

  abs_s_p <- abs(x1_p + x2_p)
  pilot_df <- data.frame(abs_s = abs_s_p)

  var_model <- lm(resid2_pilot ~ abs_s, data = pilot_df)

  x1_f <- dat$X[, 2]
  x2_f <- dat$X[, 3]
  abs_s_f <- abs(x1_f + x2_f)
  full_df <- data.frame(abs_s = abs_s_f)

  sigma2_hat <- predict(var_model, newdata = full_df)

  # 음수 방지
  bad <- !is.finite(sigma2_hat) | sigma2_hat <= 0
  if (any(bad)) {
    good_vals <- sigma2_hat[!bad]
    sigma2_hat[bad] <- if (length(good_vals) > 0) median(good_vals)
                       else median(dat$sigma2_vec)
  }

  pmax(sigma2_hat, 1e-10)
}

# ── 단일 replicate 실행 ────────────────────────────────────────────────────────

run_one_rep <- function(dat, k, X_test, mu_test, y_test) {
  N <- nrow(dat$X)
  results <- list()

  # ── Pilot: 분산 추정 전용 (모집단에 되돌려 놓음) ────────────────────────────
  idx_pilot  <- sample(N, n0, replace = FALSE)
  pi_pilot   <- rep(n0 / N, N)
  beta_pilot <- fit_ipw(dat$X, dat$y, idx_pilot, pi = pi_pilot)
  if (is.null(beta_pilot)) beta_pilot <- beta

  sigma2_hat <- estimate_sigma2_abslinear(dat, idx_pilot, beta_pilot)

  plugin_score <- sqrt(sigma2_hat * dat$ell)
  oracle_score <- sqrt(dat$sigma2_vec * dat$ell)

  diag_sigma_cor <- safe_cor(sigma2_hat, dat$sigma2_vec)
  diag_score_cor <- safe_cor(plugin_score, oracle_score)

  # ── Method 실행 함수 ────────────────────────────────────────────────────────
  # n0 then k 방식: 모든 method가 전체 N에서 k개 Bernoulli sampling
  # two-stage methods도 pilot을 모집단에 되돌려 놓고 k개 뽑음
  # → weight = 1/pi_i 통일, stage-scaled weight 불필요

  run_method <- function(method, score, is_plugin = FALSE) {
    pi_raw    <- compute_pi(score, k)
    idx_final <- which(rbinom(N, 1, pi_raw) == 1)
    n_clip    <- sum(pi_raw >= 1)

    # ESS_pi: effective sample size of sampling probabilities
    ess_pi     <- sum(pi_raw)^2 / sum(pi_raw^2)
    # max_weight: maximum IPW weight
    max_weight <- max(1 / pmax(pi_raw[pi_raw > 0], 1e-10))

    beta_hat  <- fit_ipw(dat$X, dat$y, idx_final, pi = pi_raw)

    # cond_XtWX: condition number of weighted design matrix
    cond_XtWX <- tryCatch({
      if (length(idx_final) < ncol(dat$X)) return(NA_real_)
      w   <- 1 / pmax(pi_raw[idx_final], 1e-10)
      XtW <- t(dat$X[idx_final, , drop=FALSE]) * w
      XtWX <- XtW %*% dat$X[idx_final, , drop=FALSE]
      ev <- eigen(XtWX, symmetric=TRUE, only.values=TRUE)$values
      max(ev) / max(min(ev), 1e-12)
    }, error = function(e) NA_real_)

    if (is.null(beta_hat)) {
      return(list(
        method = method, fail = TRUE,
        excess_risk = NA, pred_mse = NA, param_mse = NA,
        c1_pi = NA, n_realized = length(idx_final), n_clip = n_clip,
        ess_pi = ess_pi, max_weight = max_weight, cond_XtWX = cond_XtWX,
        diag_sigma_cor = ifelse(is_plugin, diag_sigma_cor, NA),
        diag_score_cor = ifelse(is_plugin, diag_score_cor, NA)
      ))
    }

    list(
      method      = method, fail = FALSE,
      excess_risk = mean((X_test %*% (beta_hat - beta))^2),
      pred_mse    = mean((y_test - X_test %*% beta_hat)^2),
      param_mse   = sum((beta_hat - beta)^2),
      c1_pi       = mean(dat$sigma2_vec * dat$ell / pmax(pi_raw, 1e-10)),
      n_realized  = length(idx_final),
      n_clip      = n_clip,
      ess_pi      = ess_pi,
      max_weight  = max_weight,
      cond_XtWX   = cond_XtWX,
      diag_sigma_cor = ifelse(is_plugin, diag_sigma_cor, NA),
      diag_score_cor = ifelse(is_plugin, diag_score_cor, NA)
    )
  }

  results[["SRS"]]               <- run_method("SRS",               rep(1, N))
  results[["OPT-homo"]]          <- run_method("OPT-homo",          sqrt(dat$ell))
  results[["OPT-hetero-oracle"]] <- run_method("OPT-hetero-oracle", sqrt(dat$sigma2_vec * dat$ell))
  results[["OPT-hetero-plugin"]] <- run_method("OPT-hetero-plugin", sqrt(sigma2_hat * dat$ell), is_plugin = TRUE)

  results
}

# ── 메인 루프 ─────────────────────────────────────────────────────────────────

cat("====================================================\n")
cat(" p=9 Normal-only simulation (Step 4-1)\n")
cat(sprintf(" N=%d / N_test=%d / n0=%d / n_rep=%d\n", N, N_test, n0, n_rep))
cat(sprintf(" k: %s\n", paste(ks, collapse = ", ")))
cat(" X = cbind(1, x1, x2, x3)\n")
cat(" x1, x2 ~ N(0,1); x3 ~ Bernoulli(0.5) - 0.5\n")
cat(" sigma(x) = |x1+x2| + 0.5\n")
cat(" variance model: resid^2 ~ |x1+x2|; score = sqrt(sigma2_hat * ell)\n")
cat(" sampling: n0 then k (pilot returned, final Bernoulli sampling from full N)\n")
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
        ess_pi         = rv$ess_pi,
        max_weight     = rv$max_weight,
        cond_XtWX      = rv$cond_XtWX,
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
saveRDS(results_df, "sim/results_p9.rds")
cat(sprintf("\n저장 완료: sim/results_p9.rds (%d행)\n", nrow(results_df)))

# ── 출력 ──────────────────────────────────────────────────────────────────────

method_order <- c("SRS", "OPT-homo", "OPT-hetero-oracle", "OPT-hetero-plugin")

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
print_metric(results_df, "c1_pi",       "[2] 평균 c1_pi (diagnostic)")
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

# ── Diagnostic: IPW instability ──────────────────────────────────────────────

cat("\n====================================================\n")
cat(" [6] IPW instability diagnostic (ESS_pi, max_weight, cond_XtWX)\n")
cat("====================================================\n")
cat(sprintf("%-25s %6s %10s %12s %14s\n", "method \ k", "k", "ESS_pi", "max_weight", "cond_XtWX"))
cat(strrep("-", 75), "\n")
agg_inst <- aggregate(cbind(ess_pi, max_weight, cond_XtWX) ~ method + k,
                      data = results_df[!results_df$fail, ], FUN = mean, na.rm = TRUE)
for (m in method_order) for (k_ in ks) {
  sub <- agg_inst[agg_inst$method == m & agg_inst$k == k_, ]
  if (nrow(sub) == 0) next
  cat(sprintf("%-25s %6d %10.1f %12.2f %14.2f\n",
              m, k_, sub$ess_pi, sub$max_weight, sub$cond_XtWX))
}

# ── Diagnostic: sigma2_hat 품질 ───────────────────────────────────────────────

cat("\n====================================================\n")
cat(" [7] sigma2_hat 추정 품질 (plugin)\n")
cat("====================================================\n")
cat(sprintf("%-25s %6s %14s %14s\n", "method", "k", "sigma_cor", "score_cor"))
cat(strrep("-", 65), "\n")
agg_diag <- aggregate(cbind(diag_sigma_cor, diag_score_cor) ~ method + k,
                      data = results_df[!results_df$fail &
                                        results_df$method == "OPT-hetero-plugin", ],
                      FUN = mean, na.rm = TRUE)
for (k_ in ks) {
  sub <- agg_diag[agg_diag$k == k_, ]
  if (nrow(sub) == 0) next
  cat(sprintf("%-25s %6d %14.4f %14.4f\n",
              "OPT-hetero-plugin", k_, sub$diag_sigma_cor, sub$diag_score_cor))
}

# ── Ratio 비교 ────────────────────────────────────────────────────────────────

cat("\n====================================================\n")
cat(" [8] homo 대비 비교 (ratio of means)\n")
cat("====================================================\n")
cat(sprintf("%-10s %-25s %12s %12s %10s %10s\n",
            "k", "method", "excess_diff", "excess_ratio", "c1_diff", "c1_ratio"))
cat(strrep("-", 85), "\n")

valid <- results_df[!results_df$fail, ]
for (k_ in ks) {
  sub_k   <- valid[valid$k == k_, ]
  homo_ex <- mean(sub_k$excess_risk[sub_k$method == "OPT-homo"], na.rm = TRUE)
  homo_c1 <- mean(sub_k$c1_pi[sub_k$method == "OPT-homo"], na.rm = TRUE)
  for (comp in c("OPT-hetero-oracle", "OPT-hetero-plugin")) {
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
  sub_k     <- valid[valid$k == k_, ]
  homo_ex   <- mean(sub_k$excess_risk[sub_k$method == "OPT-homo"], na.rm = TRUE)
  oracle_ex <- mean(sub_k$excess_risk[sub_k$method == "OPT-hetero-oracle"], na.rm = TRUE)
  plugin_ex <- mean(sub_k$excess_risk[sub_k$method == "OPT-hetero-plugin"], na.rm = TRUE)

  homo_c1   <- mean(sub_k$c1_pi[sub_k$method == "OPT-homo"], na.rm = TRUE)
  oracle_c1 <- mean(sub_k$c1_pi[sub_k$method == "OPT-hetero-oracle"], na.rm = TRUE)
  plugin_c1 <- mean(sub_k$c1_pi[sub_k$method == "OPT-hetero-plugin"], na.rm = TRUE)

  cat(sprintf("\nk=%d:\n", k_))
  cat(sprintf("  oracle / homo(ex):   %.4f → %s\n",
              oracle_ex / homo_ex,
              ifelse(oracle_ex < homo_ex, "oracle gain ✓", "oracle no gain")))
  cat(sprintf("  plugin / homo(ex):   %.4f → %s\n",
              plugin_ex / homo_ex,
              ifelse(plugin_ex < homo_ex, "plugin gain ✓", "plugin no gain")))
  cat(sprintf("  plugin / oracle(ex): %.4f → %s\n",
              plugin_ex / oracle_ex,
              ifelse(plugin_ex < oracle_ex * 1.10,
                     "expected plug-in cost ✓",
                     "large plug-in gap")))

  cat(sprintf("  oracle / homo(c1):   %.4f\n", oracle_c1 / homo_c1))
  cat(sprintf("  plugin / homo(c1):   %.4f\n", plugin_c1 / homo_c1))
  cat(sprintf("  plugin / oracle(c1): %.4f\n", plugin_c1 / oracle_c1))
}