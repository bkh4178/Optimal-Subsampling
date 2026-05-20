# run_sim_p9_4_3.R
# Step 4-3: p=9 Mixed Covariate + Variance Direction Grid
#
# Purpose:
#   Check whether oracle/plugin gains persist as dimension increases.
#   Only dimension changes vs p=3; everything else identical.
#
# DGP:
#   X = cbind(1, x1, ..., x9)
#   x1,...,x9 ~ iid N(0,1)
#   beta = rep(1, 10)
#   sigma(x) varies by var_case:
#     x1_x2_x4: |x1 + x2 + x4| + 0.5
#     x1_x2_x7: |x1 + x2 + x7| + 0.5
#     x1_x4_x7: |x1 + x4 + x7| + 0.5
#
# Plug-in variance model:
#   Fit resid^2 ~ abs_s on the pilot sample, where abs_s matches var_case.
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
var_cases <- c("x1_x2_x4", "x1_x2_x7", "x1_x4_x7")

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

var_case_label <- function(var_case) {
  switch(var_case,
         x1_x2_x4 = "sigma(x)=|x1+x2+x4|+0.5",
         x1_x2_x7 = "sigma(x)=|x1+x2+x7|+0.5",
         x1_x4_x7 = "sigma(x)=|x1+x4+x7|+0.5",
         stop("unknown var_case: ", var_case))
}

make_abs_s_from_X <- function(X, var_case) {
  x1 <- X[, 2]
  x2 <- X[, 3]
  x4 <- X[, 5]
  x7 <- X[, 8]

  s_var <- switch(var_case,
                  x1_x2_x4 = x1 + x2 + x4,
                  x1_x2_x7 = x1 + x2 + x7,
                  x1_x4_x7 = x1 + x4 + x7,
                  stop("unknown var_case: ", var_case))
  abs(s_var)
}

# ── DGP ───────────────────────────────────────────────────────────────────────

generate_data <- function(N, beta, seed = NULL, var_case) {
  if (!is.null(seed)) set.seed(seed)
  x1 <- rnorm(N)
  x2 <- rnorm(N)
  x3 <- rnorm(N)

  x4 <- rbinom(N, 1, 0.5) - 0.5
  x5 <- rbinom(N, 1, 0.5) - 0.5
  x6 <- rbinom(N, 1, 0.5) - 0.5

  x7 <- runif(N, -sqrt(3), sqrt(3))
  x8 <- runif(N, -sqrt(3), sqrt(3))
  x9 <- runif(N, -sqrt(3), sqrt(3))

  X <- cbind(1, x1, x2, x3, x4, x5, x6, x7, x8, x9)

  s_var <- switch(var_case,
                  x1_x2_x4 = x1 + x2 + x4,
                  x1_x2_x7 = x1 + x2 + x7,
                  x1_x4_x7 = x1 + x4 + x7,
                  stop("unknown var_case: ", var_case))
  sigma_vec  <- abs(s_var) + 0.5
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
# Fit resid^2 ~ abs_s on the pilot sample.
# This estimates sigma^2(x), not sigma(x) directly.
# The sampling score uses sqrt(sigma2_hat * ell).

estimate_sigma2_abslinear <- function(dat, idx_pilot, beta_pilot, var_case) {
  resid_pilot <- dat$y[idx_pilot] -
    dat$X[idx_pilot, , drop = FALSE] %*% beta_pilot
  resid2_pilot <- as.vector(resid_pilot)^2

  # E[resid^2 | X] = sigma^2(X), so fit resid^2 ~ abs_s for the current var_case.
  abs_s_p <- make_abs_s_from_X(dat$X[idx_pilot, , drop = FALSE], var_case)
  pilot_df <- data.frame(abs_s = abs_s_p)

  var_model <- lm(resid2_pilot ~ abs_s, data = pilot_df)

  abs_s_f <- make_abs_s_from_X(dat$X, var_case)
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

run_one_rep <- function(dat, k, X_test, mu_test, y_test, var_case) {
  N <- nrow(dat$X)
  results <- list()

  # ── Pilot: 분산 추정 전용 (모집단에 되돌려 놓음) ────────────────────────────
  idx_pilot  <- sample(N, n0, replace = FALSE)
  pi_pilot   <- rep(n0 / N, N)
  beta_pilot <- fit_ipw(dat$X, dat$y, idx_pilot, pi = pi_pilot)
  if (is.null(beta_pilot)) beta_pilot <- beta

  sigma2_hat <- estimate_sigma2_abslinear(dat, idx_pilot, beta_pilot, var_case)

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
cat(" p=9 Mixed Covariate + Variance Direction Grid simulation (Step 4-3)\n")
cat(sprintf(" N=%d / N_test=%d / n0=%d / n_rep=%d\n", N, N_test, n0, n_rep))
cat(sprintf(" k: %s\n", paste(ks, collapse = ", ")))
cat(" X = cbind(1, x1, ..., x9)\n")
cat(" x1,x2,x3 ~ N(0,1), x4,x5,x6 ~ Bernoulli(0.5)-0.5, x7,x8,x9 ~ Uniform(-sqrt(3), sqrt(3))\n")
cat(" variance cases: x1_x2_x4, x1_x2_x7, x1_x4_x7\n")
cat(" variance model: resid^2 ~ abs_s(var_case); score = sqrt(sigma2_hat * ell)\n")
cat(" sampling: n0 then k (pilot returned, final Bernoulli sampling from full N)\n")
cat("====================================================\n")

all_rows <- list()
idx_row  <- 1L

for (var_case in var_cases) {
  cat(sprintf("\nvariance case: %s, %s\n", var_case, var_case_label(var_case)))

  dat      <- generate_data(N, beta, seed = 123, var_case = var_case)
  dat_test <- generate_data(N_test, beta, seed = 999, var_case = var_case)
  X_test   <- dat_test$X
  mu_test  <- dat_test$mu
  y_test   <- dat_test$y

  for (k in ks) {
    cat(sprintf("\n[k = %d]\n", k))
    time_k <- proc.time()

    for (r in seq_len(n_rep)) {
      res <- run_one_rep(dat, k, X_test, mu_test, y_test, var_case)

      for (m in names(res)) {
        rv <- res[[m]]
        all_rows[[idx_row]] <- data.frame(
          var_case       = var_case,
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
}

results_df <- do.call(rbind, all_rows)
saveRDS(results_df, "sim/results_p9_mixed_vargrid.rds")
cat(sprintf("\n저장 완료: sim/results_p9_mixed_vargrid.rds (%d행)\n", nrow(results_df)))

# ── 출력 ──────────────────────────────────────────────────────────────────────

method_order <- c("SRS", "OPT-homo", "OPT-hetero-oracle", "OPT-hetero-plugin")

print_metric <- function(df, var, title) {
  cat(sprintf("\n====================================================\n"))
  cat(sprintf(" %s\n", title))
  cat(sprintf("====================================================\n"))
  agg <- aggregate(as.formula(paste(var, "~ var_case + method + k")),
                   data = df[!df$fail, ], FUN = mean, na.rm = TRUE)

  for (vc in var_cases) {
    cat(sprintf("\nvariance case: %s, %s\n", vc, var_case_label(vc)))
    cat(sprintf("%-25s", "method \\ k"))
    for (k in ks) cat(sprintf("%12d", k))
    cat("\n")
    cat(strrep("-", 25 + 12 * length(ks)), "\n")

    rank_tbl <- lapply(ks, function(k_) {
      sub <- agg[agg$var_case == vc & agg$k == k_, ]
      sub <- sub[order(sub[[var]]), ]
      setNames(seq_len(nrow(sub)), sub$method)
    })

    for (m in method_order) {
      cat(sprintf("%-25s", m))
      for (i in seq_along(ks)) {
        k_  <- ks[i]
        val <- agg[[var]][agg$var_case == vc & agg$method == m & agg$k == k_]
        rk  <- rank_tbl[[i]][m]
        if (length(val) == 0) cat(sprintf("%12s", "NA"))
        else cat(sprintf(" %9.5f(%d)", val, rk))
      }
      cat("\n")
    }
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
cat(sprintf("%-12s %-25s %6s %10s %10s\n", "var_case", "method \\ k", "k", "n_clip", "n_realized"))
cat(strrep("-", 70), "\n")
agg_clip <- aggregate(cbind(n_clip, n_realized) ~ var_case + method + k,
                      data = results_df[!results_df$fail, ], FUN = mean)
for (vc in var_cases) for (m in method_order) for (k_ in ks) {
  sub <- agg_clip[agg_clip$var_case == vc & agg_clip$method == m & agg_clip$k == k_, ]
  if (nrow(sub) == 0) next
  cat(sprintf("%-12s %-25s %6d %10.2f %10.2f\n",
              vc, m, k_, sub$n_clip, sub$n_realized))
}

# ── Diagnostic: IPW instability ──────────────────────────────────────────────

cat("\n====================================================\n")
cat(" [6] IPW instability diagnostic (ESS_pi, max_weight, cond_XtWX)\n")
cat("====================================================\n")
cat(sprintf("%-12s %-25s %6s %10s %12s %14s\n",
            "var_case", "method \\ k", "k", "ESS_pi", "max_weight", "cond_XtWX"))
cat(strrep("-", 90), "\n")
agg_inst <- aggregate(cbind(ess_pi, max_weight, cond_XtWX) ~ var_case + method + k,
                      data = results_df[!results_df$fail, ], FUN = mean, na.rm = TRUE)
for (vc in var_cases) for (m in method_order) for (k_ in ks) {
  sub <- agg_inst[agg_inst$var_case == vc & agg_inst$method == m & agg_inst$k == k_, ]
  if (nrow(sub) == 0) next
  cat(sprintf("%-12s %-25s %6d %10.1f %12.2f %14.2f\n",
              vc, m, k_, sub$ess_pi, sub$max_weight, sub$cond_XtWX))
}

# ── Diagnostic: sigma2_hat 품질 ───────────────────────────────────────────────

cat("\n====================================================\n")
cat(" [7] sigma2_hat 추정 품질 (plugin)\n")
cat("====================================================\n")
cat(sprintf("%-12s %-25s %6s %14s %14s\n",
            "var_case", "method", "k", "sigma_cor", "score_cor"))
cat(strrep("-", 80), "\n")
agg_diag <- aggregate(cbind(diag_sigma_cor, diag_score_cor) ~ var_case + method + k,
                      data = results_df[!results_df$fail &
                                        results_df$method == "OPT-hetero-plugin", ],
                      FUN = mean, na.rm = TRUE)
for (vc in var_cases) for (k_ in ks) {
  sub <- agg_diag[agg_diag$var_case == vc & agg_diag$k == k_, ]
  if (nrow(sub) == 0) next
  cat(sprintf("%-12s %-25s %6d %14.4f %14.4f\n",
              vc, "OPT-hetero-plugin", k_, sub$diag_sigma_cor, sub$diag_score_cor))
}

# ── Ratio 비교 ────────────────────────────────────────────────────────────────

cat("\n====================================================\n")
cat(" [8] homo 대비 비교 (ratio of means)\n")
cat("====================================================\n")
cat(sprintf("%-12s %-10s %-25s %12s %12s %10s %10s\n",
            "var_case", "k", "method", "excess_diff", "excess_ratio", "c1_diff", "c1_ratio"))
cat(strrep("-", 100), "\n")

valid <- results_df[!results_df$fail, ]
for (vc in var_cases) for (k_ in ks) {
  sub_k   <- valid[valid$var_case == vc & valid$k == k_, ]
  homo_ex <- mean(sub_k$excess_risk[sub_k$method == "OPT-homo"], na.rm = TRUE)
  homo_c1 <- mean(sub_k$c1_pi[sub_k$method == "OPT-homo"], na.rm = TRUE)
  for (comp in c("OPT-hetero-oracle", "OPT-hetero-plugin")) {
    comp_ex <- mean(sub_k$excess_risk[sub_k$method == comp], na.rm = TRUE)
    comp_c1 <- mean(sub_k$c1_pi[sub_k$method == comp], na.rm = TRUE)
    cat(sprintf("%-12s %-10d %-25s %12.6f %12.4f %10.4f %10.4f\n",
                vc, k_, comp, comp_ex - homo_ex, comp_ex / homo_ex,
                comp_c1 - homo_c1, comp_c1 / homo_c1))
  }
}

# ── 자동 판정 ─────────────────────────────────────────────────────────────────

cat("\n====================================================\n")
cat(" [자동 판정]\n")
cat("====================================================\n")

for (vc in var_cases) for (k_ in ks) {
  sub_k     <- valid[valid$var_case == vc & valid$k == k_, ]
  homo_ex   <- mean(sub_k$excess_risk[sub_k$method == "OPT-homo"], na.rm = TRUE)
  oracle_ex <- mean(sub_k$excess_risk[sub_k$method == "OPT-hetero-oracle"], na.rm = TRUE)
  plugin_ex <- mean(sub_k$excess_risk[sub_k$method == "OPT-hetero-plugin"], na.rm = TRUE)

  homo_c1   <- mean(sub_k$c1_pi[sub_k$method == "OPT-homo"], na.rm = TRUE)
  oracle_c1 <- mean(sub_k$c1_pi[sub_k$method == "OPT-hetero-oracle"], na.rm = TRUE)
  plugin_c1 <- mean(sub_k$c1_pi[sub_k$method == "OPT-hetero-plugin"], na.rm = TRUE)

  cat(sprintf("\nvariance case: %s, %s / k=%d:\n", vc, var_case_label(vc), k_))
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
