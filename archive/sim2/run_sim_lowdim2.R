# run_sim_lowdim.R
# 저차원 시뮬레이션 (p=2, intercept 포함)
#
# 확정 설계:
#   - Design matrix: cbind(1, x1, x2), x1/x2 ~ N(0,1)
#   - beta = (1, 1, 1)
#   - sigma(x) = |x1 + x2| + 0.5  (directional heteroscedasticity)
#   - N = 10000, N_test = 10000, n0 = 100, n_rep = 500
#   - k = 500, 1000, 2000
#   - Pilot n0개 고정 포함 + 나머지 k-n0개 Bernoulli sampling
#   - Clipping: iterative budget 재분배 방식 (sum(pi)=k 보장) + clip 개수 출력
#
# Methods: SRS, OPT-homo, OPT-hetero-oracle, OPT-hetero-plugin
# Metrics: excess_risk, c1_pi, pred_mse, param_mse

set.seed(42)
if (!dir.exists("sim")) dir.create("sim")

# ── 설정 ──────────────────────────────────────────────────────────────────────

N      <- 10000
N_test <- 10000
n0     <- 100
n_rep  <- 500
ks     <- c(500, 1000, 2000)
beta   <- c(1, 1, 1)   # intercept + x1 + x2

# ── 유틸리티 ──────────────────────────────────────────────────────────────────

# score → pi 변환 (iterative clipping + budget 재분배)
# sum(pi) = k, pi <= 1 동시 만족
compute_pi_simple <- function(score, k, max_iter = 100) {
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

# IPW estimator (lm with weights)
fit_ipw <- function(X_design, y, idx, pi) {
  w <- 1 / pi[idx]
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
  X  <- cbind(1, x1, x2)   # N x 3, intercept 포함

  # 분산 구조: sigma(x) = |x1 + x2| + 0.5
  sigma_vec  <- abs(x1 + x2) + 0.5
  sigma2_vec <- sigma_vec^2

  mu <- as.vector(X %*% beta)
  y  <- mu + rnorm(N, 0, sigma_vec)

  # A_hat, leverage (intercept 포함 기준)
  A_hat     <- crossprod(X) / N
  A_hat_inv <- solve(A_hat)
  ell       <- rowSums((X %*% A_hat_inv) * X)

  list(X = X, y = y, mu = mu,
       sigma_vec = sigma_vec, sigma2_vec = sigma2_vec,
       ell = ell, A_hat = A_hat, beta = beta)
}

# ── 단일 replicate 실행 ────────────────────────────────────────────────────────

run_one_rep <- function(dat, k, X_test, mu_test, y_test) {
  N   <- nrow(dat$X)
  p   <- ncol(dat$X)

  results <- list()

  # ── Pilot sampling ──────────────────────────────────────────────────────────
  idx_pilot  <- sample(N, n0, replace = FALSE)
  pi_pilot   <- rep(n0 / N, N)

  beta_pilot <- fit_ipw(dat$X, dat$y, idx_pilot, pi_pilot)
  if (is.null(beta_pilot)) beta_pilot <- beta  # fallback

  # sigma2_hat: pilot 잔차 → log-linear 모델 → 전체 복원
  resid_pilot  <- dat$y[idx_pilot] - dat$X[idx_pilot, , drop=FALSE] %*% beta_pilot
  log_resid2   <- log(as.vector(resid_pilot)^2 + 1e-10)
  pilot_df     <- as.data.frame(dat$X[idx_pilot, , drop=FALSE])
  colnames(pilot_df) <- paste0("V", seq_len(p))
  var_model    <- lm(log_resid2 ~ . - 1, data = pilot_df)
  full_df      <- as.data.frame(dat$X)
  colnames(full_df) <- paste0("V", seq_len(p))
  sigma2_hat   <- exp(predict(var_model, newdata = full_df))
  bad          <- !is.finite(sigma2_hat)
  if (any(bad)) sigma2_hat[bad] <- median(sigma2_hat[!bad])

  # ── Score 계산 ──────────────────────────────────────────────────────────────
  scores <- list(
    "SRS"               = rep(1, N),
    "OPT-homo"          = sqrt(dat$ell),
    "OPT-hetero-oracle" = sqrt(dat$sigma2_vec * dat$ell),
    "OPT-hetero-plugin" = sqrt(sigma2_hat * dat$ell)
  )

  for (method in names(scores)) {

    # pi 계산 (k-n0 budget 기준)
    budget <- k - n0
    pi_raw <- compute_pi_simple(scores[[method]], budget)

    # pilot 인덱스 제외하고 나머지에서 Bernoulli sampling
    pi_full <- pi_raw
    pi_full[idx_pilot] <- 0   # pilot은 이미 포함, 중복 방지

    idx_second <- which(rbinom(N, 1, pi_full) == 1)
    idx_final  <- union(idx_pilot, idx_second)

    # final pi 벡터 (IPW용): pilot은 n0/N, 나머지는 pi_raw
    pi_ipw <- pi_raw
    pi_ipw[idx_pilot] <- n0 / N

    # clip 개수
    n_clip <- sum(pi_raw >= 1)

    # beta 추정
    beta_hat <- fit_ipw(dat$X, dat$y, idx_final, pi_ipw)

    if (is.null(beta_hat)) {
      results[[method]] <- list(
        method      = method,
        fail        = TRUE,
        excess_risk = NA, pred_mse = NA,
        param_mse   = NA, c1_pi    = NA,
        n_realized  = length(idx_final),
        n_clip      = n_clip
      )
      next
    }

    # metrics
    excess_risk <- mean((X_test %*% (beta_hat - beta))^2)
    pred_mse    <- mean((y_test - X_test %*% beta_hat)^2)
    param_mse   <- sum((beta_hat - beta)^2)
    c1_pi       <- mean(dat$sigma2_vec * dat$ell / pmax(pi_ipw, 1e-10))

    results[[method]] <- list(
      method      = method,
      fail        = FALSE,
      excess_risk = excess_risk,
      pred_mse    = pred_mse,
      param_mse   = param_mse,
      c1_pi       = c1_pi,
      n_realized  = length(idx_final),
      n_clip      = n_clip
    )
  }

  results
}

# ── 메인 루프 ─────────────────────────────────────────────────────────────────

cat("====================================================\n")
cat(" 저차원 시뮬레이션 (p=2, intercept 포함)\n")
cat(sprintf(" N=%d / N_test=%d / n0=%d / n_rep=%d\n", N, N_test, n0, n_rep))
cat(sprintf(" k: %s\n", paste(ks, collapse=", ")))
cat(sprintf(" sigma(x) = |x1+x2| + 0.5\n"))
cat("====================================================\n")

# test data 고정
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
    dat <- generate_data(N, beta)
    res <- run_one_rep(dat, k, X_test, mu_test, y_test)

    for (m in names(res)) {
      rv <- res[[m]]
      all_rows[[idx_row]] <- data.frame(
        k           = k,
        rep         = r,
        method      = rv$method,
        fail        = rv$fail,
        excess_risk = rv$excess_risk,
        pred_mse    = rv$pred_mse,
        param_mse   = rv$param_mse,
        c1_pi       = rv$c1_pi,
        n_realized  = rv$n_realized,
        n_clip      = rv$n_clip,
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
saveRDS(results_df, "sim/results_lowdim.rds")
cat(sprintf("\n저장 완료: sim/results_lowdim.rds (%d행)\n", nrow(results_df)))

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
      k_ <- ks[i]
      val <- agg[[var]][agg$method == m & agg$k == k_]
      rk  <- rank_tbl[[i]][m]
      if (length(val) == 0) cat(sprintf("%12s", "NA"))
      else cat(sprintf(" %9.5f(%d)", val, rk))
    }
    cat("\n")
  }
}

print_metric(results_df, "excess_risk", "[1] 평균 excess_risk (낮을수록 좋음, Primary)")
print_metric(results_df, "c1_pi",       "[2] 평균 c1_pi (leading term diagnostic)")
print_metric(results_df, "pred_mse",    "[3] 평균 pred_mse (y_test 기준, Secondary)")
print_metric(results_df, "param_mse",   "[4] 평균 param_mse (beta 추정 오차)")

# ── Clip 개수 출력 ────────────────────────────────────────────────────────────

cat("\n====================================================\n")
cat(" [5] 평균 clip 개수 (n_clip) 및 실현 샘플 크기\n")
cat("====================================================\n")
cat(sprintf("%-25s %6s %10s %10s\n", "method \\ k", "k", "n_clip", "n_realized"))
cat(strrep("-", 55), "\n")

agg_clip <- aggregate(cbind(n_clip, n_realized) ~ method + k,
                      data = results_df[!results_df$fail, ], FUN = mean)

for (m in method_order) {
  for (k_ in ks) {
    sub <- agg_clip[agg_clip$method == m & agg_clip$k == k_, ]
    if (nrow(sub) == 0) next
    cat(sprintf("%-25s %6d %10.2f %10.2f\n",
                m, k_, sub$n_clip, sub$n_realized))
  }
}

# ── Ratio 비교 (homo 대비, ratio of means) ────────────────────────────────────

cat("\n====================================================\n")
cat(" [6] homo 대비 oracle/plugin 비교 (ratio of means)\n")
cat("====================================================\n")
cat(sprintf("%-10s %-22s %12s %12s %10s %10s\n",
            "k", "method", "excess_diff", "excess_ratio", "c1_diff", "c1_ratio"))
cat(strrep("-", 80), "\n")

valid <- results_df[!results_df$fail, ]
for (k_ in ks) {
  sub_k <- valid[valid$k == k_, ]
  homo_ex <- mean(sub_k$excess_risk[sub_k$method == "OPT-homo"], na.rm=TRUE)
  homo_c1 <- mean(sub_k$c1_pi[sub_k$method == "OPT-homo"], na.rm=TRUE)

  for (comp in c("OPT-hetero-oracle", "OPT-hetero-plugin")) {
    comp_ex <- mean(sub_k$excess_risk[sub_k$method == comp], na.rm=TRUE)
    comp_c1 <- mean(sub_k$c1_pi[sub_k$method == comp], na.rm=TRUE)
    cat(sprintf("%-10d %-22s %12.6f %12.4f %10.4f %10.4f\n",
                k_, comp,
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

for (k_ in ks) {
  sub_k    <- valid[valid$k == k_, ]
  homo_ex  <- mean(sub_k$excess_risk[sub_k$method == "OPT-homo"], na.rm=TRUE)
  homo_c1  <- mean(sub_k$c1_pi[sub_k$method == "OPT-homo"], na.rm=TRUE)
  oracle_ex <- mean(sub_k$excess_risk[sub_k$method == "OPT-hetero-oracle"], na.rm=TRUE)
  oracle_c1 <- mean(sub_k$c1_pi[sub_k$method == "OPT-hetero-oracle"], na.rm=TRUE)

  c1_ratio <- oracle_c1 / homo_c1
  ex_ratio <- oracle_ex / homo_ex

  if (c1_ratio < 1 && ex_ratio < 1) {
    verdict <- "oracle improves both c1_pi and excess_risk ✓"
  } else if (c1_ratio < 1 && ex_ratio >= 1) {
    verdict <- "oracle improves c1_pi but worsens excess_risk → IPW instability 의심"
  } else {
    verdict <- "c1_pi도 개선 안 됨 → pi/ell 계산 문제 의심"
  }
  cat(sprintf("k=%4d: c1_ratio=%.4f  excess_ratio=%.4f  → %s\n",
              k_, c1_ratio, ex_ratio, verdict))
}

cat("\n")
