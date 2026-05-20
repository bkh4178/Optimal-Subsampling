# run_lowdim_oracle_sanity.R
# Low-dimensional sanity check: OPT-homo vs OPT-hetero-oracle
#
# 핵심 변경: design matrix에 intercept 포함
# methods: UNI-IPW, OPT-homo, OPT-hetero-oracle (3개만)
# Stage 1: p=1 / Stage 2: p=2
# Case 1: 등분산 / Case 2: sigma2 = 1 + gamma*(X1^2 + ... + Xp^2)

# ── 설정 ──────────────────────────────────────────────────────────────────────

set.seed(42)

N       <- 10000     # full data size (50000으로 쉽게 변경 가능)
N_test  <- 20000     # test set size
n_rep   <- 500       # replicates
ps      <- c(1, 2)   # covariate dimensions (intercept 제외)
gammas  <- c(1, 2, 4)
ks      <- c(200, 500, 1000, 2000, 5000)
cases   <- c(1, 2)

if (!dir.exists("sim")) dir.create("sim")

# ── 유틸리티 함수 ─────────────────────────────────────────────────────────────

compute_pi <- function(score, k, max_iter = 100) {
  N <- length(score)
  if (any(!is.finite(score)) || any(score < 0)) stop("invalid score")
  if (sum(score) == 0) stop("sum(score) == 0")

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
    pi[free] <- score[free] / sum(score[free]) * budget
  }
  pmax(pmin(pi, 1), 0)
}

pi_diagnostics <- function(pi) {
  inv_pi <- 1 / pmax(pi, 1e-10)
  list(
    clip_rate  = mean(pi >= 1 - 1e-10),
    min_pi     = min(pi),
    q01_pi     = quantile(pi, 0.01),
    median_pi  = median(pi),
    q99_pi     = quantile(pi, 0.99),
    max_pi     = max(pi),
    q99_inv_pi = quantile(inv_pi, 0.99),
    max_inv_pi = max(inv_pi)
  )
}

fit_ipw <- function(X_design, y, idx, pi) {
  # IPW estimator: beta_hat = (X'WX)^{-1} X'Wy, W = diag(1/pi)
  X_sub <- X_design[idx, , drop = FALSE]
  y_sub <- y[idx]
  w     <- 1 / pi[idx]
  XtW   <- t(X_sub) * w
  XtWX  <- XtW %*% X_sub
  XtWy  <- XtW %*% y_sub
  tryCatch(as.vector(solve(XtWX, XtWy)), error = function(e) NULL)
}

# ── DGP ───────────────────────────────────────────────────────────────────────

generate_data <- function(N, p, case, gamma, beta_true) {
  # covariate
  X_raw <- matrix(rnorm(N * p), nrow = N, ncol = p)

  # design matrix with intercept
  X_design <- cbind(1, X_raw)   # N x (p+1)

  # A_hat & leverage (intercept 포함)
  A_hat     <- crossprod(X_design) / N
  A_hat_inv <- solve(A_hat)
  ell       <- rowSums((X_design %*% A_hat_inv) * X_design)

  # sigma2
  if (case == 1) {
    sigma2_raw <- rep(1, N)
  } else {
    sigma2_raw <- 1 + gamma * rowSums(X_raw^2)
  }
  sigma2 <- sigma2_raw / mean(sigma2_raw)

  # response
  mu      <- as.vector(X_design %*% beta_true)
  epsilon <- rnorm(N, 0, sqrt(sigma2))
  y       <- mu + epsilon

  list(X_design = X_design, y = y, mu = mu,
       sigma2 = sigma2, ell = ell, A_hat = A_hat,
       beta_true = beta_true, p = p, case = case, gamma = gamma)
}

# ── 단일 replicate 실행 ────────────────────────────────────────────────────────

run_one_rep <- function(dat, k, X_test_design, mu_test) {
  N <- nrow(dat$X_design)

  results <- list()

  methods <- list(
    list(name = "UNI-IPW",          score = rep(1, N)),
    list(name = "OPT-homo",         score = sqrt(dat$ell)),
    list(name = "OPT-hetero-oracle",score = sqrt(dat$sigma2 * dat$ell))
  )

  for (m in methods) {
    pi  <- compute_pi(m$score, k)
    idx <- which(rbinom(N, 1, pi) == 1)

    if (length(idx) < ncol(dat$X_design) + 1) {
      # singular 방지
      results[[m$name]] <- list(
        method = m$name, fail = TRUE,
        beta_hat = NULL, pi = pi, n_sub = length(idx)
      )
      next
    }

    beta_hat <- fit_ipw(dat$X_design, dat$y, idx, pi)

    if (is.null(beta_hat)) {
      results[[m$name]] <- list(
        method = m$name, fail = TRUE,
        beta_hat = NULL, pi = pi, n_sub = length(idx)
      )
      next
    }

    # metrics
    pred_diff  <- X_test_design %*% (beta_hat - dat$beta_true)
    excess_risk <- mean(pred_diff^2)
    param_mse   <- sum((beta_hat - dat$beta_true)^2)
    c1_pi       <- mean(dat$sigma2 * dat$ell / pi)
    diag        <- pi_diagnostics(pi)

    results[[m$name]] <- c(
      list(method = m$name, fail = FALSE,
           beta_hat = beta_hat, pi = pi,
           excess_risk = excess_risk,
           param_mse   = param_mse,
           c1_pi       = c1_pi,
           n_sub       = length(idx)),
      diag
    )
  }
  results
}

# ── 메인 루프 ─────────────────────────────────────────────────────────────────

all_rows <- list()
idx_row  <- 1L
time_start <- proc.time()

cat("====================================================\n")
cat(sprintf(" Low-dim Oracle Sanity Check\n"))
cat(sprintf(" N=%d / N_test=%d / n_rep=%d\n", N, N_test, n_rep))
cat(sprintf(" ps: %s / gammas: %s / ks: %s\n",
            paste(ps, collapse=","),
            paste(gammas, collapse=","),
            paste(ks, collapse=",")))
cat("====================================================\n")

for (p in ps) {

  beta_true <- rep(1, p + 1)   # intercept 포함 (p+1)개

  # test data 생성 (p, gamma, case별로 독립)
  for (gamma in gammas) {
    for (case in cases) {

      cat(sprintf("\n[p=%d / gamma=%g / Case %d]\n", p, gamma, case))

      # test set (고정)
      dat_test   <- generate_data(N_test, p, case, gamma, beta_true)
      X_test_des <- dat_test$X_design
      mu_test    <- dat_test$mu

      for (k in ks) {

        time_k <- proc.time()

        for (rep in seq_len(n_rep)) {
          dat <- generate_data(N, p, case, gamma, beta_true)
          res <- run_one_rep(dat, k, X_test_des, mu_test)

          for (m_name in names(res)) {
            r <- res[[m_name]]
            row <- data.frame(
              p           = p,
              gamma       = gamma,
              case        = case,
              k           = k,
              rep         = rep,
              method      = m_name,
              fail        = r$fail,
              excess_risk = ifelse(r$fail, NA, r$excess_risk),
              param_mse   = ifelse(r$fail, NA, r$param_mse),
              c1_pi       = ifelse(r$fail, NA, r$c1_pi),
              n_sub       = r$n_sub,
              clip_rate   = ifelse(r$fail, NA, r$clip_rate),
              min_pi      = ifelse(r$fail, NA, r$min_pi),
              q01_pi      = ifelse(r$fail, NA, r$q01_pi),
              median_pi   = ifelse(r$fail, NA, r$median_pi),
              q99_pi      = ifelse(r$fail, NA, r$q99_pi),
              max_pi      = ifelse(r$fail, NA, r$max_pi),
              q99_inv_pi  = ifelse(r$fail, NA, r$q99_inv_pi),
              max_inv_pi  = ifelse(r$fail, NA, r$max_inv_pi),
              stringsAsFactors = FALSE
            )
            all_rows[[idx_row]] <- row
            idx_row <- idx_row + 1L
          }
        }

        elapsed_k <- round((proc.time() - time_k)["elapsed"], 1)
        cat(sprintf("  k = %4d 완료 (%4.1f초)\n", k, elapsed_k))
      }
    }
  }
}

results_df <- do.call(rbind, all_rows)
saveRDS(results_df, "sim/results_lowdim_oracle.rds")

elapsed_total <- round((proc.time() - time_start)["elapsed"], 1)
cat(sprintf("\n총 소요: %.1f초 / 총 행 수: %d\n", elapsed_total, nrow(results_df)))

# ── 출력표 함수 ───────────────────────────────────────────────────────────────

method_order <- c("UNI-IPW", "OPT-homo", "OPT-hetero-oracle")

print_table1 <- function(df) {
  cat("\n====================================================\n")
  cat(" [출력표 1] 평균 excess_risk (순위 포함, 낮을수록 좋음)\n")
  cat("====================================================\n")
  agg <- aggregate(excess_risk ~ p + gamma + case + k + method,
                   data = df, FUN = mean, na.rm = TRUE)
  for (p_ in ps) for (gm in gammas) for (cs in cases) {
    cat(sprintf("\n--- p=%d / gamma=%g / Case %d ---\n", p_, gm, cs))
    cat(sprintf("%-25s", "method \\ k"))
    for (k in ks) cat(sprintf("%10d", k)); cat("\n")
    cat(strrep("-", 25 + 10 * length(ks)), "\n")
    # 순위 계산
    rank_tbl <- lapply(ks, function(k) {
      sub <- agg[agg$p==p_ & agg$gamma==gm & agg$case==cs & agg$k==k, ]
      sub <- sub[order(sub$excess_risk), ]
      setNames(seq_len(nrow(sub)), sub$method)
    })
    names(rank_tbl) <- as.character(ks)
    for (m in method_order) {
      cat(sprintf("%-25s", m))
      for (k in ks) {
        val <- agg$excess_risk[agg$p==p_ & agg$gamma==gm & agg$case==cs &
                                 agg$k==k & agg$method==m]
        rk  <- rank_tbl[[as.character(k)]][m]
        if (length(val)==0) cat(sprintf("%10s","NA"))
        else cat(sprintf(" %7.5f(%d)", val, rk))
      }; cat("\n")
    }
  }
}

print_table2 <- function(df) {
  cat("\n====================================================\n")
  cat(" [출력표 2] 평균 c1_pi (낮을수록 이론적으로 좋음)\n")
  cat("====================================================\n")
  agg <- aggregate(c1_pi ~ p + gamma + case + k + method,
                   data = df, FUN = mean, na.rm = TRUE)
  for (p_ in ps) for (gm in gammas) for (cs in cases) {
    cat(sprintf("\n--- p=%d / gamma=%g / Case %d ---\n", p_, gm, cs))
    cat(sprintf("%-25s", "method \\ k"))
    for (k in ks) cat(sprintf("%10d", k)); cat("\n")
    cat(strrep("-", 25 + 10 * length(ks)), "\n")
    for (m in method_order) {
      cat(sprintf("%-25s", m))
      for (k in ks) {
        val <- agg$c1_pi[agg$p==p_ & agg$gamma==gm & agg$case==cs &
                           agg$k==k & agg$method==m]
        if (length(val)==0) cat(sprintf("%10s","NA"))
        else cat(sprintf("%10.3f", val))
      }; cat("\n")
    }
  }
}

print_table3 <- function(df) {
  cat("\n====================================================\n")
  cat(" [출력표 3] Ratio diagnostic (hetero-oracle / homo)\n")
  cat("====================================================\n")
  for (p_ in ps) for (gm in gammas) for (cs in cases) for (k_ in ks) {
    sub <- df[df$p==p_ & df$gamma==gm & df$case==cs & df$k==k_, ]
    homo   <- sub[sub$method=="OPT-homo", ]
    oracle <- sub[sub$method=="OPT-hetero-oracle", ]
    common_reps <- intersect(homo$rep[!is.na(homo$c1_pi)],
                             oracle$rep[!is.na(oracle$c1_pi)])
    if (length(common_reps) < 10) next
    h_c1  <- homo$c1_pi[homo$rep %in% common_reps]
    o_c1  <- oracle$c1_pi[oracle$rep %in% common_reps]
    h_ex  <- homo$excess_risk[homo$rep %in% common_reps]
    o_ex  <- oracle$excess_risk[oracle$rep %in% common_reps]
    c1_ratio <- o_c1 / h_c1
    ex_ratio <- o_ex / h_ex
    cat(sprintf("p=%d gm=%g cs=%d k=%4d | c1_ratio=%.3f  excess_ratio=%.3f  prop_c1_oracle_better=%.2f  prop_excess_oracle_better=%.2f\n",
                p_, gm, cs, k_,
                mean(c1_ratio, na.rm=TRUE),
                mean(ex_ratio, na.rm=TRUE),
                mean(c1_ratio < 1, na.rm=TRUE),
                mean(ex_ratio < 1, na.rm=TRUE)))
  }
}

print_table4 <- function(df) {
  cat("\n====================================================\n")
  cat(" [출력표 4] IPW instability diagnostic\n")
  cat("====================================================\n")
  for (p_ in ps) for (gm in gammas) for (cs in cases) {
    cat(sprintf("\n--- p=%d / gamma=%g / Case %d ---\n", p_, gm, cs))
    cat(sprintf("%-25s %5s %8s %10s %10s %8s %7s\n",
                "method\\k", "k", "n_sub", "q99_inv_pi", "max_inv_pi", "clip%", "fail%"))
    for (k_ in ks) for (m in method_order) {
      sub <- df[df$p==p_ & df$gamma==gm & df$case==cs & df$k==k_ & df$method==m, ]
      if (nrow(sub)==0) next
      cat(sprintf("%-25s %5d %8.1f %10.2f %10.2f %8.3f %7.3f\n",
                  m, k_,
                  mean(sub$n_sub, na.rm=TRUE),
                  mean(sub$q99_inv_pi, na.rm=TRUE),
                  mean(sub$max_inv_pi, na.rm=TRUE),
                  mean(sub$clip_rate, na.rm=TRUE),
                  mean(sub$fail, na.rm=TRUE)))
    }
  }
}

print_table5 <- function(df) {
  cat("\n====================================================\n")
  cat(" [출력표 5] Case 1 sanity check\n")
  cat("====================================================\n")
  cs <- 1
  for (p_ in ps) for (gm in gammas) for (k_ in ks) {
    sub    <- df[df$p==p_ & df$gamma==gm & df$case==cs & df$k==k_, ]
    homo   <- sub[sub$method=="OPT-homo", ]
    oracle <- sub[sub$method=="OPT-hetero-oracle", ]
    common_reps <- intersect(homo$rep[!is.na(homo$excess_risk)],
                             oracle$rep[!is.na(oracle$excess_risk)])
    if (length(common_reps) < 10) next
    h_ex <- homo$excess_risk[homo$rep %in% common_reps]
    o_ex <- oracle$excess_risk[oracle$rep %in% common_reps]
    diff_ex  <- mean(o_ex - h_ex, na.rm=TRUE)
    ratio_ex <- mean(o_ex / h_ex, na.rm=TRUE)
    cat(sprintf("p=%d gm=%g k=%4d | diff_excess=%.6f  ratio_excess=%.4f\n",
                p_, gm, k_, diff_ex, ratio_ex))
  }
}

print_table6 <- function(df) {
  cat("\n====================================================\n")
  cat(" [출력표 6] Case 2 success check\n")
  cat("====================================================\n")
  cs <- 2
  for (p_ in ps) for (gm in gammas) for (k_ in ks) {
    sub    <- df[df$p==p_ & df$gamma==gm & df$case==cs & df$k==k_, ]
    homo   <- sub[sub$method=="OPT-homo", ]
    oracle <- sub[sub$method=="OPT-hetero-oracle", ]
    common_reps <- intersect(homo$rep[!is.na(homo$excess_risk)],
                             oracle$rep[!is.na(oracle$excess_risk)])
    if (length(common_reps) < 10) next
    h_c1 <- homo$c1_pi[homo$rep %in% common_reps]
    o_c1 <- oracle$c1_pi[oracle$rep %in% common_reps]
    h_ex <- homo$excess_risk[homo$rep %in% common_reps]
    o_ex <- oracle$excess_risk[oracle$rep %in% common_reps]
    cat(sprintf("p=%d gm=%g k=%4d | diff_excess=%.6f ratio_excess=%.4f diff_c1=%.3f ratio_c1=%.4f\n",
                p_, gm, k_,
                mean(o_ex - h_ex, na.rm=TRUE),
                mean(o_ex / h_ex, na.rm=TRUE),
                mean(o_c1 - h_c1, na.rm=TRUE),
                mean(o_c1 / h_c1, na.rm=TRUE)))
  }
}

# ── 출력 실행 ─────────────────────────────────────────────────────────────────

print_table1(results_df)
print_table2(results_df)
print_table3(results_df)
print_table4(results_df)
print_table5(results_df)
print_table6(results_df)

# ── 자동 판정 ─────────────────────────────────────────────────────────────────

cat("\n====================================================\n")
cat(" [자동 판정]\n")
cat("====================================================\n")

for (p_ in ps) for (gm in gammas) {

  # Case 1 sanity check
  sub_c1 <- results_df[results_df$p==p_ & results_df$gamma==gm & results_df$case==1, ]
  homo_c1   <- sub_c1[sub_c1$method=="OPT-homo", ]
  oracle_c1 <- sub_c1[sub_c1$method=="OPT-hetero-oracle", ]
  common <- intersect(homo_c1$rep[!is.na(homo_c1$excess_risk)],
                      oracle_c1$rep[!is.na(oracle_c1$excess_risk)])
  if (length(common) > 10) {
    ratio <- mean(oracle_c1$excess_risk[oracle_c1$rep %in% common] /
                  homo_c1$excess_risk[homo_c1$rep %in% common], na.rm=TRUE)
    verdict <- if (ratio > 0.95 & ratio < 1.05) "PASSED" else "FAILED"
    cat(sprintf("Case 1 sanity check [p=%d, gamma=%g]: ratio=%.4f → %s\n",
                p_, gm, ratio, verdict))
  }

  # Case 2 check
  sub_c2 <- results_df[results_df$p==p_ & results_df$gamma==gm & results_df$case==2, ]
  for (k_ in ks) {
    homo_k   <- sub_c2[sub_c2$method=="OPT-homo"          & sub_c2$k==k_, ]
    oracle_k <- sub_c2[sub_c2$method=="OPT-hetero-oracle" & sub_c2$k==k_, ]
    common2  <- intersect(homo_k$rep[!is.na(homo_k$excess_risk)],
                          oracle_k$rep[!is.na(oracle_k$excess_risk)])
    if (length(common2) < 10) next
    c1_ratio <- mean(oracle_k$c1_pi[oracle_k$rep %in% common2] /
                     homo_k$c1_pi[homo_k$rep %in% common2], na.rm=TRUE)
    ex_ratio <- mean(oracle_k$excess_risk[oracle_k$rep %in% common2] /
                     homo_k$excess_risk[homo_k$rep %in% common2], na.rm=TRUE)
    if (c1_ratio < 1 & ex_ratio < 1) {
      msg <- "theory and finite-sample both support hetero-oracle"
    } else if (c1_ratio < 1 & ex_ratio >= 1) {
      msg <- "leading objective improves but finite-sample excess risk worsens; likely IPW instability"
    } else {
      msg <- "possible implementation issue in pi/ell/objective"
    }
    cat(sprintf("Case 2 [p=%d, gamma=%g, k=%4d]: c1_ratio=%.3f excess_ratio=%.3f → %s\n",
                p_, gm, k_, c1_ratio, ex_ratio, msg))
  }
}

# ── summary 저장 ──────────────────────────────────────────────────────────────

summary_list <- list(
  agg_excess = aggregate(excess_risk ~ p + gamma + case + k + method,
                         data = results_df, FUN = mean, na.rm = TRUE),
  agg_c1     = aggregate(c1_pi ~ p + gamma + case + k + method,
                         data = results_df, FUN = mean, na.rm = TRUE)
)
saveRDS(summary_list, "sim/summary_lowdim_oracle.rds")
cat("\nsummary 저장 완료: sim/summary_lowdim_oracle.rds\n")
