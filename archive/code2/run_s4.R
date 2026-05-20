# run_s4.R
# Simulation Section v3 — S4: Continuous Separated Heteroskedasticity
#
# DGP: p=9 (Normal 3 / Discrete Uniform 3 / Bernoulli 3), 기존 구조 유지
#      sigma2(X) = exp(0.5 * X[,1]) — Normal 첫 번째 열
#      X[,1] ~ N(0,1): leverage 기여가 있지만 DUnif/Bernoulli와 혼재
#      → S3(Bernoulli 열 기반 group-wise)와 달리 continuous covariate 기반 이분산
#
# Goal: S3 결과의 robustness 확인
#       "Bernoulli 열 결과가 continuous setting에서도 유지되는가?"
#
# N = 10000 / k ∈ {200, 500, 1000, 2000} / n_rep = 300 / n0 = 200

source("code2/methods.R")
source("code2/metrics.R")

# ── S4 전용 DGP ───────────────────────────────────────────────────────────────
generate_data_s4 <- function(N = 10000, beta0 = NULL) {

  p <- 9
  if (is.null(beta0)) beta0 <- rep(1, p)

  # 공변량 생성 (기존 3/3/3 구조 동일)
  X_norm  <- matrix(rnorm(N * 3), N, 3)
  X_dunif <- matrix(sample(1:5, N * 3, replace = TRUE), N, 3)
  X_dunif <- sweep(X_dunif, 2, colMeans(X_dunif), FUN = "-")
  X_bern  <- matrix(rbinom(N * 3, 1, 0.5), N, 3)
  X_bern  <- sweep(X_bern, 2, colMeans(X_bern), FUN = "-")

  X <- cbind(X_norm, X_dunif, X_bern)

  # A_hat 및 leverage score
  A_hat     <- crossprod(X) / N
  A_hat_inv <- solve(A_hat)
  ell       <- rowSums((X %*% A_hat_inv) * X)
  h         <- ell / N

  ell_rank  <- rank(ell) / (N + 1)
  ell_score <- qnorm(ell_rank)
  ell_score <- pmax(pmin(ell_score, 2), -2)

  # S4 sigma2: continuous, Normal 첫 번째 열(X[,1]) 기반
  # X[,1] ~ N(0,1) — leverage와 부분적 연관이 있지만
  # DUnif/Bernoulli 열들과 혼재하므로 S2보다 훨씬 약한 alignment
  sigma2_raw <- exp(0.5 * X[, 1])
  sigma2     <- as.vector(sigma2_raw / mean(sigma2_raw))

  mu      <- as.vector(X %*% beta0)
  epsilon <- rnorm(N, mean = 0, sd = sqrt(sigma2))
  y       <- mu + epsilon

  list(
    X = X, y = y, mu = mu,
    sigma2 = sigma2, ell = ell, h = h,
    A_hat = A_hat, beta0 = beta0,
    case = "4"
  )
}

# ── 실험 설정 ─────────────────────────────────────────────────────────────────
k_vals <- c(200, 500, 1000, 2000)
n_rep  <- 300
N      <- 10000
n0     <- 200

results <- list()
idx     <- 1

for (k_val in k_vals) {
  for (b in 1:n_rep) {

    set.seed(b * 1000 + 40 + k_val)   # case=4 고정

    dat      <- generate_data_s4(N = N)
    dat_test <- generate_data_s4(N = N)

    X_test <- dat_test$X
    y_test <- dat_test$y

    methods_list <- list(
      tryCatch(run_full(dat),                                          error = function(e) NULL),
      tryCatch(run_uni_ipw(dat, k = k_val),                           error = function(e) NULL),
      tryCatch(run_lev_ipw(dat, k = k_val),                           error = function(e) NULL),
      tryCatch(run_opt_homo(dat, k = k_val),                          error = function(e) NULL),
      tryCatch(run_opt_hetero_oracle(dat, k = k_val),                 error = function(e) NULL),
      tryCatch(run_opt_hetero_plugin(dat, k = k_val, n0 = n0),        error = function(e) NULL),
      tryCatch(run_opt_hetero_tempering(dat, k = k_val, gamma = 0.5), error = function(e) NULL)
    )

    for (res in methods_list) {
      if (is.null(res)) next

      m        <- compute_metrics(res, dat, X_test, y_test)
      m$case   <- 4
      m$k      <- k_val
      m$rep    <- b

      results[[idx]] <- m
      idx <- idx + 1
    }
  }
  cat(sprintf("S4 | k = %4d done\n", k_val))
}

sim_s4 <- do.call(rbind, results)

saveRDS(sim_s4, "results/sim_s4_results.rds")
write.csv(sim_s4, "results/sim_s4_results.csv", row.names = FALSE)

cat("Done. Total rows:", nrow(sim_s4), "\n")
# 예상 행 수: 4 k × 300 reps × 7 methods = 8,400

# ── 터미널 1차 확인 ───────────────────────────────────────────────────────────

method_order <- c("FULL", "UNI-IPW", "LEV-IPW", "OPT-homo",
                  "OPT-hetero-oracle", "OPT-hetero-plugin",
                  "OPT-hetero-tempering-0.5")
sim_s4$method <- factor(sim_s4$method, levels = method_order)

# 1. k별 method 순위 (mean excess risk 기준, FULL 제외)
cat("\n\n========== [S4] 1. k별 순위 (mean excess risk, FULL 제외) ==========\n")
agg <- aggregate(excess_risk ~ k + method, data = sim_s4[sim_s4$method != "FULL", ],
                 FUN = mean)
agg <- agg[order(agg$k, agg$excess_risk), ]
for (kv in k_vals) {
  sub <- agg[agg$k == kv, ]
  cat(sprintf("\n  k = %4d\n", kv))
  for (i in seq_len(nrow(sub))) {
    cat(sprintf("    %d. %-35s %.5f\n", i, as.character(sub$method[i]), sub$excess_risk[i]))
  }
}

# 2. Tail diagnostic — mean vs. median
cat("\n\n========== [S4] 2. Tail diagnostic — mean vs. median ==========\n")
cat(sprintf("  %-35s %6s %8s %8s %8s\n", "method", "k", "mean", "median", "mean/med"))
for (kv in k_vals) {
  for (m in levels(sim_s4$method)[-1]) {
    er <- sim_s4$excess_risk[sim_s4$method == m & sim_s4$k == kv]
    if (length(er) == 0) next
    cat(sprintf("  %-35s %6d %8.5f %8.5f %8.3f\n",
                m, kv, mean(er), median(er), mean(er) / median(er)))
  }
}

# 3. π concentration diagnostic
cat("\n\n========== [S4] 3. π concentration diagnostic ==========\n")
cat(sprintf("  %-35s %6s %8s %10s %10s %10s\n",
            "method", "k", "clip_rt", "ESS_pi", "top1%_mass", "max/min_pi"))

for (kv in k_vals) {
  set.seed(1 * 1000 + 40 + kv)
  dat_diag <- generate_data_s4(N = N)

  pi_list <- list(
    "UNI-IPW"                  = rep(kv / N, N),
    "LEV-IPW"                  = compute_pi(dat_diag$ell, k = kv),
    "OPT-homo"                 = compute_pi(sqrt(dat_diag$ell), k = kv),
    "OPT-hetero-oracle"        = compute_pi(sqrt(dat_diag$sigma2 * dat_diag$ell), k = kv),
    "OPT-hetero-tempering-0.5" = compute_pi(sqrt(dat_diag$sigma2 * dat_diag$ell)^0.5, k = kv)
  )

  for (nm in names(pi_list)) {
    pv      <- pi_list[[nm]]
    ess     <- sum(pv)^2 / sum(pv^2)
    top_idx <- pv > quantile(pv, 0.99)
    if (sum(top_idx) == 0) {
      top_mass <- 0.01
    } else {
      top_mass <- sum(pv[top_idx]) / sum(pv)
    }
    mm_ratio <- max(pv) / (min(pv[pv > 0]) + 1e-15)
    clip_rt  <- mean(pv >= 1 - 1e-10)
    cat(sprintf("  %-35s %6d %8.4f %10.1f %10.4f %10.1f\n",
                nm, kv, clip_rt, ess, top_mass, mm_ratio))
  }
}

# 4. 분포 진단 — mean / median / sd / q90
cat("\n\n========== [S4] 4. 분포 진단 (excess risk) ==========\n")
for (kv in k_vals) {
  cat(sprintf("\n  k = %4d\n", kv))
  sub <- sim_s4[sim_s4$method != "FULL" & sim_s4$k == kv, ]
  cat(sprintf("  %-35s %8s %8s %8s %8s\n", "method", "mean", "median", "sd", "q90"))
  for (m in levels(sim_s4$method)[-1]) {
    er <- sub$excess_risk[sub$method == m]
    if (length(er) == 0) next
    cat(sprintf("  %-35s %8.5f %8.5f %8.5f %8.5f\n",
                m, mean(er), median(er), sd(er), quantile(er, 0.90)))
  }
}