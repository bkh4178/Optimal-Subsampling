# run_sim.R
# Sanity-check simulation
# Cases: 1, 2, 3 / k: 200, 500, 1000 / n = 10000 / n0 = 200 / reps = 100

source("code/dgp.R")
source("code/methods.R")
source("code/metrics.R")

# 실험 설정
cases  <- c(1, 2, 3)
k_vals <- c(200, 500, 1000)
n_rep  <- 100
N      <- 10000
n0     <- 200

results <- list()
idx     <- 1

for (case_id in cases) {
  for (k_val in k_vals) {
    for (b in 1:n_rep) {
      
      set.seed(b * 1000 + case_id * 10 + k_val)
      
      # train / test 데이터 생성
      dat      <- generate_data(N = N, case = case_id)
      dat_test <- generate_data(N = N, case = case_id)
      
      X_test <- dat_test$X
      y_test <- dat_test$y
      
      # 각 method 실행 + 지표 계산
      methods_list <- list(
        tryCatch(run_full(dat),                                error = function(e) NULL),
        tryCatch(run_uni_ipw(dat, k = k_val),                 error = function(e) NULL),
        tryCatch(run_lev_ipw(dat, k = k_val),                 error = function(e) NULL),
        tryCatch(run_opt_homo(dat, k = k_val),                error = function(e) NULL),
        tryCatch(run_opt_hetero_oracle(dat, k = k_val),       error = function(e) NULL),
        tryCatch(run_opt_hetero_plugin(dat, k = k_val, n0 = n0), error = function(e) NULL)
      )
      
      for (res in methods_list) {
        if (is.null(res)) next  # 실패한 method skip
        
        m       <- compute_metrics(res, dat, X_test, y_test)
        m$case  <- case_id
        m$k     <- k_val
        m$rep   <- b
        
        results[[idx]] <- m
        idx <- idx + 1
      }
    }
    cat(sprintf("case %d | k = %4d done\n", case_id, k_val))
  }
}

# 결과 합치기
sim_results <- do.call(rbind, results)

# 저장
saveRDS(sim_results, "results/sanity_check_results.rds")
write.csv(sim_results, "results/sanity_check_results.csv", row.names = FALSE)

cat("Done. Total rows:", nrow(sim_results), "\n")




library(dplyr)
sim_results %>%
  group_by(case, k, method) %>%
  summarise(
    mean_excess_risk = mean(excess_risk),
    mean_param_mse   = mean(param_mse),
    mean_n_realized  = mean(n_realized),
    .groups = "drop"
  ) %>%
  arrange(case, k, mean_excess_risk) %>%
  print(n = 54)


set.seed(1)
dat <- generate_data(case = 2)

score_oracle <- sqrt(dat$sigma2 * dat$ell)
summary(score_oracle)
sd(score_oracle)

pi_oracle <- compute_pi(score_oracle, k = 500)
summary(pi_oracle)
mean(pi_oracle >= 1 - 1e-10)  # clip rate





set.seed(1)
dat <- generate_data(case = 2)

# oracle vs homo pi 분포 비교
pi_oracle <- compute_pi(sqrt(dat$sigma2 * dat$ell), k = 500)
pi_homo   <- compute_pi(sqrt(dat$ell), k = 500)

# pi랑 sigma2 상관관계
cat("cor(pi_oracle, sigma2):", cor(pi_oracle, dat$sigma2), "\n")
cat("cor(pi_homo,   sigma2):", cor(pi_homo,   dat$sigma2), "\n")

# 높은 sigma2 관측치에 높은 pi 주는지
cat("cor(pi_oracle, ell):", cor(pi_oracle, dat$ell), "\n")
cat("cor(pi_homo,   ell):", cor(pi_homo,   dat$ell), "\n")

# beta_hat 비교
idx_oracle <- which(rbinom(10000, 1, pi_oracle) == 1)
idx_homo   <- which(rbinom(10000, 1, pi_homo) == 1)

beta_oracle <- fit_wls(dat$X, dat$y, idx_oracle, pi_oracle)
beta_homo   <- fit_wls(dat$X, dat$y, idx_homo,   pi_homo)

cat("param_mse oracle:", sum((beta_oracle - dat$beta0)^2), "\n")
cat("param_mse homo:  ", sum((beta_homo   - dat$beta0)^2), "\n")





set.seed(42)
n_rep <- 200
dat <- generate_data(case = 2)

pm_oracle <- numeric(n_rep)
pm_homo   <- numeric(n_rep)

for (b in 1:n_rep) {
  pi_oracle <- compute_pi(sqrt(dat$sigma2 * dat$ell), k = 500)
  pi_homo   <- compute_pi(sqrt(dat$ell), k = 500)
  
  idx_oracle <- which(rbinom(10000, 1, pi_oracle) == 1)
  idx_homo   <- which(rbinom(10000, 1, pi_homo) == 1)
  
  beta_oracle <- fit_wls(dat$X, dat$y, idx_oracle, pi_oracle)
  beta_homo   <- fit_wls(dat$X, dat$y, idx_homo,   pi_homo)
  
  pm_oracle[b] <- sum((beta_oracle - dat$beta0)^2)
  pm_homo[b]   <- sum((beta_homo   - dat$beta0)^2)
}

cat("oracle: mean =", mean(pm_oracle), "| sd =", sd(pm_oracle), "\n")
cat("homo:   mean =", mean(pm_homo),   "| sd =", sd(pm_homo),   "\n")





# fit_wls 내부 확인
# w = 1 / pi[idx] 로 되어 있는지
# pi는 N × 1 전체 벡터인지

set.seed(1)
dat <- generate_data(case = 2)
pi_oracle <- compute_pi(sqrt(dat$sigma2 * dat$ell), k = 500)
idx <- which(rbinom(10000, 1, pi_oracle) == 1)

# subsample된 관측치의 weight 분포
w <- 1 / pi_oracle[idx]
cat("weight summary:\n")
print(summary(w))
cat("cor(w, sigma2[idx]):", cor(w, dat$sigma2[idx]), "\n")




sim_results %>%
  filter(case == 2) %>%
  group_by(k, method) %>%
  summarise(mean_excess_risk = mean(excess_risk), .groups = "drop") %>%
  arrange(k, mean_excess_risk)


debug_score <- function(dat, k) {
  scores <- list(
    LEV               = dat$ell,
    OPT_homo          = sqrt(dat$ell),
    OPT_hetero_oracle = sqrt(dat$sigma2 * dat$ell)
  )
  out <- lapply(names(scores), function(nm) {
    s  <- scores[[nm]]
    pi <- compute_pi(s, k)
    data.frame(
      method        = nm,
      score_min     = min(s),
      score_q1      = quantile(s, 0.25),
      score_med     = median(s),
      score_mean    = mean(s),
      score_q3      = quantile(s, 0.75),
      score_max     = max(s),
      pi_min        = min(pi),
      pi_q1         = quantile(pi, 0.25),
      pi_med        = median(pi),
      pi_mean       = mean(pi),
      pi_q3         = quantile(pi, 0.75),
      pi_max        = max(pi),
      clip_rate     = mean(pi >= 1 - 1e-10),
      eff_n         = sum(pi)^2 / sum(pi^2),
      cor_pi_sigma2 = cor(pi, dat$sigma2),
      cor_pi_ell    = cor(pi, dat$ell)
    )
  })
  do.call(rbind, out)
}


compare_objective <- function(dat, k) {
  N <- nrow(dat$X)
  pi_uni    <- rep(k / N, N)
  pi_lev    <- compute_pi(dat$ell, k)
  pi_homo   <- compute_pi(sqrt(dat$ell), k)
  pi_oracle <- compute_pi(sqrt(dat$sigma2 * dat$ell), k)
  
  obj <- function(pi) sum(dat$sigma2 * dat$ell / pi)
  
  data.frame(
    method    = c("UNI", "LEV", "OPT-homo", "OPT-hetero-oracle"),
    objective = c(obj(pi_uni), obj(pi_lev), obj(pi_homo), obj(pi_oracle))
  )
}

set.seed(1)
dat2 <- generate_data(case = 2)
dat3 <- generate_data(case = 3)

cat("=== Case 2 objective ===\n")
print(compare_objective(dat2, k = 500))

cat("=== Case 2 score/pi debug ===\n")
print(debug_score(dat2, k = 500))

cat("=== Case 3 objective ===\n")
print(compare_objective(dat3, k = 500))

cat("=== Case 3 score/pi debug ===\n")
print(debug_score(dat3, k = 500))


set.seed(1)
dat2 <- generate_data(case = 2)

for (k_val in c(500, 1000, 2000, 5000)) {
  pi_homo   <- compute_pi(sqrt(dat2$ell), k = k_val)
  pi_oracle <- compute_pi(sqrt(dat2$sigma2 * dat2$ell), k = k_val)
  
  idx_homo   <- which(rbinom(10000, 1, pi_homo) == 1)
  idx_oracle <- which(rbinom(10000, 1, pi_oracle) == 1)
  
  beta_homo   <- fit_wls(dat2$X, dat2$y, idx_homo,   pi_homo)
  beta_oracle <- fit_wls(dat2$X, dat2$y, idx_oracle, pi_oracle)
  
  cat(sprintf("k=%4d | homo: %.4f | oracle: %.4f\n",
              k_val,
              sum((beta_homo   - dat2$beta0)^2),
              sum((beta_oracle - dat2$beta0)^2)))
}




set.seed(42)
dat2 <- generate_data(case = 2)
n_rep <- 500

pm_homo   <- numeric(n_rep)
pm_oracle <- numeric(n_rep)

for (b in 1:n_rep) {
  pi_homo   <- compute_pi(sqrt(dat2$ell), k = 500)
  pi_oracle <- compute_pi(sqrt(dat2$sigma2 * dat2$ell), k = 500)
  
  idx_homo   <- which(rbinom(10000, 1, pi_homo) == 1)
  idx_oracle <- which(rbinom(10000, 1, pi_oracle) == 1)
  
  pm_homo[b]   <- sum((fit_wls(dat2$X, dat2$y, idx_homo,   pi_homo)   - dat2$beta0)^2)
  pm_oracle[b] <- sum((fit_wls(dat2$X, dat2$y, idx_oracle, pi_oracle) - dat2$beta0)^2)
}

cat("homo:   mean =", mean(pm_homo),   "| sd =", sd(pm_homo),   "| median =", median(pm_homo),   "\n")
cat("oracle: mean =", mean(pm_oracle), "| sd =", sd(pm_oracle), "| median =", median(pm_oracle), "\n")


set.seed(1)
dat2 <- generate_data(case = 2)

pi_homo   <- compute_pi(sqrt(dat2$ell), k = 500)
pi_oracle <- compute_pi(sqrt(dat2$sigma2 * dat2$ell), k = 500)

# GLS weight 방향과의 상관
gls_weight <- 1 / dat2$sigma2

cat("cor(1/pi_homo,   gls_weight):", cor(1/pi_homo,   gls_weight), "\n")
cat("cor(1/pi_oracle, gls_weight):", cor(1/pi_oracle, gls_weight), "\n")



set.seed(42)
dat2  <- generate_data(case = 2)
n_rep <- 500
N     <- nrow(dat2$X)
k     <- 500

beta_mat_homo   <- matrix(NA, n_rep, 12)
beta_mat_oracle <- matrix(NA, n_rep, 12)

pi_homo   <- compute_pi(sqrt(dat2$ell), k = k)
pi_oracle <- compute_pi(sqrt(dat2$sigma2 * dat2$ell), k = k)

for (b in 1:n_rep) {
  idx_homo   <- which(rbinom(N, 1, pi_homo) == 1)
  idx_oracle <- which(rbinom(N, 1, pi_oracle) == 1)
  
  beta_mat_homo[b, ]   <- fit_wls(dat2$X, dat2$y, idx_homo,   pi_homo)
  beta_mat_oracle[b, ] <- fit_wls(dat2$X, dat2$y, idx_oracle, pi_oracle)
}

# empirical variance of beta_hat (trace of cov matrix)
emp_var_homo   <- sum(diag(cov(beta_mat_homo)))
emp_var_oracle <- sum(diag(cov(beta_mat_oracle)))

# 이론 예측값: (1/N) * tr(A^{-1} B_pi A^{-1})
A_inv <- solve(dat2$A_hat)
B_homo   <- t(dat2$X) %*% diag(dat2$sigma2 / pi_homo)   %*% dat2$X / N
B_oracle <- t(dat2$X) %*% diag(dat2$sigma2 / pi_oracle) %*% dat2$X / N

theory_var_homo   <- sum(diag(A_inv %*% B_homo   %*% A_inv)) / N
theory_var_oracle <- sum(diag(A_inv %*% B_oracle %*% A_inv)) / N

cat("homo:   empirical =", emp_var_homo,   "| theory =", theory_var_homo,   "\n")
cat("oracle: empirical =", emp_var_oracle, "| theory =", theory_var_oracle, "\n")