run_opt_hetero_oracle_shrink <- function(dat, k, lambda) {
  N <- nrow(dat$X)
  
  pi_oracle <- compute_pi(sqrt(dat$sigma2 * dat$ell), k = k)
  pi_shrink <- (1 - lambda) * pi_oracle + lambda * (k / N)
  
  # shrink 후에도 sum = k 유지 확인 (선형 결합이라 자동으로 성립)
  idx <- which(rbinom(N, 1, pi_shrink) == 1)
  
  list(
    method     = paste0("OPT-hetero-oracle-shrink-", lambda),
    beta_hat   = fit_wls(dat$X, dat$y, idx = idx, pi = pi_shrink),
    pi         = pi_shrink,
    idx        = idx,
    n_realized = length(idx),
    clip_rate  = mean(pi_shrink >= 1 - 1e-10)
  )
}


set.seed(42)
dat2  <- generate_data(case = 2)
n_rep <- 500
N     <- nrow(dat2$X)
k     <- 500

lambdas <- c(0, 0.1, 0.2, 0.5)  # 0 = 원래 oracle
results <- list()

for (lam in lambdas) {
  pm <- numeric(n_rep)
  er <- numeric(n_rep)
  
  dat_test <- generate_data(case = 2)
  
  for (b in 1:n_rep) {
    res <- run_opt_hetero_oracle_shrink(dat2, k = k, lambda = lam)
    pm[b] <- compute_param_mse(res$beta_hat, dat2$beta0)
    er[b] <- compute_excess_risk(res$beta_hat, dat2$beta0, dat_test$X)
  }
  
  results[[as.character(lam)]] <- data.frame(
    lambda      = lam,
    mean_pm     = mean(pm),
    median_pm   = median(pm),
    sd_pm       = sd(pm),
    mean_er     = mean(er),
    median_er   = median(er)
  )
  
  cat(sprintf("lambda = %.1f | param_mse mean = %.4f | excess_risk mean = %.4f\n",
              lam, mean(pm), mean(er)))
}

do.call(rbind, results)





set.seed(42)
dat2 <- generate_data(case = 2)
dat_test <- generate_data(case = 2)
n_rep <- 500

pm_homo <- numeric(n_rep)
for (b in 1:n_rep) {
  pi_homo <- compute_pi(sqrt(dat2$ell), k = 500)
  idx     <- which(rbinom(nrow(dat2$X), 1, pi_homo) == 1)
  beta    <- fit_wls(dat2$X, dat2$y, idx, pi_homo)
  pm_homo[b] <- compute_excess_risk(beta, dat2$beta0, dat_test$X)
}
cat("OPT-homo excess_risk mean:", mean(pm_homo), "\n")







set.seed(42)
N      <- 10000
k_vals <- c(200, 500, 1000, 2000)
gammas <- c(0.2, 0.4, 0.7, 1.0)
n_rep  <- 300

gamma_results <- list()
idx_g <- 1

for (gam in gammas) {
  
  dat      <- generate_data(case = 2, gamma = gam)
  dat_test <- generate_data(case = 2, gamma = gam)
  
  # diagnostic: sigma2 구조 확인
  cat(sprintf("\n--- gamma = %.1f ---\n", gam))
  cat("cor(sigma2, ell):", round(cor(dat$sigma2, dat$ell), 4), "\n")
  cat("sd(sigma2):", round(sd(dat$sigma2), 4), "\n")
  
  for (k_val in k_vals) {
    
    er_homo   <- numeric(n_rep)
    er_oracle <- numeric(n_rep)
    
    pi_homo   <- compute_pi(sqrt(dat$ell), k = k_val)
    pi_oracle <- compute_pi(sqrt(dat$sigma2 * dat$ell), k = k_val)
    
    eff_n_homo   <- sum(pi_homo)^2   / sum(pi_homo^2)
    eff_n_oracle <- sum(pi_oracle)^2 / sum(pi_oracle^2)
    
    for (b in 1:n_rep) {
      idx_h <- which(rbinom(N, 1, pi_homo) == 1)
      idx_o <- which(rbinom(N, 1, pi_oracle) == 1)
      
      beta_h <- fit_wls(dat$X, dat$y, idx_h, pi_homo)
      beta_o <- fit_wls(dat$X, dat$y, idx_o, pi_oracle)
      
      er_homo[b]   <- compute_excess_risk(beta_h, dat$beta0, dat_test$X)
      er_oracle[b] <- compute_excess_risk(beta_o, dat$beta0, dat_test$X)
    }
    
    gamma_results[[idx_g]] <- data.frame(
      gamma        = gam,
      k            = k_val,
      cor_sigma_ell = round(cor(dat$sigma2, dat$ell), 4),
      eff_n_homo   = round(eff_n_homo, 1),
      eff_n_oracle = round(eff_n_oracle, 1),
      homo_mean    = mean(er_homo),
      homo_median  = median(er_homo),
      homo_q90     = quantile(er_homo, 0.90),
      oracle_mean  = mean(er_oracle),
      oracle_median = median(er_oracle),
      oracle_q90   = quantile(er_oracle, 0.90),
      oracle_wins  = mean(er_oracle) < mean(er_homo)
    )
    idx_g <- idx_g + 1
  }
}

res_gamma <- do.call(rbind, gamma_results)
print(res_gamma)






