# shrinkage_experiment.R
# Case 2, gamma sensitivity + shrinkage diagnostic
# methods: OPT-homo, OPT-hetero-oracle, oracle-shrink (lambda = 0.1, 0.2, 0.5)
# metrics: mean/median/q90/q95/q99 of excess risk
# diagnostics: eff_n, rel_norm_As_minus_A, condition number

source("dgp.R")
source("methods.R")
source("metrics.R")

set.seed(42)
N      <- 10000
k_vals <- c(200, 500, 1000, 2000)
gammas <- c(0.2, 0.4, 0.7, 1.0)
lambdas <- c(0.1, 0.2, 0.5)
n_rep  <- 300

# 공통 유틸: design error 계산
compute_design_diagnostics <- function(dat, idx, pi) {
  X     <- dat$X
  N     <- nrow(X)
  A_hat <- dat$A_hat
  w     <- 1 / pi[idx]
  Xs    <- X[idx, , drop = FALSE]
  A_s   <- crossprod(Xs * sqrt(w)) / N
  list(
    rel_norm  = norm(A_s - A_hat, type = "F") / norm(A_hat, type = "F"),
    cond_As   = kappa(A_s, exact = FALSE)
  )
}

all_results <- list()
idx_r <- 1

for (gam in gammas) {
  
  cat(sprintf("\n===== gamma = %.1f =====\n", gam))
  
  dat      <- generate_data(case = 2, gamma = gam)
  dat_test <- generate_data(case = 2, gamma = gam)
  
  # pi 미리 계산 (데이터 고정, sampling만 반복)
  pi_homo   <- compute_pi(sqrt(dat$ell), k = 500)  # k별로 아래서 재계산
  
  for (k_val in k_vals) {
    
    cat(sprintf("  k = %d\n", k_val))
    
    # pi 계산
    pi_homo   <- compute_pi(sqrt(dat$ell), k = k_val)
    pi_oracle <- compute_pi(sqrt(dat$sigma2 * dat$ell), k = k_val)
    pi_shrink <- lapply(lambdas, function(lam)
      (1 - lam) * pi_oracle + lam * (k_val / N)
    )
    names(pi_shrink) <- paste0("shrink_", lambdas)
    
    # eff_n
    eff_n_homo   <- sum(pi_homo)^2   / sum(pi_homo^2)
    eff_n_oracle <- sum(pi_oracle)^2 / sum(pi_oracle^2)
    eff_n_shrink <- sapply(pi_shrink, function(p) sum(p)^2 / sum(p^2))
    
    # 결과 저장용
    er_list <- list(
      homo   = numeric(n_rep),
      oracle = numeric(n_rep)
    )
    for (lam in lambdas) {
      er_list[[paste0("shrink_", lam)]] <- numeric(n_rep)
    }
    
    dd_homo <- dd_oracle <- vector("list", n_rep)
    
    for (b in 1:n_rep) {
      
      # homo
      idx_h  <- which(rbinom(N, 1, pi_homo) == 1)
      beta_h <- fit_wls(dat$X, dat$y, idx_h, pi_homo)
      er_list$homo[b] <- compute_excess_risk(beta_h, dat$beta0, dat_test$X)
      dd_homo[[b]]    <- compute_design_diagnostics(dat, idx_h, pi_homo)
      
      # oracle
      idx_o  <- which(rbinom(N, 1, pi_oracle) == 1)
      beta_o <- fit_wls(dat$X, dat$y, idx_o, pi_oracle)
      er_list$oracle[b] <- compute_excess_risk(beta_o, dat$beta0, dat_test$X)
      dd_oracle[[b]]    <- compute_design_diagnostics(dat, idx_o, pi_oracle)
      
      # shrink versions
      for (lam in lambdas) {
        ps     <- pi_shrink[[paste0("shrink_", lam)]]
        idx_s  <- which(rbinom(N, 1, ps) == 1)
        beta_s <- fit_wls(dat$X, dat$y, idx_s, ps)
        er_list[[paste0("shrink_", lam)]][b] <-
          compute_excess_risk(beta_s, dat$beta0, dat_test$X)
      }
    }
    
    # design diagnostic 평균
    mean_rel_homo   <- mean(sapply(dd_homo,   function(d) d$rel_norm))
    mean_cond_homo  <- mean(sapply(dd_homo,   function(d) d$cond_As))
    mean_rel_oracle <- mean(sapply(dd_oracle, function(d) d$rel_norm))
    mean_cond_oracle<- mean(sapply(dd_oracle, function(d) d$cond_As))
    
    # 결과 저장
    method_names <- c("homo", "oracle", paste0("shrink_", lambdas))
    
    for (mn in method_names) {
      er <- er_list[[mn]]
      
      if (mn == "homo") {
        eff_n_val  <- eff_n_homo
        rel_norm   <- mean_rel_homo
        cond_val   <- mean_cond_homo
        lam_val    <- NA
      } else if (mn == "oracle") {
        eff_n_val  <- eff_n_oracle
        rel_norm   <- mean_rel_oracle
        cond_val   <- mean_cond_oracle
        lam_val    <- 0
      } else {
        lam_val    <- as.numeric(sub("shrink_", "", mn))
        ps         <- pi_shrink[[mn]]
        eff_n_val  <- sum(ps)^2 / sum(ps^2)
        rel_norm   <- NA  # shrink는 별도 측정 생략
        cond_val   <- NA
        lam_val    <- lam_val
      }
      
      all_results[[idx_r]] <- data.frame(
        gamma         = gam,
        k             = k_val,
        method        = mn,
        lambda        = lam_val,
        cor_sigma_ell = round(cor(dat$sigma2, dat$ell), 4),
        sd_sigma2     = round(sd(dat$sigma2), 4),
        eff_n         = round(eff_n_val, 1),
        mean_er       = mean(er),
        median_er     = median(er),
        q90_er        = quantile(er, 0.90),
        q95_er        = quantile(er, 0.95),
        q99_er        = quantile(er, 0.99),
        mean_rel_norm = rel_norm,
        mean_cond     = cond_val
      )
      idx_r <- idx_r + 1
    }
  }
}

shrink_results <- do.call(rbind, all_results)
saveRDS(shrink_results, "shrinkage_experiment_results.rds")
write.csv(shrink_results, "shrinkage_experiment_results.csv", row.names = FALSE)
cat("\nDone. Total rows:", nrow(shrink_results), "\n")