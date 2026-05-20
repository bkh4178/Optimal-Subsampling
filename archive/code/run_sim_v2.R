# run_sim_v2.R
# Full simulation v2
# Cases: 1, 2 only / k: 200, 500, 1000, 2000 / N = 10000 / n0 = 200 / reps = 300
# DGP: p=9, Normal 3 / Discrete Uniform 3 / Bernoulli 3 (lognormal 제거)
# Case 3 제외 (variance ⊥ leverage 설계 논리 문제로 보류)

source("code/dgp.R")
source("code/methods.R")
source("code/metrics.R")

# 실험 설정
cases  <- c(1, 2)
k_vals <- c(200, 500, 1000, 2000)
n_rep  <- 300
N      <- 10000
n0     <- 200

results <- list()
idx     <- 1

for (case_id in cases) {
  for (k_val in k_vals) {
    for (b in 1:n_rep) {
      
      set.seed(b * 1000 + case_id * 10 + k_val)
      
      dat      <- generate_data(N = N, case = case_id)
      dat_test <- generate_data(N = N, case = case_id)
      
      X_test <- dat_test$X
      y_test <- dat_test$y
      
      methods_list <- list(
        tryCatch(run_full(dat),                                        error = function(e) NULL),
        tryCatch(run_uni_ipw(dat, k = k_val),                         error = function(e) NULL),
        tryCatch(run_lev_ipw(dat, k = k_val),                         error = function(e) NULL),
        tryCatch(run_opt_homo(dat, k = k_val),                        error = function(e) NULL),
        tryCatch(run_opt_hetero_oracle(dat, k = k_val),               error = function(e) NULL),
        tryCatch(run_opt_hetero_plugin(dat, k = k_val, n0 = n0),      error = function(e) NULL),
        tryCatch(run_opt_hetero_shrink(dat, k = k_val, lambda = 0.5), error = function(e) NULL)
      )
      
      for (res in methods_list) {
        if (is.null(res)) next
        
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

sim_results <- do.call(rbind, results)

saveRDS(sim_results, "results/sim_v2_results.rds")
write.csv(sim_results, "results/sim_v2_results.csv", row.names = FALSE)

cat("Done. Total rows:", nrow(sim_results), "\n")
# 예상 행 수: 2 cases × 4 k × 300 reps × 7 methods = 16,800