# run_large_N_sim.R
# Large N simulation
# N = 100,000 / Cases: 1, 2, 3 / k: 500, 1000, 2000, 5000 / n0 = 500 / reps = 50 (sanity)
# n_rep = 300으로 올려서 full run 예정

source("code/dgp.R")
source("code/methods.R")
source("code/metrics.R")

# 실험 설정
cases  <- c(1, 2, 3)
k_vals <- c(500, 1000, 2000, 5000)
n_rep  <- 50       # sanity check; full run 시 300으로 변경
N      <- 100000
n0     <- 500      # N 커졌으니 pilot도 비례해서 키움

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
        tryCatch(run_full(dat),                                       error = function(e) NULL),
        tryCatch(run_uni_ipw(dat, k = k_val),                        error = function(e) NULL),
        tryCatch(run_lev_ipw(dat, k = k_val),                        error = function(e) NULL),
        tryCatch(run_opt_homo(dat, k = k_val),                       error = function(e) NULL),
        tryCatch(run_opt_hetero_oracle(dat, k = k_val),              error = function(e) NULL),
        tryCatch(run_opt_hetero_plugin(dat, k = k_val, n0 = n0),     error = function(e) NULL),
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
    cat(sprintf("case %d | k = %5d done\n", case_id, k_val))
  }
}

sim_results <- do.call(rbind, results)

saveRDS(sim_results, "results/large_N_sim_results.rds")
write.csv(sim_results, "results/large_N_sim_results.csv", row.names = FALSE)

cat("Done. Total rows:", nrow(sim_results), "\n")
# 예상 행 수 (sanity): 3 cases × 4 k × 50 reps × 7 methods = 4,200
# 예상 행 수 (full):   3 cases × 4 k × 300 reps × 7 methods = 25,200