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




