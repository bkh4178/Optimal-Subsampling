# 06_runtime_utils.R
# 런타임 측정 유틸리티
# 의존: 04_metrics.R (run_one_rep)

# 단일 rep 전체 경과 시간 (초)
time_full_rep <- function(dat, k, X_test, mu_test, y_test, n0) {
  t0 <- proc.time()
  run_one_rep(dat, k, X_test, mu_test, y_test, n0)
  (proc.time() - t0)[["elapsed"]]
}

# n_bench 번 반복 측정 → 평균 / SD / median 반환
benchmark_rep <- function(dat, k, X_test, mu_test, y_test, n0,
                          n_bench = 50L, seed_start = 2001L) {
  # 워밍업 (측정 제외)
  for (w in 1:3) time_full_rep(dat, k, X_test, mu_test, y_test, n0)
  
  times <- numeric(n_bench)
  for (i in seq_len(n_bench)) {
    set.seed(seed_start + i)
    times[i] <- time_full_rep(dat, k, X_test, mu_test, y_test, n0)
    if (i %% 10L == 0L) cat(sprintf("  bench %d/%d: %.3fs\n", i, n_bench, times[i]))
  }
  list(
    mean   = mean(times),
    sd     = sd(times),
    median = median(times),
    times  = times,
    n_bench = n_bench
  )
}
