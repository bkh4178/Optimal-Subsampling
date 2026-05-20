# analyze_results.R
# Full simulation 결과 분석
# 직접 확인용 — excess_risk 기준 mean / median 요약

library(dplyr)
library(tidyr)

# ── 데이터 로드 ──────────────────────────────────────────────
# df <- read.csv("results/full_sim_results.csv")
df <- read.csv("results/large_N_sim_results.csv")

# method 순서 고정 (표 출력용)
method_order <- c("FULL", "UNI-IPW", "LEV-IPW", "OPT-homo",
                  "OPT-hetero-oracle", "OPT-hetero-plugin", "OPT-hetero-shrink-0.5")
df$method <- factor(df$method, levels = method_order)

# ── 요약 함수 ─────────────────────────────────────────────────
summarize_by <- function(data, metric = "excess_risk") {
  data %>%
    group_by(case, k, method) %>%
    summarise(
      mean   = round(mean(.data[[metric]]), 5),
      median = round(median(.data[[metric]]), 5),
      q90    = round(quantile(.data[[metric]], 0.90), 5),
      sd     = round(sd(.data[[metric]]), 5),
      .groups = "drop"
    )
}

# ── Case별 출력 ───────────────────────────────────────────────
for (c_id in c(1, 2, 3)) {
  cat(sprintf("\n%s\nCASE %d — Excess Risk\n%s\n", strrep("=", 60), c_id, strrep("=", 60)))
  sub <- df %>% filter(case == c_id)
  print(summarize_by(sub, "excess_risk"), n = 100)
}

# ── k별 1위 method (mean excess_risk 기준) ────────────────────
cat("\n\n== 1위 method (mean excess_risk, FULL 제외) ==\n")
df %>%
  filter(method != "FULL") %>%
  group_by(case, k, method) %>%
  summarise(mean_er = mean(excess_risk), .groups = "drop") %>%
  group_by(case, k) %>%
  slice_min(mean_er, n = 1) %>%
  arrange(case, k) %>%
  print(n = 50)

# ── oracle vs plugin vs shrink 비교 (Case 2, 3) ───────────────
cat("\n\n== Oracle / Plugin / Shrink 비교 (Case 2, 3) ==\n")
df %>%
  filter(case %in% c(2, 3),
         method %in% c("OPT-hetero-oracle", "OPT-hetero-plugin", "OPT-hetero-shrink-0.5")) %>%
  group_by(case, k, method) %>%
  summarise(mean_er = round(mean(excess_risk), 5), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = method, values_from = mean_er) %>%
  arrange(case, k) %>%
  print(n = 50)

df <- read.csv("results/large_N_sim_results.csv")

df %>%
  group_by(method, k) %>%
  summarise(
    mean_realized = mean(n_realized),
    sd_realized   = sd(n_realized),
    min_realized  = min(n_realized),
    max_realized  = max(n_realized),
    mean_clip     = mean(clip_rate),
    .groups = "drop"
  ) %>%
  arrange(k, method) %>%
  print(n = Inf)
