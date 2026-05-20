res <- readRDS("sim/results_lowdim7_homo.rds")
sub <- res[res$k == 500 & !res$fail, ]
homo   <- sub$excess_risk[sub$method == "OPT-homo"]
oracle <- sub$excess_risk[sub$method == "OPT-hetero-oracle"]
cat(mean(homo), mean(oracle), "\n")
cat(sd(homo)/sqrt(length(homo)), "\n")
cat(sd(oracle)/sqrt(length(oracle)), "\n")



# run_sim_lowdim7_homo.R에서 oracle score가 실제로 homo랑 같은지
set.seed(1)
dat <- generate_data(10000, c(1,1,1), seed=123)
score_homo   <- sqrt(dat$ell)
score_oracle <- sqrt(dat$sigma2_vec * dat$ell)
all.equal(score_homo, score_oracle)
summary(dat$sigma2_vec)

# n_rep를 2000으로 늘리면 차이가 줄어드는지
mean(oracle) / mean(homo)




res <- readRDS("sim/results_lowdim7_homo.rds")
sub <- res[res$k == 500 & !res$fail, ]

homo   <- sub$excess_risk[sub$method == "OPT-homo"]
oracle <- sub$excess_risk[sub$method == "OPT-hetero-oracle"]

cat("mean homo:  ", mean(homo), "\n")
cat("mean oracle:", mean(oracle), "\n")
cat("ratio:      ", mean(oracle)/mean(homo), "\n")
cat("SE homo:    ", sd(homo)/sqrt(length(homo)), "\n")
cat("SE oracle:  ", sd(oracle)/sqrt(length(oracle)), "\n")
cat("diff / SE:  ", (mean(oracle)-mean(homo)) / sqrt(sd(homo)^2/length(homo) + sd(oracle)^2/length(oracle)), "\n")






res <- readRDS("sim/results_p3_mixed.rds")
names(res)
str(res)


library(dplyr)
install.packages('dplyr')
res %>%
  filter(!fail) %>%
  group_by(method) %>%
  summarise(
    mean_excess = mean(excess_risk, na.rm = TRUE),
    sd_excess   = sd(excess_risk, na.rm = TRUE),
    mcse_excess = sd_excess / sqrt(n()),
    mean_c1     = mean(c1_pi, na.rm = TRUE),
    sd_c1       = sd(c1_pi, na.rm = TRUE),
    mcse_c1     = sd_c1 / sqrt(n()),
    .groups = "drop"
  )



library(dplyr)
library(tidyr)

wide <- res %>%
  filter(!fail) %>%
  select(rep, method, excess_risk, c1_pi) %>%
  pivot_wider(
    names_from = method,
    values_from = c(excess_risk, c1_pi)
  )

wide %>%
  summarise(
    diff_plugin_homo_ex =
      mean(`excess_risk_OPT-hetero-plugin` - `excess_risk_OPT-homo`, na.rm = TRUE),
    mcse_plugin_homo_ex =
      sd(`excess_risk_OPT-hetero-plugin` - `excess_risk_OPT-homo`, na.rm = TRUE) / sqrt(n()),
    
    diff_oracle_homo_ex =
      mean(`excess_risk_OPT-hetero-oracle` - `excess_risk_OPT-homo`, na.rm = TRUE),
    mcse_oracle_homo_ex =
      sd(`excess_risk_OPT-hetero-oracle` - `excess_risk_OPT-homo`, na.rm = TRUE) / sqrt(n()),
    
    diff_plugin_oracle_ex =
      mean(`excess_risk_OPT-hetero-plugin` - `excess_risk_OPT-hetero-oracle`, na.rm = TRUE),
    mcse_plugin_oracle_ex =
      sd(`excess_risk_OPT-hetero-plugin` - `excess_risk_OPT-hetero-oracle`, na.rm = TRUE) / sqrt(n())
  )



wide_ratio <- wide %>%
  mutate(
    plugin_homo_ex_ratio =
      `excess_risk_OPT-hetero-plugin` / `excess_risk_OPT-homo`,
    oracle_homo_ex_ratio =
      `excess_risk_OPT-hetero-oracle` / `excess_risk_OPT-homo`,
    plugin_oracle_ex_ratio =
      `excess_risk_OPT-hetero-plugin` / `excess_risk_OPT-hetero-oracle`,
    
    plugin_homo_c1_ratio =
      `c1_pi_OPT-hetero-plugin` / `c1_pi_OPT-homo`,
    oracle_homo_c1_ratio =
      `c1_pi_OPT-hetero-oracle` / `c1_pi_OPT-homo`,
    plugin_oracle_c1_ratio =
      `c1_pi_OPT-hetero-plugin` / `c1_pi_OPT-hetero-oracle`
  )

summary(wide_ratio$plugin_homo_ex_ratio)
summary(wide_ratio$plugin_homo_c1_ratio)




res %>%
  filter(method == "OPT-hetero-plugin", !fail) %>%
  summarise(
    mean_sigma_cor = mean(diag_sigma_cor, na.rm = TRUE),
    sd_sigma_cor   = sd(diag_sigma_cor, na.rm = TRUE),
    mean_score_cor = mean(diag_score_cor, na.rm = TRUE),
    sd_score_cor   = sd(diag_score_cor, na.rm = TRUE),
    min_score_cor  = min(diag_score_cor, na.rm = TRUE),
    q05_score_cor  = quantile(diag_score_cor, 0.05, na.rm = TRUE),
    q95_score_cor  = quantile(diag_score_cor, 0.95, na.rm = TRUE)
  )


install.packages('tidyr')
library(dplyr)
library(tidyr)

res <- readRDS("sim/results_p9.rds")

wide <- res %>%
  filter(!fail) %>%
  select(rep, method, excess_risk, c1_pi) %>%
  pivot_wider(
    names_from = method,
    values_from = c(excess_risk, c1_pi)
  )

wide %>%
  summarise(
    diff_plugin_homo_ex =
      mean(`excess_risk_OPT-hetero-plugin` - `excess_risk_OPT-homo`, na.rm = TRUE),
    mcse_plugin_homo_ex =
      sd(`excess_risk_OPT-hetero-plugin` - `excess_risk_OPT-homo`, na.rm = TRUE) / sqrt(n()),

    diff_oracle_homo_ex =
      mean(`excess_risk_OPT-hetero-oracle` - `excess_risk_OPT-homo`, na.rm = TRUE),
    mcse_oracle_homo_ex =
      sd(`excess_risk_OPT-hetero-oracle` - `excess_risk_OPT-homo`, na.rm = TRUE) / sqrt(n()),

    diff_plugin_oracle_ex =
      mean(`excess_risk_OPT-hetero-plugin` - `excess_risk_OPT-hetero-oracle`, na.rm = TRUE),
    mcse_plugin_oracle_ex =
      sd(`excess_risk_OPT-hetero-plugin` - `excess_risk_OPT-hetero-oracle`, na.rm = TRUE) / sqrt(n())
  )




library(dplyr)
library(tidyr)

res <- readRDS("sim/results_p9_mixed.rds")

wide <- res %>%
  filter(!fail) %>%
  select(rep, method, excess_risk, c1_pi) %>%
  pivot_wider(
    names_from = method,
    values_from = c(excess_risk, c1_pi)
  )

wide %>%
  summarise(
    diff_plugin_homo_ex =
      mean(`excess_risk_OPT-hetero-plugin` - `excess_risk_OPT-homo`, na.rm = TRUE),
    mcse_plugin_homo_ex =
      sd(`excess_risk_OPT-hetero-plugin` - `excess_risk_OPT-homo`, na.rm = TRUE) / sqrt(n()),

    diff_oracle_homo_ex =
      mean(`excess_risk_OPT-hetero-oracle` - `excess_risk_OPT-homo`, na.rm = TRUE),
    mcse_oracle_homo_ex =
      sd(`excess_risk_OPT-hetero-oracle` - `excess_risk_OPT-homo`, na.rm = TRUE) / sqrt(n()),

    diff_plugin_oracle_ex =
      mean(`excess_risk_OPT-hetero-plugin` - `excess_risk_OPT-hetero-oracle`, na.rm = TRUE),
    mcse_plugin_oracle_ex =
      sd(`excess_risk_OPT-hetero-plugin` - `excess_risk_OPT-hetero-oracle`, na.rm = TRUE) / sqrt(n())
  )
