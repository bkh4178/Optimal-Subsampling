# 06_figures.R
# Generate figures from simulation_2plugin results
# Run from: simulation_2plugin/
# Requires: ggplot2

library(ggplot2)

dir_main <- file.path("figures", "main")
if (!dir.exists(dir_main)) dir.create(dir_main, recursive = TRUE)

compute_mcse <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))

# ── Method aesthetics (전체 figure 공통) ──────────────────────────────────────
METHOD_LINES <- c("SRS", "LEV-IPW", "IBOSS", "OPT-homo",
                   "OPT-hetero-oracle",
                   "OPT-hetero-plugin-param", "OPT-hetero-plugin-nonparam")

METHOD_COLORS <- c(
  "SRS"                        = "#999999",
  "LEV-IPW"                    = "#0072B2",
  "IBOSS"                      = "#56B4E9",
  "OPT-homo"                   = "#009E73",
  "OPT-hetero-oracle"          = "#D55E00",
  "OPT-hetero-plugin-param"    = "#CC79A7",
  "OPT-hetero-plugin-nonparam" = "#7570B3"
)

METHOD_SHAPES <- c(
  "SRS"                        = 16L,
  "LEV-IPW"                    = 17L,
  "IBOSS"                      = 4L,
  "OPT-homo"                   = 15L,
  "OPT-hetero-oracle"          = 18L,
  "OPT-hetero-plugin-param"    = 8L,
  "OPT-hetero-plugin-nonparam" = 3L
)

METHOD_LABELS <- c(
  "SRS"                        = "SRS",
  "LEV-IPW"                    = "LEV-IPW",
  "IBOSS"                      = "IBOSS",
  "OPT-homo"                   = "OPT-homo",
  "OPT-hetero-oracle"          = "OPT-oracle",
  "OPT-hetero-plugin-param"    = "OPT-plugin (param)",
  "OPT-hetero-plugin-nonparam" = "OPT-plugin (RF)"
)

# 공통 요약 헬퍼: method × xvar 기준 mean/CI
summarize_metric <- function(df, methods, xvar, yvar) {
  do.call(rbind, lapply(methods, function(m) {
    sub <- df[df$method == m, ]
    if (nrow(sub) == 0L) return(NULL)
    xs <- sort(unique(sub[[xvar]]))
    do.call(rbind, lapply(xs, function(xv) {
      x  <- sub[[yvar]][sub[[xvar]] == xv]
      x  <- x[is.finite(x)]
      if (length(x) == 0L) return(NULL)
      mn <- mean(x); se <- compute_mcse(x)
      data.frame(method = m, x = xv, mean = mn,
                 lo = mn - 1.96 * se, hi = mn + 1.96 * se,
                 stringsAsFactors = FALSE)
    }))
  }))
}

# ═════════════════════════════════════════════════════════════════════════
# Figure 1 — Main performance: k vs pred_mse_N (전 방법)
# ═════════════════════════════════════════════════════════════════════════
cat("Building Figure 1 (main performance) ...\n")

exp1 <- readRDS(file.path("results", "raw", "exp1", "exp1_k_comparison.rds"))
exp1v <- exp1[!exp1$fail, ]

full_mean <- mean(exp1v$pred_mse_N[exp1v$method == "FULL"], na.rm = TRUE)

sum1 <- summarize_metric(exp1v, METHOD_LINES, "k", "pred_mse_N")
sum1$method <- factor(sum1$method, levels = METHOD_LINES)
k_vals <- sort(unique(sum1$x))

p1 <- ggplot(sum1, aes(x = x, y = mean, color = method, shape = method)) +
  geom_hline(yintercept = full_mean, linetype = "dashed",
             color = "black", linewidth = 0.5) +
  annotate("text", x = min(k_vals), y = full_mean,
           label = "FULL OLS", hjust = 0, vjust = -0.6, size = 3) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = lo, ymax = hi),
                width = diff(range(k_vals)) * 0.02, linewidth = 0.45) +
  scale_color_manual(values = METHOD_COLORS[METHOD_LINES],
                      labels = METHOD_LABELS[METHOD_LINES]) +
  scale_shape_manual(values = METHOD_SHAPES[METHOD_LINES],
                      labels = METHOD_LABELS[METHOD_LINES]) +
  scale_x_continuous(breaks = k_vals) +
  labs(x = "Subsample size k", y = "Predictive Excess Risk (mean ± 1.96 MCSE)",
       color = "Method", shape = "Method") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", legend.title = element_text(size = 10),
        panel.grid.minor = element_blank())

ggsave(file.path(dir_main, "fig1_main_performance.pdf"), p1, width = 7.5, height = 4.8)
cat("  Saved: figures/main/fig1_main_performance.pdf\n")

# ═════════════════════════════════════════════════════════════════════════
# Figure 2 — Leading term (Metric 1) vs Actual ER (Metric 2)
# ═════════════════════════════════════════════════════════════════════════
cat("Building Figure 2 (leading vs actual ER) ...\n")

sum2a <- summarize_metric(exp1v, METHOD_LINES, "k", "er_leading")
sum2a$metric <- "Leading term (theoretical)"
sum2b <- summarize_metric(exp1v, METHOD_LINES, "k", "pred_mse_N")
sum2b$metric <- "Actual ER (realized)"

sum2 <- rbind(sum2a, sum2b)
sum2$method <- factor(sum2$method, levels = METHOD_LINES)
sum2$metric <- factor(sum2$metric,
                       levels = c("Leading term (theoretical)", "Actual ER (realized)"))

p2 <- ggplot(sum2, aes(x = x, y = mean, color = method, shape = method)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.2) +
  geom_errorbar(aes(ymin = lo, ymax = hi),
                width = diff(range(k_vals)) * 0.02, linewidth = 0.4) +
  facet_wrap(~metric, nrow = 1) +
  scale_color_manual(values = METHOD_COLORS[METHOD_LINES],
                      labels = METHOD_LABELS[METHOD_LINES]) +
  scale_shape_manual(values = METHOD_SHAPES[METHOD_LINES],
                      labels = METHOD_LABELS[METHOD_LINES]) +
  scale_x_continuous(breaks = k_vals) +
  labs(x = "Subsample size k", y = "Excess Risk (mean ± 1.96 MCSE)",
       color = "Method", shape = "Method") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        strip.text = element_text(size = 10, face = "bold"))

ggsave(file.path(dir_main, "fig2_leading_vs_actual.pdf"), p2, width = 10, height = 4.8)
cat("  Saved: figures/main/fig2_leading_vs_actual.pdf\n")

# ═════════════════════════════════════════════════════════════════════════
# Figure 3 — n0 sensitivity: oracle / plugin-param / plugin-nonparam
# ═════════════════════════════════════════════════════════════════════════
cat("Building Figure 3 (n0 sensitivity) ...\n")

exp2 <- readRDS(file.path("results", "raw", "exp2", "exp2_n0_sensitivity.rds"))
exp2v <- exp2[!exp2$fail, ]

FIG3_METHODS <- c("OPT-hetero-oracle",
                   "OPT-hetero-plugin-param", "OPT-hetero-plugin-nonparam")

sum3 <- summarize_metric(exp2v, FIG3_METHODS, "n0", "pred_mse_N")
sum3$method <- factor(sum3$method, levels = FIG3_METHODS)
n0_vals <- sort(unique(sum3$x))

p3 <- ggplot(sum3, aes(x = x, y = mean, color = method, shape = method)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = lo, ymax = hi),
                width = diff(range(n0_vals)) * 0.02, linewidth = 0.45) +
  scale_color_manual(values = METHOD_COLORS[FIG3_METHODS],
                      labels = METHOD_LABELS[FIG3_METHODS]) +
  scale_shape_manual(values = METHOD_SHAPES[FIG3_METHODS],
                      labels = METHOD_LABELS[FIG3_METHODS]) +
  scale_x_continuous(breaks = n0_vals) +
  labs(x = "Pilot size n0", y = "Predictive Excess Risk (mean ± 1.96 MCSE)",
       color = "Method", shape = "Method") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave(file.path(dir_main, "fig3_n0_sensitivity.pdf"), p3, width = 7, height = 4.5)
cat("  Saved: figures/main/fig3_n0_sensitivity.pdf\n")

# ═════════════════════════════════════════════════════════════════════════
# Figure 4 — Plugin estimation quality: sigma_cor & score_cor vs n0 (2 subplots)
# ═════════════════════════════════════════════════════════════════════════
cat("Building Figure 4 (plugin estimation quality) ...\n")

FIG4_METHODS <- c("OPT-hetero-plugin-param", "OPT-hetero-plugin-nonparam")

sum4_sigma <- summarize_metric(exp2v, FIG4_METHODS, "n0", "diag_sigma_cor")
sum4_sigma$metric <- "sigma^2~correlation"
sum4_score <- summarize_metric(exp2v, FIG4_METHODS, "n0", "diag_score_cor")
sum4_score$metric <- "score~correlation"

sum4 <- rbind(sum4_sigma, sum4_score)
sum4$method <- factor(sum4$method, levels = FIG4_METHODS)
sum4$metric <- factor(sum4$metric, levels = c("sigma^2~correlation", "score~correlation"))

METRIC_LABELS4 <- c(
  "sigma^2~correlation" = expression(sigma^2~correlation),
  "score~correlation"   = expression(score~correlation)
)

p4 <- ggplot(sum4, aes(x = x, y = mean, color = method, shape = method)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.8) +
  geom_errorbar(aes(ymin = lo, ymax = hi),
                width = diff(range(n0_vals)) * 0.02, linewidth = 0.4) +
  facet_wrap(~metric, nrow = 1, labeller = label_parsed) +
  scale_color_manual(values = METHOD_COLORS[FIG4_METHODS],
                      labels = METHOD_LABELS[FIG4_METHODS]) +
  scale_shape_manual(values = METHOD_SHAPES[FIG4_METHODS],
                      labels = METHOD_LABELS[FIG4_METHODS]) +
  scale_x_continuous(breaks = n0_vals) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(x = "Pilot size n0", y = "Correlation", color = "Plugin", shape = "Plugin") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        strip.text = element_text(size = 10, face = "bold"))

ggsave(file.path(dir_main, "fig4_plugin_quality.pdf"), p4, width = 9, height = 4.5)
cat("  Saved: figures/main/fig4_plugin_quality.pdf\n")

# ═════════════════════════════════════════════════════════════════════════
# Figure 5 — Effect of variance-model misspecification
# ═════════════════════════════════════════════════════════════════════════
cat("Building Figure 5 (misspecification) ...\n")

exp5 <- readRDS(file.path("results", "raw", "exp5", "exp5_misspec_n0.rds"))
exp5v <- exp5[!exp5$fail, ]

FIG5_N0 <- 200  # ⚠ n0_default와 일치하는지 확인
exp5_sub <- exp5v[exp5v$n0 == FIG5_N0, ]

# 기준선: oracle, homo (setting/family와 무관한 상수)
oracle_mean <- mean(exp5_sub$pred_mse_N[exp5_sub$method == "OPT-hetero-oracle"], na.rm = TRUE)
homo_mean   <- mean(exp5_sub$pred_mse_N[exp5_sub$method == "OPT-homo"],          na.rm = TRUE)

setting_map <- data.frame(
  method  = c("param-correct", "param-mild", "param-severe",
              "rf-full",       "rf-mild",    "rf-severe"),
  family  = c("param", "param", "param", "RF", "RF", "RF"),
  setting = c("Correct", "Mild", "Severe", "Correct", "Mild", "Severe"),
  stringsAsFactors = FALSE
)

sum5 <- do.call(rbind, lapply(setting_map$method, function(m) {
  x <- exp5_sub$pred_mse_N[exp5_sub$method == m]
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NULL)
  mn <- mean(x); se <- compute_mcse(x)
  data.frame(method = m, mean = mn, lo = mn - 1.96 * se, hi = mn + 1.96 * se,
             stringsAsFactors = FALSE)
}))
sum5 <- merge(sum5, setting_map, by = "method")
sum5$setting <- factor(sum5$setting, levels = c("Correct", "Mild", "Severe"))
sum5$family  <- factor(sum5$family, levels = c("param", "RF"))

FAMILY_COLORS <- c("param" = "#CC79A7", "RF" = "#7570B3")

p5 <- ggplot(sum5, aes(x = setting, y = mean, color = family, group = family)) +
  geom_hline(yintercept = oracle_mean, linetype = "dashed",
             color = "#D55E00", linewidth = 0.5) +
  geom_hline(yintercept = homo_mean, linetype = "dotted",
             color = "#009E73", linewidth = 0.5) +
  geom_line(linewidth = 0.3, alpha = 0.5,
            position = position_dodge(width = 0.15)) +
  geom_point(size = 3, position = position_dodge(width = 0.15)) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.1,
                linewidth = 0.45, position = position_dodge(width = 0.15)) +
  annotate("text", x = 3.05, y = oracle_mean, label = "Oracle",
           hjust = 0, vjust = -0.5, size = 3, color = "#D55E00") +
  annotate("text", x = 3.05, y = homo_mean, label = "Homo",
           hjust = 0, vjust = 1.3, size = 3, color = "#009E73") +
  scale_color_manual(values = FAMILY_COLORS,
                      labels = c("param" = "Parametric plugin",
                                 "RF"    = "RF (nonparam) plugin")) +
  coord_cartesian(xlim = c(0.8, 3.5), clip = "off") +
  labs(x = "Specification", y = "Predictive excess risk",
       color = "Plugin family",
       title = "Effect of variance-model misspecification") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.margin = margin(5.5, 30, 5.5, 5.5))

ggsave(file.path(dir_main, "fig5_misspecification.pdf"), p5, width = 6.5, height = 4.5)
cat("  Saved: figures/main/fig5_misspecification.pdf\n")

# ═════════════════════════════════════════════════════════════════════════
# Figure 6 — Runtime
#   ⚠ mtry=2/nodesize=2 (재튜닝 전) 데이터 기준 — 잠정본.
#     runtime 재실행 완료되면 이 블록만 그대로 재실행하면 됨.
# ═════════════════════════════════════════════════════════════════════════
cat("Building Figure 6 (runtime) ...\n")

RUNTIME_METHOD_LINES <- c(
  "FULL", "SRS", "LEV-IPW", "OPT-homo", "OPT-hetero-oracle",
  "OPT-hetero-plugin-param", "OPT-hetero-plugin-nonparam"
)
RUNTIME_METHOD_LABELS <- c(
  "FULL"                       = "FULL OLS",
  "SRS"                        = "SRS",
  "LEV-IPW"                    = "LEV-IPW",
  "OPT-homo"                   = "OPT-homo",
  "OPT-hetero-oracle"          = "OPT-oracle",
  "OPT-hetero-plugin-param"    = "OPT-plugin (param)",
  "OPT-hetero-plugin-nonparam" = "OPT-plugin (RF)"
)
RUNTIME_METHOD_COLORS <- c(
  "FULL"                       = "#333333",
  "SRS"                        = "#999999",
  "LEV-IPW"                    = "#0072B2",
  "OPT-homo"                   = "#009E73",
  "OPT-hetero-oracle"          = "#D55E00",
  "OPT-hetero-plugin-param"    = "#CC79A7",
  "OPT-hetero-plugin-nonparam" = "#7570B3"
)
RUNTIME_METHOD_SHAPES <- c(
  "FULL"                       = 3L,
  "SRS"                        = 16L,
  "LEV-IPW"                    = 17L,
  "OPT-homo"                   = 15L,
  "OPT-hetero-oracle"          = 18L,
  "OPT-hetero-plugin-param"    = 8L,
  "OPT-hetero-plugin-nonparam" = 3L
)

# ── 6a. Stage decomposition stacked bar (k=1000 고정) ──────────────────────
decomp_path <- file.path("results", "raw", "exp3", "exp3_runtime_decomposed.rds")

if (!file.exists(decomp_path)) {
  cat("  ⚠ Decomposition 파일 없음 (", decomp_path, ") — 6a 스킵\n")
} else {
  exp3_decomp <- readRDS(decomp_path)
  fig6a_dat <- exp3_decomp[exp3_decomp$k == 1000, ]

  STAGE_ORDER_FIG <- c(
    "pilot_draw", "pilot_fit", "variance_est",
    "moment_matrix", "score_compute", "subsample_draw", "final_fit"
  )
  STAGE_COLORS <- c(
    "pilot_draw"     = "#CC79A7", "pilot_fit"      = "#E69F00",
    "variance_est"   = "#D55E00", "moment_matrix"  = "#0072B2",
    "score_compute"  = "#56B4E9", "subsample_draw" = "#009E73",
    "final_fit"      = "#999999"
  )
  STAGE_LABELS <- c(
    "pilot_draw" = "Pilot draw", "pilot_fit" = "Pilot fit",
    "variance_est" = "Variance est.", "moment_matrix" = "Moment matrix",
    "score_compute" = "Score compute", "subsample_draw" = "Subsample draw",
    "final_fit" = "Final fit"
  )

  fig6a_dat$stage  <- factor(fig6a_dat$stage, levels = STAGE_ORDER_FIG)
  fig6a_dat$method <- factor(fig6a_dat$method, levels = RUNTIME_METHOD_LINES)
  fig6a_dat <- fig6a_dat[!is.na(fig6a_dat$stage) & !is.na(fig6a_dat$method), ]

  p6a <- ggplot(fig6a_dat, aes(x = method, y = mean_time, fill = stage)) +
    geom_bar(stat = "identity", position = "stack", width = 0.6) +
    scale_fill_manual(values = STAGE_COLORS, labels = STAGE_LABELS, drop = FALSE) +
    scale_x_discrete(labels = RUNTIME_METHOD_LABELS[RUNTIME_METHOD_LINES], drop = FALSE) +
    labs(x = "Method", y = "Mean runtime (seconds)", fill = "Stage",
         title = "Runtime decomposition by stage (k = 1000, N = 100,000)",
         subtitle = "Provisional: pre-retuning RF hyperparameters (mtry=2, nodesize=2)") +
    theme_bw(base_size = 12) +
    theme(legend.position = "right", panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 25, hjust = 1),
          plot.subtitle = element_text(size = 8, color = "gray40"))

  ggsave(file.path(dir_main, "fig6a_runtime_stacked.pdf"), p6a, width = 9, height = 5)
  cat("  Saved: figures/main/fig6a_runtime_stacked.pdf\n")
}

# ── 6b. N scaling line plot (log-log) ──────────────────────────────────────
scale_path <- file.path("results", "raw", "exp3", "exp3_runtime_scaling.rds")

if (!file.exists(scale_path)) {
  cat("  ⚠ Scaling 파일 없음 (", scale_path, ") — 6b 스킵\n")
} else {
  exp3_scale <- readRDS(scale_path)

  scale_total_rows <- list()
  for (m in RUNTIME_METHOD_LINES) {
    for (nv in sort(unique(exp3_scale$N))) {
      sub <- exp3_scale[exp3_scale$method == m & exp3_scale$N == nv, ]
      if (nrow(sub) == 0L) next
      scale_total_rows[[length(scale_total_rows) + 1L]] <- data.frame(
        method = m, N = nv, total_mean = sum(sub$mean_time),
        stringsAsFactors = FALSE
      )
    }
  }
  scale_total_df <- do.call(rbind, scale_total_rows)
  scale_total_df$method <- factor(scale_total_df$method, levels = RUNTIME_METHOD_LINES)

  n_breaks <- sort(unique(scale_total_df$N))
  n_labels <- format(n_breaks, big.mark = ",", scientific = FALSE, trim = TRUE)

  p6b <- ggplot(scale_total_df, aes(x = N, y = total_mean, color = method, shape = method)) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 2.5) +
    scale_color_manual(values = RUNTIME_METHOD_COLORS[RUNTIME_METHOD_LINES],
                        labels = RUNTIME_METHOD_LABELS[RUNTIME_METHOD_LINES]) +
    scale_shape_manual(values = RUNTIME_METHOD_SHAPES[RUNTIME_METHOD_LINES],
                        labels = RUNTIME_METHOD_LABELS[RUNTIME_METHOD_LINES]) +
    scale_x_log10(breaks = n_breaks, labels = n_labels) +
    scale_y_log10() +
    labs(x = "Full-data size N (log scale)", y = "Mean total runtime, seconds (log scale)",
         color = "Method", shape = "Method",
         title = "Runtime scaling with N",
         subtitle = "Provisional: pre-retuning RF hyperparameters (mtry=2, nodesize=2)") +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 20, hjust = 1),
          plot.subtitle = element_text(size = 8, color = "gray40"))

  ggsave(file.path(dir_main, "fig6b_runtime_scaling.pdf"), p6b, width = 7.5, height = 4.8)
  cat("  Saved: figures/main/fig6b_runtime_scaling.pdf\n")
}

cat("\nAll figures done.\n")