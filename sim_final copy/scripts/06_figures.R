# 06_figures.R
# Generate figures from simulation results
# Run from: sim_final/
# Requires: ggplot2

library(ggplot2)

# ── Output directories ────────────────────────────────────────────────────────
dir_main     <- file.path("figures", "main")
dir_appendix <- file.path("figures", "appendix")
for (d in c(dir_main, dir_appendix)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# ── Helpers ───────────────────────────────────────────────────────────────────
compute_mcse <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))

# ── Method aesthetics (consistent across all figures) ─────────────────────────
METHOD_LINES <- c("SRS", "LEV-IPW", "OPT-homo", "OPT-hetero-oracle", "OPT-hetero-plugin")

METHOD_COLORS <- c(
  "SRS"               = "#999999",
  "LEV-IPW"           = "#0072B2",
  "OPT-homo"          = "#009E73",
  "OPT-hetero-oracle" = "#D55E00",
  "OPT-hetero-plugin" = "#CC79A7"
)

METHOD_SHAPES <- c(
  "SRS"               = 16L,
  "LEV-IPW"           = 17L,
  "OPT-homo"          = 15L,
  "OPT-hetero-oracle" = 18L,
  "OPT-hetero-plugin" = 8L
)

METHOD_LABELS <- c(
  "SRS"               = "SRS",
  "LEV-IPW"           = "LEV-IPW",
  "OPT-homo"          = "OPT-homo",
  "OPT-hetero-oracle" = "OPT-oracle",
  "OPT-hetero-plugin" = "OPT-plugin"
)

# ─────────────────────────────────────────────────────────────────────────────
# Figure 1 — Exp1: k vs Excess Risk
# ─────────────────────────────────────────────────────────────────────────────
cat("Building Figure 1 ...\n")

exp1 <- readRDS(file.path("results", "raw", "exp1", "exp1_k_comparison.rds"))
exp1 <- exp1[!exp1$fail, ]

# Summarise per (method, k)
sum1 <- do.call(rbind, lapply(METHOD_LINES, function(m) {
  sub <- exp1[exp1$method == m, ]
  do.call(rbind, lapply(sort(unique(sub$k)), function(kv) {
    x  <- sub$excess_risk[sub$k == kv]
    mn <- mean(x, na.rm = TRUE)
    se <- compute_mcse(x)
    data.frame(method = m, k = kv,
               mean   = mn,
               lo     = mn - 1.96 * se,
               hi     = mn + 1.96 * se,
               stringsAsFactors = FALSE)
  }))
}))
sum1$method <- factor(sum1$method, levels = METHOD_LINES)

# FULL reference (averaged across all k, since FULL does not depend on k)
full_mean <- mean(exp1$excess_risk[exp1$method == "FULL"], na.rm = TRUE)

k_vals <- sort(unique(sum1$k))

p1 <- ggplot(sum1, aes(x = k, y = mean, color = method, shape = method)) +
  geom_hline(yintercept = full_mean,
             linetype = "dashed", color = "black", linewidth = 0.5) +
  annotate("text",
           x = min(k_vals), y = full_mean,
           label = "FULL OLS", hjust = 0, vjust = -0.6, size = 3) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = lo, ymax = hi),
                width = diff(range(k_vals)) * 0.02,
                linewidth = 0.45) +
  scale_color_manual(
    values = METHOD_COLORS[METHOD_LINES],
    labels = METHOD_LABELS[METHOD_LINES]
  ) +
  scale_shape_manual(
    values = METHOD_SHAPES[METHOD_LINES],
    labels = METHOD_LABELS[METHOD_LINES]
  ) +
  scale_x_continuous(breaks = k_vals) +
  labs(
    x     = "Subsample size k",
    y     = "Excess Risk (mean ± 1.96 MCSE)",
    color = "Method",
    shape = "Method"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position  = "bottom",
    legend.title     = element_text(size = 10),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(dir_main, "fig1_exp1_k_excess_risk.pdf"),
       p1, width = 7, height = 4.5)
cat("  Saved: figures/main/fig1_exp1_k_excess_risk.pdf\n")

# ─────────────────────────────────────────────────────────────────────────────
# Figure 2 — Exp2: n0 vs sigma_cor / score_cor
# ─────────────────────────────────────────────────────────────────────────────
cat("Building Figure 2 ...\n")

exp2    <- readRDS(file.path("results", "raw", "exp2", "exp2_n0_sensitivity.rds"))
plugin2 <- exp2[exp2$method == "OPT-hetero-plugin" & !exp2$fail, ]

sum2 <- do.call(rbind, lapply(sort(unique(plugin2$n0)), function(n0v) {
  sub <- plugin2[plugin2$n0 == n0v, ]
  data.frame(
    n0        = n0v,
    sigma_cor = mean(sub$diag_sigma_cor, na.rm = TRUE),
    score_cor = mean(sub$diag_score_cor, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))

# Reshape to long
sum2_long <- rbind(
  data.frame(n0     = sum2$n0,
             value  = sum2$sigma_cor,
             metric = "sigma^2~correlation",
             stringsAsFactors = FALSE),
  data.frame(n0     = sum2$n0,
             value  = sum2$score_cor,
             metric = "score~correlation",
             stringsAsFactors = FALSE)
)
sum2_long$metric <- factor(sum2_long$metric,
                            levels = c("sigma^2~correlation", "score~correlation"))

COR_COLORS   <- c("sigma^2~correlation" = "#D55E00", "score~correlation" = "#0072B2")
COR_LTYPES   <- c("sigma^2~correlation" = "solid",   "score~correlation" = "dashed")
COR_SHAPES   <- c("sigma^2~correlation" = 16L,       "score~correlation" = 17L)
COR_LABELS   <- c("sigma^2~correlation" = expression(sigma^2~correlation),
                  "score~correlation"   = expression(score~correlation))

p2 <- ggplot(sum2_long, aes(x = n0, y = value,
                             color = metric, linetype = metric, shape = metric)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 3) +
  scale_color_manual(values = COR_COLORS,  labels = COR_LABELS) +
  scale_linetype_manual(values = COR_LTYPES, labels = COR_LABELS) +
  scale_shape_manual(values = COR_SHAPES,  labels = COR_LABELS) +
  scale_x_continuous(breaks = sort(unique(sum2$n0))) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    x        = "Pilot size n0",
    y        = "Correlation",
    color    = NULL,
    linetype = NULL,
    shape    = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position  = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave(file.path(dir_main, "fig2_exp2_n0_cor.pdf"),
       p2, width = 6, height = 4)
cat("  Saved: figures/main/fig2_exp2_n0_cor.pdf\n")

# ─────────────────────────────────────────────────────────────────────────────
# Figure 3 — Exp4: score_cor vs plugin–oracle gap
# ─────────────────────────────────────────────────────────────────────────────
cat("Building Figure 3 ...\n")

exp4 <- readRDS(file.path("results", "raw", "exp4", "exp4_misspecification.rds"))

PLUGIN_SETTINGS <- c("correct_scale", "mild_omitted", "severe_omitted")
PLUGIN_METHODS  <- c("plugin-correct", "plugin-mild", "plugin-severe")

SETTING_LABELS <- c(
  "correct_scale"   = "Correct\n(sigma-scale)",
  "mild_omitted"    = "Mild\nomitted",
  "severe_omitted"  = "Severe\nomitted"
)

SETTING_COLORS <- c(
  "correct_scale"   = "#009E73",
  "mild_omitted"    = "#E69F00",
  "severe_omitted"  = "#D55E00"
)

valid4 <- exp4[!exp4$fail, ]

# Mean score_cor per plugin setting
sc_df <- do.call(rbind, lapply(seq_along(PLUGIN_SETTINGS), function(i) {
  sub <- valid4[valid4$method == PLUGIN_METHODS[i], ]
  data.frame(
    misspec_setting = PLUGIN_SETTINGS[i],
    mean_score_cor  = mean(sub$score_cor, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))

# Paired diff: plugin - oracle (merge by rep × k)
pd_df <- do.call(rbind, lapply(seq_along(PLUGIN_SETTINGS), function(i) {
  pm <- PLUGIN_METHODS[i]
  d1 <- exp4[exp4$method == pm,                  c("rep", "k", "excess_risk")]
  d2 <- exp4[exp4$method == "OPT-hetero-oracle",  c("rep", "k", "excess_risk")]
  mg <- merge(d1, d2, by = c("rep", "k"), suffixes = c(".plugin", ".oracle"))
  d  <- mg$excess_risk.plugin - mg$excess_risk.oracle
  d  <- d[is.finite(d)]
  md <- mean(d)
  se <- compute_mcse(d)
  data.frame(
    misspec_setting = PLUGIN_SETTINGS[i],
    mean_diff       = md,
    lower           = md - 1.96 * se,
    upper           = md + 1.96 * se,
    stringsAsFactors = FALSE
  )
}))

sum3 <- merge(sc_df, pd_df, by = "misspec_setting")
sum3 <- sum3[complete.cases(sum3), ]
sum3$misspec_setting <- factor(sum3$misspec_setting, levels = PLUGIN_SETTINGS)
sum3$label <- SETTING_LABELS[as.character(sum3$misspec_setting)]

# Nudge labels upward by a fraction of the y-range
y_range  <- range(c(sum3$lower, sum3$upper))
nudge_up <- diff(y_range) * 0.08

p3 <- ggplot(sum3, aes(x = mean_score_cor, y = mean_diff, color = misspec_setting)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = 0.008, linewidth = 0.6, show.legend = FALSE) +
  geom_point(size = 4) +
  geom_text(aes(label = label),
            vjust = 0, nudge_y = nudge_up,
            size = 3, lineheight = 0.85, show.legend = FALSE) +
  scale_color_manual(
    values = SETTING_COLORS,
    labels = c("Correct\n(sigma-scale)", "Mild omitted", "Severe omitted")
  ) +
  labs(
    x     = "Mean score correlation (plugin vs oracle)",
    y     = "Excess Risk: plugin - oracle (mean diff)",
    color = "Misspecification"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position  = "bottom",
    legend.text      = element_text(size = 8),
    panel.grid.minor = element_blank(),
    plot.margin      = margin(5.5, 20, 5.5, 5.5)
  )

ggsave(file.path(dir_main, "fig3_exp4_scatter.pdf"),
       p3, width = 8, height = 5)
cat("  Saved: figures/main/fig3_exp4_scatter.pdf\n")

# ─────────────────────────────────────────────────────────────────────────────
# Figure 4 — Exp3-a: Stage decomposition stacked bar (k=1000, N=100,000)
# ─────────────────────────────────────────────────────────────────────────────
cat("Building Figure 4 ...\n")

exp3_decomp <- readRDS(file.path("results", "raw", "exp3",
                                 "exp3_runtime_decomposed.rds"))
fig4_dat <- exp3_decomp[exp3_decomp$k == 1000, ]

RUNTIME_METHOD_LINES <- c(
  "FULL", "SRS", "LEV-IPW", "OPT-homo",
  "OPT-hetero-oracle", "OPT-hetero-plugin"
)
RUNTIME_METHOD_LABELS <- c(
  "FULL"              = "FULL OLS",
  "SRS"               = "SRS",
  "LEV-IPW"           = "LEV-IPW",
  "OPT-homo"          = "OPT-homo",
  "OPT-hetero-oracle" = "OPT-oracle",
  "OPT-hetero-plugin" = "OPT-plugin"
)

STAGE_ORDER_FIG <- c(
  "pilot_draw", "pilot_fit", "variance_est",
  "moment_matrix", "score_compute", "subsample_draw", "final_fit"
)
STAGE_COLORS <- c(
  "pilot_draw"     = "#CC79A7",
  "pilot_fit"      = "#E69F00",
  "variance_est"   = "#D55E00",
  "moment_matrix"  = "#0072B2",
  "score_compute"  = "#56B4E9",
  "subsample_draw" = "#009E73",
  "final_fit"      = "#999999"
)
STAGE_LABELS <- c(
  "pilot_draw"     = "Pilot draw",
  "pilot_fit"      = "Pilot fit",
  "variance_est"   = "Variance est.",
  "moment_matrix"  = "Moment matrix",
  "score_compute"  = "Score compute",
  "subsample_draw" = "Subsample draw",
  "final_fit"      = "Final fit"
)

fig4_dat$stage <- factor(fig4_dat$stage, levels = STAGE_ORDER_FIG)
fig4_dat$method <- factor(fig4_dat$method, levels = RUNTIME_METHOD_LINES)
fig4_dat <- fig4_dat[!is.na(fig4_dat$stage) & !is.na(fig4_dat$method), ]

p4 <- ggplot(fig4_dat,
             aes(x = method, y = mean_time, fill = stage)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6) +
  scale_fill_manual(
    values = STAGE_COLORS,
    labels = STAGE_LABELS,
    drop   = FALSE
  ) +
  scale_x_discrete(
    labels = RUNTIME_METHOD_LABELS[RUNTIME_METHOD_LINES],
    drop   = FALSE
  ) +
  labs(
    x    = "Method",
    y    = "Mean runtime (seconds)",
    fill = "Stage",
    title = "Runtime decomposition by stage (k = 1000, N = 100,000)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position  = "right",
    panel.grid.minor = element_blank(),
    axis.text.x      = element_text(angle = 25, hjust = 1)
  )

ggsave(file.path(dir_main, "fig4_runtime_stacked.pdf"),
       p4, width = 8, height = 5)
cat("  Saved: figures/main/fig4_runtime_stacked.pdf\n")

# ─────────────────────────────────────────────────────────────────────────────
# Figure 5 — Exp3-b: N scaling line plot
# ─────────────────────────────────────────────────────────────────────────────
cat("Building Figure 5 ...\n")

exp3_scale <- readRDS(file.path("results", "raw", "exp3",
                                "exp3_runtime_scaling.rds"))

RUNTIME_METHOD_COLORS <- c(
  "FULL"              = "#333333",
  "SRS"               = "#999999",
  "LEV-IPW"           = "#0072B2",
  "OPT-homo"          = "#009E73",
  "OPT-hetero-oracle" = "#D55E00",
  "OPT-hetero-plugin" = "#CC79A7"
)
RUNTIME_METHOD_SHAPES <- c(
  "FULL"              = 3L,
  "SRS"               = 16L,
  "LEV-IPW"           = 17L,
  "OPT-homo"          = 15L,
  "OPT-hetero-oracle" = 18L,
  "OPT-hetero-plugin" = 8L
)

scale_total_rows <- list()
for (m in RUNTIME_METHOD_LINES) {
  for (nv in sort(unique(exp3_scale$N))) {
    sub <- exp3_scale[exp3_scale$method == m & exp3_scale$N == nv, ]
    if (nrow(sub) == 0L) next
    scale_total_rows[[length(scale_total_rows) + 1L]] <- data.frame(
      method     = m,
      N          = nv,
      total_mean = sum(sub$mean_time),
      stringsAsFactors = FALSE
    )
  }
}
scale_total_df <- do.call(rbind, scale_total_rows)
scale_total_df$method <- factor(
  scale_total_df$method, levels = RUNTIME_METHOD_LINES
)

n_breaks <- sort(unique(scale_total_df$N))
n_labels <- format(n_breaks, big.mark = ",", scientific = FALSE, trim = TRUE)

p5 <- ggplot(scale_total_df,
             aes(x = N, y = total_mean,
                 color = method, shape = method)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.5) +
  scale_color_manual(
    values = RUNTIME_METHOD_COLORS[RUNTIME_METHOD_LINES],
    labels = RUNTIME_METHOD_LABELS[RUNTIME_METHOD_LINES]
  ) +
  scale_shape_manual(
    values = RUNTIME_METHOD_SHAPES[RUNTIME_METHOD_LINES],
    labels = RUNTIME_METHOD_LABELS[RUNTIME_METHOD_LINES]
  ) +
  scale_x_log10(breaks = n_breaks, labels = n_labels) +
  scale_y_log10() +
  labs(
    x     = "Full-data size N (log scale)",
    y     = "Mean total runtime, seconds (log scale)",
    color = "Method",
    shape = "Method"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position  = "bottom",
    panel.grid.minor = element_blank(),
    axis.text.x      = element_text(angle = 20, hjust = 1)
  )

ggsave(file.path(dir_main, "fig5_runtime_scaling.pdf"),
       p5, width = 7, height = 4.5)
cat("  Saved: figures/main/fig5_runtime_scaling.pdf\n")

cat("\nAll figures done.\n")
