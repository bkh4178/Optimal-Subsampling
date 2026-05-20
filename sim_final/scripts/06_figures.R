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

PLUGIN_SETTINGS <- c("correct_scale", "variance_linear", "mild_omitted", "severe_omitted")
PLUGIN_METHODS  <- c("plugin-correct", "plugin-varlinear", "plugin-mild", "plugin-severe")

SETTING_LABELS <- c(
  "correct_scale"   = "Correct\n(sigma-scale)",
  "variance_linear" = "Var-linear",
  "mild_omitted"    = "Mild\nomitted",
  "severe_omitted"  = "Severe\nomitted"
)

SETTING_COLORS <- c(
  "correct_scale"   = "#009E73",
  "variance_linear" = "#0072B2",
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
    labels = c("Correct (sigma-scale)", "Var-linear", "Mild omitted", "Severe omitted")
  ) +
  labs(
    x     = "Mean score correlation (plugin vs oracle)",
    y     = "Excess Risk: plugin - oracle (mean diff)",
    color = "Misspecification"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position  = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave(file.path(dir_main, "fig3_exp4_scatter.pdf"),
       p3, width = 6, height = 5)
cat("  Saved: figures/main/fig3_exp4_scatter.pdf\n")

cat("\nAll figures done.\n")
