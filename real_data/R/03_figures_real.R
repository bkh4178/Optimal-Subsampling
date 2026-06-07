# =============================================================================
# 03_figures_real.R
# NYC Yellow Taxi — Pseudo Excess Risk figure
# 실행 위치: real_data/R/
# =============================================================================

res <- readRDS("../results/real_results.rds")

# 시뮬레이션 figure 스타일 참고
# 색상 및 선형 설정
method_colors <- c(
  "UNI-IPW"           = "#999999",
  "LEV-IPW"           = "#E69F00",
  "OPT-homo"          = "#56B4E9",
  "OPT-hetero-plugin" = "#D55E00"
)
method_lty <- c(
  "UNI-IPW"           = 2,
  "LEV-IPW"           = 4,
  "OPT-homo"          = 3,
  "OPT-hetero-plugin" = 1
)
method_pch <- c(
  "UNI-IPW"           = 1,
  "LEV-IPW"           = 2,
  "OPT-homo"          = 0,
  "OPT-hetero-plugin" = 16
)

methods_order <- c("UNI-IPW", "LEV-IPW", "OPT-homo", "OPT-hetero-plugin")
k_vals <- sort(unique(res$k[res$method != "FULL"]))

# 요약 계산
summary_list <- list()
for (met in methods_order) {
  for (k in k_vals) {
    v <- res$pseudo_er[res$method == met & res$k == k & !res$fail]
    summary_list[[paste(met, k)]] <- data.frame(
      method = met,
      k      = k,
      mean   = mean(v),
      mcse   = sd(v) / sqrt(length(v))
    )
  }
}
df_sum <- do.call(rbind, summary_list)
rownames(df_sum) <- NULL

# =============================================================================
# Figure: k별 Pseudo Excess Risk (log scale)
# =============================================================================
if (!dir.exists("../figures")) dir.create("../figures", recursive = TRUE)

pdf("figures/fig_real_excess_risk.pdf", width = 6.5, height = 5)

par(mar = c(4.5, 4.5, 2, 1.5))

y_all <- c(df_sum$mean - 1.96 * df_sum$mcse,
           df_sum$mean + 1.96 * df_sum$mcse)
ylim  <- range(y_all[y_all > 0], na.rm = TRUE) * c(0.8, 1.3)

plot(NULL,
     xlim = range(k_vals),
     ylim = ylim,
     log  = "y",
     xlab = "Subsample size k",
     ylab = "Pseudo Excess Risk",
     main = "Real Data: Pseudo Excess Risk by Subsample Size",
     xaxt = "n",
     las  = 1,
     cex.axis = 0.85,
     cex.lab  = 0.95)

axis(1, at = k_vals, labels = format(k_vals, big.mark = ","), cex.axis = 0.85)
grid(nx = NA, ny = NULL, col = "gray85", lty = 1)

for (met in methods_order) {
  df_m <- df_sum[df_sum$method == met, ]
  df_m <- df_m[order(df_m$k), ]

  # error bar
  arrows(df_m$k,
         df_m$mean - 1.96 * df_m$mcse,
         df_m$k,
         df_m$mean + 1.96 * df_m$mcse,
         angle = 90, code = 3, length = 0.04,
         col = method_colors[met], lwd = 1.2)

  lines(df_m$k, df_m$mean,
        col = method_colors[met],
        lty = method_lty[met],
        lwd = 1.8)

  points(df_m$k, df_m$mean,
         col = method_colors[met],
         pch = method_pch[met],
         cex = 1.2)
}

legend("topright",
       legend = c("UNI-IPW", "LEV-IPW", "OPT-homo", "OPT-hetero-plugin"),
       col    = method_colors[methods_order],
       lty    = method_lty[methods_order],
       pch    = method_pch[methods_order],
       lwd    = 1.8,
       cex    = 0.82,
       bg     = "white",
       box.lty = 1)

dev.off()
cat("✓ 저장 완료: ../figures/fig_real_excess_risk.pdf\n")
