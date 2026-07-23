# 05_mcse_summary.R
# Monte Carlo 표준오차, 요약 테이블 출력, 페어 차이 분석

# Monte Carlo 표준오차: sd(x) / sqrt(n)
compute_mcse <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(NA_real_)
  sd(x) / sqrt(length(x))
}

# method × k 기준 요약 (mean + MCSE)
summarize_metric <- function(df, metric) {
  df_v <- df[!df$fail & is.finite(df[[metric]]), ]
  rows <- list()
  for (m in unique(df_v$method)) {
    sub_m <- df_v[df_v$method == m, ]
    for (k_ in sort(unique(sub_m$k))) {
      vals <- sub_m[[metric]][sub_m$k == k_]
      rows[[length(rows) + 1L]] <- data.frame(
        method = m, k = k_,
        mean   = mean(vals),
        mcse   = compute_mcse(vals),
        n_ok   = length(vals),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

# 페어 차이: rep 단위 (m1 metric - m2 metric), MCSE 포함
paired_diff_summary <- function(df, m1, m2, metric = "pred_mse_N") {
  cols <- c("rep", "k", metric)
  d1   <- df[df$method == m1 & !df$fail, cols]
  d2   <- df[df$method == m2 & !df$fail, cols]
  mg   <- merge(d1, d2, by = c("rep", "k"), suffixes = c(".m1", ".m2"))
  mg$diff <- mg[[paste0(metric, ".m1")]] - mg[[paste0(metric, ".m2")]]

  rows <- list()
  for (k_ in sort(unique(mg$k))) {
    d <- mg$diff[mg$k == k_]
    d <- d[is.finite(d)]
    md <- mean(d)
    se <- compute_mcse(d)
    rows[[length(rows) + 1L]] <- data.frame(
      m1 = m1, m2 = m2, k = k_, metric = metric,
      mean_diff = md, mcse_diff = se,
      n_pairs   = length(d),
      lower     = md - 1.96 * se,
      upper     = md + 1.96 * se,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

# metric 테이블 콘솔 출력 (mean ± mcse)
print_metric_table <- function(df, metric, ks, title,
                               method_order = METHOD_ORDER) {
  agg <- summarize_metric(df[df$k %in% ks, ], metric)

  cat(sprintf("\n%s\n", strrep("=", 65)))
  cat(sprintf(" %s\n", title))
  cat(sprintf("%s\n", strrep("=", 65)))
  cat(sprintf("%-25s", "method \\ k"))
  for (k_ in ks) cat(sprintf("%16d", k_))
  cat("\n")
  cat(strrep("-", 25 + 16 * length(ks)), "\n")

  for (m in method_order) {
    if (!m %in% agg$method) next
    cat(sprintf("%-25s", m))
    for (k_ in ks) {
      row <- agg[agg$method == m & agg$k == k_, ]
      if (nrow(row) == 0L) cat(sprintf("%16s", "NA"))
      else cat(sprintf("  %8.5f±%4.2e", row$mean, row$mcse))
    }
    cat("\n")
  }
}

# 페어 차이 (m1 - m2) 콘솔 출력, 기본 쌍: SRS-plugin, homo-plugin, plugin-oracle
print_paired_diffs <- function(df, ks, metric = "pred_mse_N",
                               pairs = list(
                                 c("SRS",               "OPT-hetero-plugin-param"),
                                 c("OPT-homo",          "OPT-hetero-plugin-param"),
                                 c("OPT-hetero-plugin-param", "OPT-hetero-oracle"),
                                 c('OPT-hetero-plugin-param', 'OPT-hetero-plugin-nonparam')
                               )) {
  cat(sprintf("\n%s\n", strrep("=", 78)))
  cat(sprintf(" Paired differences [m1 - m2]  metric: %s\n", metric))
  cat(sprintf("%s\n", strrep("=", 78)))
  cat(sprintf("%-25s %-25s %6s %12s %12s %12s %12s\n",
              "m1", "m2", "k", "mean_diff", "mcse_diff", "lower", "upper"))
  cat(strrep("-", 108), "\n")

  for (pr in pairs) {
    pd <- tryCatch(
      paired_diff_summary(df[df$k %in% ks, ], pr[1], pr[2], metric),
      error = function(e) NULL
    )
    if (is.null(pd) || nrow(pd) == 0L) next
    for (i in seq_len(nrow(pd))) {
      cat(sprintf("%-25s %-25s %6d %12.6f %12.6f %12.6f %12.6f\n",
                  pd$m1[i], pd$m2[i], pd$k[i],
                  pd$mean_diff[i], pd$mcse_diff[i],
                  pd$lower[i], pd$upper[i]))
    }
  }
}

# clip / realized sample size 콘솔 출력
print_clip_table <- function(df, ks, method_order = METHOD_ORDER) {
  cat(sprintf("\n%s\n", strrep("=", 65)))
  cat(" Clip 및 실현 샘플 크기 (평균)\n")
  cat(sprintf("%s\n", strrep("=", 65)))
  cat(sprintf("%-25s %6s %10s %10s\n", "method", "k", "n_clip", "n_realized"))
  cat(strrep("-", 55), "\n")
  df_v <- df[!df$fail, ]
  for (m in method_order) {
    for (k_ in ks) {
      sub <- df_v[df_v$method == m & df_v$k == k_, ]
      if (nrow(sub) == 0L) next
      cat(sprintf("%-25s %6d %10.2f %10.2f\n",
                  m, k_, mean(sub$n_clip), mean(sub$n_realized)))
    }
  }
}

# IPW 불안정성 진단 출력 (ESS_pi, max_weight, cond_XtWX)
print_instability_table <- function(df, ks, method_order = METHOD_ORDER) {
  cat(sprintf("\n%s\n", strrep("=", 78)))
  cat(" IPW 불안정성 진단 (ESS_pi, max_weight, cond_XtWX)\n")
  cat(sprintf("%s\n", strrep("=", 78)))
  cat(sprintf("%-25s %6s %10s %12s %14s\n",
              "method", "k", "ESS_pi", "max_weight", "cond_XtWX"))
  cat(strrep("-", 72), "\n")
  df_v <- df[!df$fail, ]
  for (m in method_order) {
    for (k_ in ks) {
      sub <- df_v[df_v$method == m & df_v$k == k_, ]
      if (nrow(sub) == 0L) next
      cat(sprintf("%-25s %6d %10.1f %12.2f %14.2f\n",
                  m, k_,
                  mean(sub$ess_pi,     na.rm = TRUE),
                  mean(sub$max_weight, na.rm = TRUE),
                  mean(sub$cond_XtWX,  na.rm = TRUE)))
    }
  }
}
