# 04_metrics.R
# run_one_rep  : 단일 replicate 실행 – 전 method 가 공통 pilot 을 공유
# results_to_rows : 결과 list → data.frame 변환

run_one_rep <- function(dat, k, X_test, mu_test, y_test, n0 = n0_default) {
  N    <- nrow(dat$X)
  beta <- dat$beta

  # ── Pilot (모든 method 공통) ─────────────────────────────────────────────────
  idx_pilot  <- sample.int(N, n0, replace = FALSE)
  pi_pilot   <- rep(n0 / N, N)
  beta_pilot <- fit_ipw(dat$X, dat$y, idx_pilot, pi = pi_pilot)
  if (is.null(beta_pilot)) beta_pilot <- beta

  # ── Plugin sigma^2 추정 ───────────────────────────────────────────────────────
  sigma2_hat     <- estimate_sigma2_plugin(dat, idx_pilot, beta_pilot)
  oracle_score   <- sqrt(dat$sigma2_vec * dat$ell)
  plugin_score   <- sqrt(sigma2_hat * dat$ell)
  diag_sigma_cor <- safe_cor(sigma2_hat, dat$sigma2_vec)
  diag_score_cor <- safe_cor(plugin_score, oracle_score)

  # ── Sampling scores ──────────────────────────────────────────────────────────
  scores <- get_scores(dat, sigma2_hat)

  # ── 단일 method 실행 (내부 함수) ─────────────────────────────────────────────
  run_method <- function(method, score) {
    samp      <- sample_bernoulli(score, k)
    idx_final <- samp$idx
    pi_raw    <- samp$pi
    n_clip    <- sum(pi_raw >= 1 - 1e-12)

    nz         <- pi_raw[pi_raw > 0]
    ess_pi     <- if (length(nz) > 0L) sum(pi_raw)^2 / sum(pi_raw^2) else NA_real_
    max_weight <- if (length(nz) > 0L) max(1 / pmax(nz, 1e-10)) else NA_real_

    beta_hat <- fit_ipw(dat$X, dat$y, idx_final, pi = pi_raw)

    cond_XtWX <- tryCatch({
      if (length(idx_final) < ncol(dat$X)) return(NA_real_)
      w    <- 1 / pmax(pi_raw[idx_final], 1e-10)
      XtWX <- crossprod(dat$X[idx_final, , drop = FALSE] * sqrt(w))
      ev   <- eigen(XtWX, symmetric = TRUE, only.values = TRUE)$values
      ev   <- ev[is.finite(ev) & ev > 0]
      if (length(ev) == 0L) return(NA_real_)
      max(ev) / max(min(ev), 1e-12)
    }, error = function(e) NA_real_)

    is_plugin <- (method == "OPT-hetero-plugin")

    if (is.null(beta_hat)) {
      return(list(
        method = method, fail = TRUE,
        excess_risk = NA_real_, pred_mse = NA_real_, param_mse = NA_real_,
        c1_pi = NA_real_,
        n_realized = length(idx_final), n_clip = n_clip,
        ess_pi = ess_pi, max_weight = max_weight, cond_XtWX = cond_XtWX,
        diag_sigma_cor = if (is_plugin) diag_sigma_cor else NA_real_,
        diag_score_cor = if (is_plugin) diag_score_cor else NA_real_
      ))
    }

    list(
      method      = method, fail = FALSE,
      excess_risk = as.numeric(t(beta_hat - beta) %*% dat$A_hat %*% (beta_hat - beta)),
      pred_mse    = mean(as.numeric(mu_test - X_test %*% beta_hat)^2),
      param_mse   = mean((beta_hat - beta)^2),
      c1_pi       = mean(dat$sigma2_vec * dat$ell / pmax(pi_raw, 1e-10)),
      n_realized  = length(idx_final), n_clip = n_clip,
      ess_pi = ess_pi, max_weight = max_weight, cond_XtWX = cond_XtWX,
      diag_sigma_cor = if (is_plugin) diag_sigma_cor else NA_real_,
      diag_score_cor = if (is_plugin) diag_score_cor else NA_real_
    )
  }

  # ── FULL: 전체 데이터 OLS (pilot / subsampling 없음) ─────────────────────────
  full_samp <- sample_full(dat)
  bhat_full <- full_samp$beta_hat
  full_result <- if (is.null(bhat_full)) {
    list(
      method = "FULL", fail = TRUE,
      excess_risk = NA_real_, pred_mse = NA_real_, param_mse = NA_real_,
      c1_pi = NA_real_,
      n_realized = full_samp$n_realized, n_clip = full_samp$n_clip,
      ess_pi = NA_real_, max_weight = NA_real_, cond_XtWX = NA_real_,
      diag_sigma_cor = NA_real_, diag_score_cor = NA_real_
    )
  } else {
    list(
      method      = "FULL", fail = FALSE,
      excess_risk = as.numeric(t(bhat_full - beta) %*% dat$A_hat %*% (bhat_full - beta)),
      pred_mse    = mean(as.numeric(mu_test - X_test %*% bhat_full)^2),
      param_mse   = mean((bhat_full - beta)^2),
      c1_pi       = NA_real_,
      n_realized  = full_samp$n_realized, n_clip = full_samp$n_clip,
      ess_pi = NA_real_, max_weight = NA_real_, cond_XtWX = NA_real_,
      diag_sigma_cor = NA_real_, diag_score_cor = NA_real_
    )
  }

  results <- list()
  results[["FULL"]]              <- full_result
  results[["SRS"]]               <- run_method("SRS",               scores[["SRS"]])
  results[["LEV-IPW"]]           <- run_method("LEV-IPW",           scores[["LEV-IPW"]])
  results[["OPT-homo"]]          <- run_method("OPT-homo",          scores[["OPT-homo"]])
  results[["OPT-hetero-oracle"]] <- run_method("OPT-hetero-oracle", scores[["OPT-hetero-oracle"]])
  results[["OPT-hetero-plugin"]] <- run_method("OPT-hetero-plugin", scores[["OPT-hetero-plugin"]])
  results
}

results_to_rows <- function(rep_results, rep_id, k) {
  do.call(rbind, lapply(names(rep_results), function(m) {
    rv <- rep_results[[m]]
    data.frame(
      k              = k,
      rep            = rep_id,
      method         = rv$method,
      fail           = rv$fail,
      excess_risk    = rv$excess_risk,
      pred_mse       = rv$pred_mse,
      param_mse      = rv$param_mse,
      c1_pi          = rv$c1_pi,
      n_realized     = rv$n_realized,
      n_clip         = rv$n_clip,
      ess_pi         = rv$ess_pi,
      max_weight     = rv$max_weight,
      cond_XtWX      = rv$cond_XtWX,
      diag_sigma_cor = rv$diag_sigma_cor,
      diag_score_cor = rv$diag_score_cor,
      stringsAsFactors = FALSE
    )
  }))
}
