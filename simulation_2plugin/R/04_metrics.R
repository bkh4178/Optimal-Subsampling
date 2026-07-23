# 04_metrics.R
# run_one_rep  : 단일 replicate 실행 – 전 method 가 공통 pilot 을 공유
# results_to_rows : 결과 list → data.frame 변환

run_one_rep <- function(dat, k, X_test, mu_test, y_test, n0 = n0_default) {
  N    <- nrow(dat$X)
  beta <- dat$beta

  # Metric 2 baseline: MSE0 = 실현된 test set의 irreducible noise
  # (mu_test, y_test 둘 다 rep/method 무관하게 고정이므로 상수)
  mse0 <- mean(as.numeric(y_test - mu_test)^2)

  # ── Pilot (모든 method 공통) ─────────────────────────────────────────────────
  idx_pilot  <- sample.int(N, n0, replace = FALSE)
  pi_pilot   <- rep(n0 / N, N)
  beta_pilot <- fit_ipw(dat$X, dat$y, idx_pilot, pi = pi_pilot)
  if (is.null(beta_pilot)) beta_pilot <- beta

  # ── Plugin sigma^2 추정 (param / nonparam 분리) ────────────────────────────────
  sigma2_hat_param    <- estimate_sigma2_plugin_param(dat, idx_pilot, beta_pilot)
  # ── RNG 격리: RF가 내부에서 소비하는 난수가 이후 sampling에 영향 주지 않도록 ──
  rng_snapshot <- .Random.seed
  sigma2_hat_nonparam <- estimate_sigma2_plugin_nonparam(dat, idx_pilot, beta_pilot)
  .Random.seed <<- rng_snapshot

  oracle_score          <- sqrt(dat$sigma2_vec * dat$ell)
  plugin_score_param    <- sqrt(sigma2_hat_param * dat$ell)
  plugin_score_nonparam <- sqrt(sigma2_hat_nonparam * dat$ell)

  diag_sigma_cor_param    <- safe_cor(sigma2_hat_param, dat$sigma2_vec)
  diag_score_cor_param    <- safe_cor(plugin_score_param, oracle_score)
  diag_sigma_cor_nonparam <- safe_cor(sigma2_hat_nonparam, dat$sigma2_vec)
  diag_score_cor_nonparam <- safe_cor(plugin_score_nonparam, oracle_score)

  diag_sigma_cor_param_sp    <- safe_cor(sigma2_hat_param,    dat$sigma2_vec, method = "spearman")
  diag_sigma_cor_nonparam_sp <- safe_cor(sigma2_hat_nonparam, dat$sigma2_vec, method = "spearman")

  # ── Sampling scores ──────────────────────────────────────────────────────────
  scores <- get_scores(dat, sigma2_hat_param, sigma2_hat_nonparam)


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

    # Metric 1: leading term tr(A^-1 B_pi) -- 항상 TRUE sigma2_vec 사용, beta_hat 불필요
    er_leading <- mean(dat$sigma2_vec / pmax(pi_raw, 1e-10) * dat$ell)

    is_plugin_param    <- (method == "OPT-hetero-plugin-param")
    is_plugin_nonparam <- (method == "OPT-hetero-plugin-nonparam")

    diag_sigma_cor <- if (is_plugin_param) diag_sigma_cor_param
                       else if (is_plugin_nonparam) diag_sigma_cor_nonparam
                       else NA_real_
    diag_score_cor <- if (is_plugin_param) diag_score_cor_param
                       else if (is_plugin_nonparam) diag_score_cor_nonparam
                       else NA_real_
    diag_sigma_cor_spearman <- if (is_plugin_param) diag_sigma_cor_param_sp
                            else if (is_plugin_nonparam) diag_sigma_cor_nonparam_sp
                            else NA_real_

    if (is.null(beta_hat)) {
      return(list(
        method = method, fail = TRUE,
        er_leading = er_leading,
        pred_mse = NA_real_, pred_mse_N = NA_real_,
        param_mse = NA_real_, param_mse_N = NA_real_,
        c1_pi = NA_real_,
        n_realized = length(idx_final), n_clip = n_clip,
        ess_pi = ess_pi, max_weight = max_weight, cond_XtWX = cond_XtWX,
        diag_sigma_cor = diag_sigma_cor,
        diag_score_cor = diag_score_cor,
        diag_sigma_cor_spearman = diag_sigma_cor_spearman
      ))
    }

    # Metric 2: N*(MSE - MSE0), 실현된 test set (y_test) 기준 -- real data와 동일 로직
    mse_hat       <- mean(as.numeric(y_test - X_test %*% beta_hat)^2)
    pred_mse_raw  <- mse_hat - mse0
    param_mse_raw <- mean((beta_hat - beta)^2)

    list(
      method      = method, fail = FALSE,
      er_leading  = er_leading,                 # Metric 1 (이미 N-scale)
      pred_mse    = pred_mse_raw,
      pred_mse_N  = N * pred_mse_raw,            # Metric 2
      param_mse   = param_mse_raw,
      param_mse_N = N * param_mse_raw,           # Metric 3
      c1_pi       = mean(dat$sigma2_vec * dat$ell / pmax(pi_raw, 1e-10)),
      n_realized  = length(idx_final), n_clip = n_clip,
      ess_pi = ess_pi, max_weight = max_weight, cond_XtWX = cond_XtWX,
      diag_sigma_cor = diag_sigma_cor,
      diag_score_cor = diag_score_cor,
      diag_sigma_cor_spearman = diag_sigma_cor_spearman
    )
  }

  # ── FULL: 전체 데이터 OLS (pilot / subsampling 없음, er_leading 정의 안 됨) ────
  full_samp <- sample_full(dat)
  bhat_full <- full_samp$beta_hat
  full_result <- if (is.null(bhat_full)) {
    list(
      method = "FULL", fail = TRUE,
      er_leading = NA_real_,
      pred_mse = NA_real_, pred_mse_N = NA_real_,
      param_mse = NA_real_, param_mse_N = NA_real_,
      c1_pi = NA_real_,
      n_realized = full_samp$n_realized, n_clip = full_samp$n_clip,
      ess_pi = NA_real_, max_weight = NA_real_, cond_XtWX = NA_real_,
      diag_sigma_cor = NA_real_, diag_score_cor = NA_real_, diag_sigma_cor_spearman = NA_real_
    )
  } else {
    mse_hat_full       <- mean(as.numeric(y_test - X_test %*% bhat_full)^2)
    pred_mse_raw_full  <- mse_hat_full - mse0
    param_mse_raw_full <- mean((bhat_full - beta)^2)
    list(
      method      = "FULL", fail = FALSE,
      er_leading  = NA_real_,
      pred_mse    = pred_mse_raw_full,
      pred_mse_N  = N * pred_mse_raw_full,
      param_mse   = param_mse_raw_full,
      param_mse_N = N * param_mse_raw_full,
      c1_pi       = NA_real_,
      n_realized  = full_samp$n_realized, n_clip = full_samp$n_clip,
      ess_pi = NA_real_, max_weight = NA_real_, cond_XtWX = NA_real_,
      diag_sigma_cor = NA_real_, diag_score_cor = NA_real_, diag_sigma_cor_spearman = NA_real_
    )
  }

  results <- list()
  results[["FULL"]]                       <- full_result
  results[["SRS"]]                        <- run_method("SRS",                        scores[["SRS"]])
  results[["LEV-IPW"]]                    <- run_method("LEV-IPW",                    scores[["LEV-IPW"]])
  results[["OPT-homo"]]                   <- run_method("OPT-homo",                   scores[["OPT-homo"]])
  results[["OPT-hetero-oracle"]]          <- run_method("OPT-hetero-oracle",          scores[["OPT-hetero-oracle"]])
  results[["OPT-hetero-plugin-param"]]    <- run_method("OPT-hetero-plugin-param",    scores[["OPT-hetero-plugin-param"]])
  results[["OPT-hetero-plugin-nonparam"]] <- run_method("OPT-hetero-plugin-nonparam", scores[["OPT-hetero-plugin-nonparam"]])
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
      er_leading     = rv$er_leading,
      pred_mse       = rv$pred_mse,
      pred_mse_N     = rv$pred_mse_N,
      param_mse      = rv$param_mse,
      param_mse_N    = rv$param_mse_N,
      c1_pi          = rv$c1_pi,
      n_realized     = rv$n_realized,
      n_clip         = rv$n_clip,
      ess_pi         = rv$ess_pi,
      max_weight     = rv$max_weight,
      cond_XtWX      = rv$cond_XtWX,
      diag_sigma_cor = rv$diag_sigma_cor,
      diag_score_cor = rv$diag_score_cor,
      diag_sigma_cor_spearman = rv$diag_sigma_cor_spearman,
      stringsAsFactors = FALSE
    )
  }))
}